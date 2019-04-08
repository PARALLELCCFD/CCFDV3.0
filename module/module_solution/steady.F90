!   _______________________________________________________________________________
!    _______/\\\\\\\\\_______/\\\\\\\\\______/\\\\\\\\\\\\\\__/\\\\\\\\\\\\\\\______
!     _____/\\\////////_____/\\\////////_____\/\\\//////////__\/\/////////////________
!      ___/\\\/____________/\\\/______________\/\\\____________\/\\\________/\\\\_______
!       __/\\\____________/\\\_________________\/\\\\\\\\\\\\\__\/\\\_______\////________
!        _\/\\\___________\/\\\_________________\/\\\/////////___\/\\\________/\\\\_______
!         _\//\\\__________\//\\\________________\/\\\____________\/\\\_______\////________
!          __\///\\\_________\///\\\______________\/\\\____________\/\\\______/\\\\_________
!           ____\////\\\\\\\\\____\////\\\\\\\\\___\/\\\____________\/\\\\\\\\\////__________
!            _______\/////////_______\/////////_____\///_____________\/////////_______________
!            __________________________________________________________________________________
!----------------------------------------------------------------------------------------------
!>  subroutine input_files_cgns
!>  last edit 2016-05-24
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
    subroutine steady()
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
        use bc_module
        use metric_module
        use variable_module
        use nodes_paras
        use nodes_var
        use nodes_mg
        use nodes_var_bc
		!>
		!>
		implicit none

#if defined PMPI
        include "mpif.h"
        integer istat(mpi_status_size)
#endif
        !>
        !>
        !>
        integer :: iprint
        integer :: n_time
        integer :: ntt,nit
        integer :: level
        integer :: isize_level,size_b
        integer :: num_grid_coarse
        integer :: idim,jdim,kdim,&
                   i,j,k,&
                   ii,&
                   irange,jrange,krange,&
                   istore,jstore,kstore
        integer :: nbl,ierr,ip,mytag,nd_src,i_size1to1_out
        !>
        !>
        real(kind = dprec)::error_l2,error,&
                            epson,&
                            r_tmp,&
                            p_tmp,error_l2_out
        real(kind = dprec),dimension(:),allocatable::ii_send_host
        real(kind = dprec):: start_t,&
                             finish_t,&
                             resmax
        !>
        !>
        type(blocks_type),pointer  :: iblock
        type(bc_types),pointer     :: ibc
        type(overlap_type),pointer :: mesh
        mesh => grids(imesh)
        !>
        !>
        !limiter = 0 !=0,no limiter,=1,minmod limiter
        ntt    = 0
        epson  = 1.e-10
        iprint = 0
        error  = 0.0
        i_size1to1_out  = 0

        if(restart .eq. 0) icyc_restart = 0
        !>
        !>
        !> set the multigrid max level
        !> if not used mg methods the global_level will eq 1
        if(meshseqflags .ne. 0 )then
            level  = imeshseque
        else
            level  = global_level
        end if
#if defined PMPI
        !>
        !>
        !> allocate the array for 1-1 blocks's boundary conditions
        !> and set the mg methods array for the leg of downward and upward
        !>
        call bc_cut_array_allocate()
        !>
        !>
        call mpi_allreduce(i_size1to1,i_size1to1_out,1,mpi_integer,mpi_max,mycomm,ierr)
        i_size1to1 = i_size1to1_out
#endif
        !>
        call initial
        !>
        !>
        !> begin the iterations for computing
iterations:do icyc=1+icyc_restart,ncycle(imeshseque)+icyc_restart
            !>
            !>
            error_l2 = 0.0
            error_l2_out    = 0.0
            !>
            !>
#if defined PMPI
            if(myid .eq. myhost)then
#endif
                call cpu_time(start_t)
#if defined PMPI
            end if
#endif
            !>
            !>
            nit = 0
            !>
            !>
            if( mgflags .ne. 0 )then
                nit = 2
            else
                nit = 1
            end if
            if(meshseqflags .ne. 0 )then
                if(imeshseque .eq. 1)then
                    nit = 1
                end if
            end if
            !>
            !>
            !> kode = 1 at finest mesh
            !> kode = 2 at no-finest mesh
            kode = 1
            !>
            !>
            !> mode = 0 at finest mesh
            !> mode = 1 at no-finest mesh
            mode = 0
            999  continue
            !> when mgflags=1,then nit=2
            !> it means the roe method will call twice
            !> first roe is computing roe(i(n,n-1)q(n))
            !> second roe is computing roe(q(n-1)) this
            !>
            do n_time =1,nit
                !>
                !>
                !> update the boundary conditions
                !> shared the conservative variable at the 1-1 boundary condition
                !> after the restrict,and the variable will be used to computing
                !> the residual.
                !> send the information to the neibour blocks's faces
                !> and the information about r,ru,rv,rw,re
                !>
                do nbl =1,nblocks
                    if(level_mg(nbl) .eq. level)then
                        call boundary(nbl,level)
                    end if
                end do
#if defined PMPI
                !>
                !> the 1-1 bc for parallel
                !> used the
                if(myid .ne. myhost)then
                    call bc_1to1_parallel(level)
                end if
                !>
                !>
#endif
                !>
                !>
                !>the cycle for blocks
                do nbl =1,nblocks

#if defined PMPI
                    if(myid .eq. n2p(nbl) .and. level_mg(nbl) .eq. level)then
#endif
                        !>
                        !>
                        !>
                        if(n_time .eq. 1)then
                            if(meshseqflags .ne. 0)then
                                if(level .ne. imeshseque)then
                                    !> kode = 1 at finest mesh
                                    !> kode = 2 at no-finest mesh
                                    kode = 2
                                    !> mode = 0 at finest mesh
                                    !> mode = 1 at no-finest mesh
                                    mode = 1
                                endif
                            else
                                if(level .ne. global_level)then
                                    !> kode = 1 at finest mesh
                                    !> kode = 2 at no-finest mesh
                                    kode = 2
                                    !> mode = 0 at finest mesh
                                    !> mode = 1 at no-finest mesh
                                    mode = 1
                                endif
                            end if
                        endif
                        !>
                        !> begin the downward leg of multigrid cycle
                        !>
                        !> from the finer mesh to coarser mesh
                        !>
                        if(level_mg(nbl) .ne. level ) goto 10001
                        !>
                        !>
                        iblock => mesh%blocks(nbl)
                        idim = iblock%idim
                        jdim = iblock%jdim
                        kdim = iblock%kdim
                        !>
                        !> initial the rhs term and the viscous term
                        do k =1,kdim+1
                            do j=1,jdim+1
                                do i=1,idim+1
                                    !>
                                    !>
                                    !>
                                    iblock%variable(i,j,k)%res_1 = 0.0
                                    iblock%variable(i,j,k)%res_2 = 0.0
                                    iblock%variable(i,j,k)%res_3 = 0.0
                                    iblock%variable(i,j,k)%res_4 = 0.0
                                    iblock%variable(i,j,k)%res_5 = 0.0
                                    if(ivisc_i .ge. 1 .or. ivisc_j .ge. 1 .or. ivisc_k .ge. 1)then
                                        iblock%cell(i,j,k)%viscous          = 0.0
                                    end if
                                    !>
                                    !>
                                    !>
                                end do
                            end do
                        end do
                        !
                        !> computing the roe u v w(tree-dimension)
                        !> from the rn run rvn rwn ren
                        !>
                        !>
                        !>
                        do k = 0, iblock%kdim + 2
                            do j = 0, iblock%jdim + 2
                                do i = 0, iblock%idim + 2
                                    !> at the wall bondary conditions
                                    !> the value is gradient so it maybe
                                    !> a negative value
                                    if(iblock%cell(i,j,k)%wall_blank .lt. 0)then
                                        !>
                                        if(iblock%cell(i,j,k)%r <= epson .or. iblock%cell(i,j,k)%p <= epson) then
                                        !>
                                            if(iblock%cell(i,j,k)%r < epson)then
#if defined DEBUG
                                                open(unit=69,file='ccfd.error',status='unknown')
                                                write(69,'("ntime=",i1," roe at block=",i3," iter=",i4," x=",i3," y=",i3,&
                                                " z=",i3," r=",e20.10)') n_time,nbl,icyc,i,j,k,iblock%cell(i,j,k)%r
#endif
                                                iblock%cell(i,j,k)%r  = epson
                                            end if
                                            if( iblock%cell(i,j,k)%p  .lt. epson)then
#if defined DEBUG
                                                open(unit=69,file='ccfd.error',status='unknown')
                                                write(69,'("ntime=",i1,"   p at block=",i3," iter=",i4," x=",i3," y=",i3,&
                                                " z=",i3," p=",e20.10)') n_time,nbl,icyc,i,j,k,iblock%cell(i,j,k)%p
#endif
                                                iblock%cell(i,j,k)%p  =  epson
                                            end if
                                            !>
                                        end if
                                    end if
                                        !>
                                end do
                            end do
                        end do
                        !>
                        !> compute the laminar flow viscous
                        !> used the sutherland's methods
                        !>
                        !>
                        call laminar_viscous_coefficient(nbl)
                        !>
                        !> compute the time step delta_t
                        !> each icy of iteration must re-computing
                        !> the delta t of time for computing
                        !>
                        if(n_time .eq. 1)  call time_step(nbl)
                        !>
                        !>
                        !> calculate viscosities coefficient
                        !> when nvisc not equal to  0 and the filed have viscous
                        !>
                        if(nvisc .gt. 1 .and. n_time .eq. 1) then
                            call turbulent_viscous_coefficient(nbl,level)
                        end if
                        !>
                        !>
                        !> compute the inviscid flux used the split-flux methods like roe scheme etc.
                        !> use the roe's scheme or Van Leer's scheme
                        !>
                        !>
                        !>
                        call muscl(nbl)
                        !>
                        !>
                        !> compute the viscous flux
                        !>
                        !>
                        if(ivisc_i .ge. 1 .or. ivisc_j .ge. 1 .or. ivisc_k .ge. 1 ) then
                            call flux_viscous(nbl)
                        end if
                        !>
                        !>
                        !>
                        !>
                        if(meshseqflags .ne. 0)then
                            if(level_mg(nbl) .eq. imeshseque .and. n_time .eq. 1)then
                                !>
                                !>
                                !> use L2-norm residual
                                !>
                                error  = 0.0
                                do k=2,iblock%kdim
                                    do j=2,iblock%jdim
                                        do i=2,iblock%idim
                                            error   = error+iblock%variable(i,j,k)%res_1**2
                                        end do
                                    end do
                                end do
                                !>
                                !>
                                if( error .ne. error .or. (real(error) .lt. 0) )then
                                    open(unit=69,file='ccfd.error',status='unknown')
                                    write(69,'("NaN detected after residual evaluation, at block",i6," cycle=",i8)') nbl,icyc
                                    stop
                                end if
                                error_l2 = error_l2 + error
                            end if
                        else
                            if(level_mg(nbl) .eq. global_level .and. n_time .eq. 1)then
                                !>
                                !>
                                !> use L2-norm residual
                                !>
                                error  = 0.0
                                resmax = 0.0
                                do k=2,iblock%kdim
                                    do j=2,iblock%jdim
                                        do i=2,iblock%idim
                                            error   = error + iblock%variable(i,j,k)%res_1**2
                                        end do
                                    end do
                                end do
                                !>
                                !>
                                if( error .ne. error  .or. (real(error) .lt. 0) )then
                                    open(unit=69,file='ccfd.error',status='unknown')
                                    write(69,'("NaN detected after residual evaluation, at block",i6," cycle=",i8)') nbl,icyc
                                    stop
                                end if
                                error_l2 = error_l2 + error
                                !>
                                !>
                            end if
                        end if
                        !
                        !> the subroutine  tau_term
                        !>
                        !> computing the tau term for correct value of n+1 level mesh
                        !> add correction from finer levels
                        !>
                        if(meshseqflags .ne. 0 )then
                            if( mgflags .ne. 0  .and. level .lt. imeshseque)then
                                call tau_term(nbl,level)
                            end if
                        else
                            if( mgflags .ne. 0  .and. level .lt. global_level)then
                                call tau_term(nbl,level)
                            end if

                        end if
                        !>the subroutine pre_smoothting
                        !>
                        !> call the implicit residual smoothing method to smooting the resiual
                        !> if desired by the value irs =1
                        if(irs .ne. 0  )then
                            call pre_smoothing(nbl)
                        end if
                        !>
                        !>
                        !>  at the first ntime cycle,will computing the q
                        !>  this value is q(n)
                        if(n_time .eq. 1)then
                            !> output the residual infomation at the first ntime
                            !> and the reisual used ls methods or 2-norm
                            !>
                            !> implicit time-solution ( lu-adi,diag algorithms)
                            !> the time schemes will add the krylov subspace method such as GMRES,CG etc
                            call lu_adi(nbl)
                            !>
                            !> computing the time solution for chimea mesh
                            !> the method from ccfdv2.0
                            !> and update by *
                            if(overlap .eq. 1)then
                                !> update the overlapped values if chimea scheme is used
                                !> null
                                !> null
                                !> null
                            end if
                            !> modify the solutions
                            !>
                            !> update solution
                            !>
                            !>
                            !>
                            do k = 2,iblock%kdim
                                do j = 2,iblock%jdim
                                    do i = 2,iblock%idim
                                        r_tmp  =  iblock%cell(i,j,k)%r + iblock%variable(i,j,k)%res_1*iblock%cell(i,j,k)%blank
                                        p_tmp  =  iblock%cell(i,j,k)%p + iblock%variable(i,j,k)%res_5*iblock%cell(i,j,k)%blank
                                        !>
                                        !>
                                        call modify_primitive(r_tmp,p_tmp,iblock%cell(i,j,k)%r,iblock%cell(i,j,k)%p)
                                        !>
                                        !>
                                        iblock%cell(i,j,k)%u = iblock%cell(i,j,k)%u + iblock%variable(i,j,k)%res_2*iblock%cell(i,j,k)%blank
                                        iblock%cell(i,j,k)%v = iblock%cell(i,j,k)%v + iblock%variable(i,j,k)%res_3*iblock%cell(i,j,k)%blank
                                        iblock%cell(i,j,k)%w = iblock%cell(i,j,k)%w + iblock%variable(i,j,k)%res_4*iblock%cell(i,j,k)%blank
                                        !>
                                        !>
                                    end do
                                end do
                            end do
                            !>
                            !>
                        endif
                        !>
                        !> restrict the residual to coarser mesh from finer mesh
                        !> used the volume weight methods
                        !> and the coarerest mesh have not restrict operator
                        !>
                        if(level .gt. 1 .and. n_time .eq. nit)then
                            call restrict(nbl)
                        end if
                        !>
                        !>
                        !>
                        10001   continue
#if defined PMPI
                    endif !blocks to processor
#endif
                end do ! nblocks cycle
                !>
                !>
                !> at the coarser mesh level=1，not used restrict the residual
                !> so it can done here and upward coarser mesh to finer mesh
#if defined PMPI
                if(myid .ne. myhost)then
#endif
                    if(level .eq. 1) goto 1000
                    if(meshseqflags .ne. 0)then
                        if(level .ne. imeshseque)then
                            kode = 1
                            mode = 1
                        end if
                    else
                        if(level .ne. global_level)then
                            kode = 1
                            mode = 1
                        end if
                    end if
#if defined PMPI
                end if
#endif
            end do !n_time cycle body
            !>
            !>
            1000   continue
            !>
            !>
#if defined PMPI
            if(myid .ne. myhost)then
#endif
                !>
                !>
                !>
                level = level - 1
                !> change the kode and mode
                !> kode = 2
                !> mode = 1
                if(level .ge. 1) goto 999
                !>
                !>  end the downward log of multigrid methods
                !>
                !>end do multigrid_level
                !>
                !> begin or continue upward leg of multigrid ccyle
                !> at upward of mg methods the solution will be correct by delta_e
                !> and the prolonged method is area rule
                !>
                !>
                level = level + 1
#if defined PMPI
            end if
#endif
            !>
            !>
            if(mgflags .ne. 0)then
                1999 continue

                do nbl=1,nblocks
#if defined PMPI
                    if(myid .eq. n2p(nbl))then
#endif
                        if(level_mg(nbl) .ne. level) goto 2000
                        if(meshseqflags .ne. 0)then
                            if(level .eq. imeshseque) goto 2001
                        else
                            if(level .eq. global_level) goto 2001
                        end if
                        !> the subroutine prolonged
                        !>
                        !> used the area rule prolongation the correct value from coarser mesh
                        !> to finer meshs and correct the value at levle finer meshs
                        !> q(n)= q(n)+i(n-1,n)(q(n-1) - i(n,n-1)q(n)) at n level
                        !>
                        call prolonged(nbl)
                        !> share the correcteed delta_q information with the processores
                        !> if used the overlap mesh
                        if(overlap .eq. 1)then
                            !> if no-used the overlap mesh
                            !> it reserve for overlap mesh
                            !> null
                            !> null
                            !> null
                        endif
                        !> end the upward and downward for one-w-cycle
                        !>
                        if(vcycle .ne. 0 .and. wcycle .ne. 0 )then
                            !> have a error for multigrid cycle
                            !> will set the cycle of mg for default
                            write(*,*) 'have a error ,v/w-cycle cannot used in the same time !'
                            write(*,*) 'ccfdv3.0 set the default cycle of mg is v-cycle    !  '
                            wcycle = 0
                        end if
                        if(vcycle .eq. 0 .and. wcycle .eq. 0)then
                            write(*,*) 'have a error, and you set the v/w-cycle  all zero !'
                            write(*,*) 'ccfdv3.0 set the default cycle of mg is v-cycle !  '
                            vcycle = 1
                        end if
                        !> send and recv the boundary condition information from neighbor blocks
                        !> synchronization before the next cycle computing
                        2000 continue
#if defined PMPI
                    end if
#endif
                end do !end the blocks cycle
                !>
                !>
                !>
#if defined PMPI
                if(myid .ne. myhost)then
#endif
                    level = level + 1
                    !>
                    !>
                    !> goto for one-w cylce of multigrid methods
                    if(vcycle .eq. 0)then
                        if(wcycle .eq. 1  .and. level .eq. 2 .and. ntt .lt. 2)then
                             ntt = ntt + 1
                             goto 999
                        end if
                        !> goto for double-w-cycle of multigrid methods
                        !>
                        if(wcycle .eq. 2 .and. level .eq. 2 .and. ntt .lt. 2)then
                            ntt = ntt + 1
                            goto 999
                        end if
                        !>
                    end if

                    if(meshseqflags .ne. 0)then
                        if(level .lt. imeshseque) goto 1999
                    else
                        if(level .lt. global_level) goto 1999
                    end if
#if defined PMPI
                end if
#endif
                !>
                !> end upward leg of multi-grid cycle
                !>
            end if !end the mg methods
            !>
            !>
            2001 continue
            !>
            !>
            iprint = iprint + 1
#if defined PMPI
            call mpi_barrier(mycomm,ierr)
#endif
            !>
            !>
            !>output the force coefficient
            call force
            !>
            !>
            !>reduce the residual from the processor to the main host
            !>
            !>
#if defined PMPI
            call mpi_allreduce(error_l2,error_l2_out,1,mpi_double_precision,mpi_sum,mycomm,ierr)
            error_l2 = error_l2_out
            !>
            !>
            !>**********send r,u,v,w,p to host  processor!***************
            !> when the flow or the solution must be store or output by the host processor
            !> send the flow and the solution to the host processor
            !> used the sen/recv mpi api
            !>
            if(iprint .eq. ntres .or. (icyc .eq. ncycle(meshseque) .and. imeshseque .eq. meshseque) ) then
                !>
                !>
                do nbl=1,nblocks,global_level
                    !>
                    !>
                    iblock => mesh%blocks(nbl)
                    nd_src = n2p(nbl)
                    mytag  = nbl
                    size_b = (iblock%idim+1)*(iblock%jdim+1)*(iblock%kdim+1)*5
                    !>
                    !>
                    if(myid .eq. nd_src .or. myid .eq. myhost)then
                        allocate(ii_send_host(size_b))
                    endif
                    !>
                    !>
                    !>
                    if(myid .eq. nd_src)then
                        ip = 0
                        do k=1,iblock%kdim+1
                            do j=1,iblock%jdim+1
                                do i=1,iblock%idim+1
                                    ip               = ip + 1
                                    ii_send_host(ip) = iblock%cell(i,j,k)%r
                                    ip               = ip + 1
                                    ii_send_host(ip) = iblock%cell(i,j,k)%u
                                    ip               = ip + 1
                                    ii_send_host(ip) = iblock%cell(i,j,k)%v
                                    ip               = ip + 1
                                    ii_send_host(ip) = iblock%cell(i,j,k)%w
                                    ip               = ip + 1
                                    ii_send_host(ip) = iblock%cell(i,j,k)%p
                                end do
                            end do
                        enddo
                        call mpi_send(ii_send_host,size_b,mpi_double_precision,myhost,mytag,mycomm,istat,ierr)
                    end if
                    !>
                    !>
                    !>
                    if(myid .eq. myhost)then
                        call mpi_recv(ii_send_host,size_b,mpi_double_precision,nd_src,mytag,mycomm,istat,ierr)
                        ip = 0
                        do k=1,iblock%kdim+1
                            do j=1,iblock%jdim+1
                                do i=1,iblock%idim+1
                                    ip                    = ip + 1
                                    iblock%cell(i,j,k)%r  = ii_send_host(ip)
                                    if(iblock%cell(i,j,k)%r < 0.)then
                                        print *,'r <0',i,j,k,nbl,iblock%cell(i,j,k)%r
                                    endif
                                    ip                    = ip + 1
                                    iblock%cell(i,j,k)%u  = ii_send_host(ip)
                                    ip                    = ip + 1
                                    iblock%cell(i,j,k)%v  = ii_send_host(ip)
                                    ip                    = ip + 1
                                    iblock%cell(i,j,k)%w  = ii_send_host(ip)
                                    ip                    = ip + 1
                                    iblock%cell(i,j,k)%p  = ii_send_host(ip)
                                end do
                            end do
                        enddo
                    end if
                    !>
                    !>
                    !>
                    if(myid == myhost .or. myid == nd_src)then
                        deallocate(ii_send_host)
                    endif
                    !>
                    !>
                end do !>do nbl=1,nblocks
                !>
                !>
            end if
            !**********send r,u,v,w,p to host processor!***************
#endif

#if defined PMPI

            if(myid == myhost)then
#endif


                if(meshseqflags .ne. 0)then
                    num_grid_coarse = num_grid / (8**(meshseque - imeshseque))
                    error_l2   = sqrt(error_l2/float(num_grid_coarse))
                else
                    error_l2   = sqrt(error_l2/float(num_grid))
                end if
                !>
                !> the finish time to this iterations
                !>
                call cpu_time(finish_t)
                !>
                !>
                finish_t = finish_t - start_t
                !>
                call residual_show(error_l2,finish_t)
                !>
                !> output the variables
                if(iprint .eq. ntres ) then
                    call restart_output()
                end if
#if defined PMPI
            endif
            !>
            if(iprint .eq. ntres ) iprint =0
            call mpi_barrier(mycomm,ierr)
            !>
            !>
#endif
        end do iterations
        !>
        !> mesh sequence for FMG methods
        !> not the end of mg
        !>
        !>
#if defined PMPI
        if(myid == myhost .and. imeshseque .eq. meshseque)then
            !>
            !> output the flow solutions and restart files
            call output
            !>
            !>
        endif
        !>
        call mpi_barrier(mycomm,ierr)
        call deallocate_nodes_var_bc()
#endif

        return
    end subroutine steady


