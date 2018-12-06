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
    subroutine input_grids_cgns
        !>
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use nodes_var
        use nodes_mg
        use nodes_paras
        use nodes_var_bc
        !>
        !>
        !>
        implicit none
        include "cgnslib_f.h"
        !>
        !>
        !>
#if defined PMPI
        include "mpif.h"
#endif
        !>
        !>
        !***********************************************************************
        !
        !     function: groups those boundaries requiring non-local information
        !               this needs to be done when computing in parallel where
        !               a boundary may be treated by a few of the processes
        !               this is e.g. required to integrate a quantity along
        !               the boundary
        !
        !     inputs:    name    type           description
        !
        !     1          pglob   pointer        global pointer to root
        !     2          blocks    pointer        blocks pointer
        !
        !     outputs :  name    type           description
        !
        !     called from: cgns grid
        !chime
        !     creation by: liuxz
        !
        !     creation date: 2016-04-21
        !
        !----------------------------------------------------------------------
        !>
        !>
        integer :: nm,&
                   idimc,jdimc,kdimc,&
                   immax,jmmax,kmmax,&
                   xyz_buf,nd_dest,mytag
        integer :: ii,jj,kk
        integer :: nbl
        integer :: i,j,k
        integer :: iii,jjj,kkk
        integer :: ierr,iccg,ier
        !>
        !>
        !>
        type(blocks_type),pointer :: iblock
        type(blocks_type),pointer :: iiblock
        type(overlap_type),pointer:: mesh
        !>
        !>
        !>
#if defined PMPI
        integer::istat(mpi_status_size)
        real(kind=dprec),dimension(:),allocatable::nobl_send_xyz
#endif
        !>
        !>
        !>
        !>
        allocate(grids(meshes))
        !>
        !>
        !>
        !>
        call get_grid_num()
        !>
        !>
        !>
        !>
        !>
        mesh => grids(imesh)
        if(myid .ne. myhost)then
            !>
            !>
            allocate(mesh%blocks(nblocks))
        end if
        !>
        !>read each blocks from grid file
        !>
        if(myid == myhost)then
            write(*,*) '*********************************************'
            write(*,*) '      ccfd.3.0:  read grid files ......    '
            write(*,*) '*********************************************'
        end if
        !>
        !>
        !>
        !>
        !>
#if defined PMPI
        !>
        !> read the coordinate of the grid by sequence
        !> when have a large-scale case want to parallel computing
        !>
        if(myid .eq. myhost)then
#endif
            !>
            !>
            !> open the cgns grid file
            open (11,file =grid_files,status='old')
            call cg_open_f(grid_files,mode_read,iccg,ier)
            if (ier .ne. 0) call cg_error_exit_f
            !>
            !>
            nbl = 1
            do nm=1,nobl
                if(mgflags .ne. 0 )then
                    !>
                    !>
                    call get_blocks_coor(nbl,iccg)
                    !>
                    !>
                    level_mg(nbl) = global_level
                    !>
                    !>
                    write(61,'("nbl=",i3," and level_mg=",i3)') nbl,level_mg(nbl)
                    write(*,'("nbl=",i3," and level_mg=",i3)') nbl,level_mg(nbl)
                    !>
                    !>
                    do ncg = 1,global_level-1
                        nbl = nbl + 1
                        kk  = 1
                        !>
                        !>
                        iiblock => mesh%blocks(nbl-1)
                        iblock  => mesh%blocks(nbl)
                        !>
                        do k=2,iiblock%kdim+1,2
                            kk = kk + 1
                            jj = 1
                            do j =2,iiblock%jdim+1,2
                                jj = jj + 1
                                ii = 1
                                do i =2,iiblock%idim+1,2
                                    ii                              =  ii + 1
                                    iblock%coordinate(ii,jj,kk)%x   =  iiblock%coordinate(i,j,k)%x
                                    iblock%coordinate(ii,jj,kk)%y   =  iiblock%coordinate(i,j,k)%y
                                    iblock%coordinate(ii,jj,kk)%z   =  iiblock%coordinate(i,j,k)%z
                                    !>
                                    !>
                                end do
                            end do
                        end do
                        !>
                        !>
                        level_mg(nbl) = global_level - ncg
                        write(61,'("block_id=",i6," idim=",i4," jdim=",i4," kdim=",i4," and level_mg=",i3)') nbl,ii,jj,kk,level_mg(nbl)
                        write(*,'("block_id=",i6," idim=",i4," jdim=",i4," kdim=",i4," and level_mg=",i3)') nbl,ii,jj,kk,level_mg(nbl)
                    end do
                    nbl = nbl + 1
                else
                    !>
                    !> read the not multi-grid grid just only 1-level
                    !>
                    call get_blocks_coor(nm,iccg)
                    level_mg(nm) = 1
                    write(61,'("block_id=",i6," and level_mg=",i6)') nm,level_mg(nm)
                    write(*,'("block_id=",i6," and level_mg=",i6)') nm,level_mg(nm)
                endif
            end do
            !>
            !> close the cgns grid files
            call cg_close_f(iccg,ier)
            if (ier .ne. 0) call cg_error_exit_f
#if defined PMPI
        end if !only the myhost
        !>
        !>
        call mpi_bcast(level_mg,nblocks,mpi_integer,myhost,mycomm,ierr)
        !>
        !>
#endif
        !> set the multi-block's boundary's information
        !>
        !>
        call get_cgns_blocks_boundary
        !>
        !>
        !>
        !> mapping the blocks to processor
        !> can use DD methods or call api form metis
        !> mapping the blocks to processor by grids of each blocks
        !> it's a load balance methods
#if defined PMPI
        if(myid .eq. myhost)then
#endif
            !>
            !>
            !> load balance methods for mapping the blocks to processors and
            !> ensure that blocks are evenly distributed to the processor
            !>
            !>
            call load_balance
            !>
            !>
#if defined PMPI
        end if
#endif
        !>
        !>
#if defined PMPI
        !>
        !>
        call mpi_bcast(n2p,nblocks,mpi_integer,myhost,mycomm,ierr)
        !>
        !>
        !>
#endif
        !>
        !>
#if defined DEBUG
        if(myid == 0)then
            open(unit=61,file='ccfd.out',form='formatted',status='unknown')
            write(61,'("-------------------------------------------------")')
            write(61,'("*the information of blocks mapping to processor**")')
            write(61,'("---------------------ccfdv3.0------------------------")')
            do nbl=1,nblocks
                write(61,'("block_id=",i4,"  to processor=",i4)') nbl,n2p(nbl)
            end do
            write(61,'("*************************************************")')
        endif
#endif
        !>
        !>
        !> check boundary conditions before computing
        !>
        call boundary_check
        !>
        !>
        !>
#if defined PMPI
        !>
        !>
        nbl = 0
        if(myid .ne. myhost)then
            if(mgflags .ne. 0)then
                do nm =1,nblocks,global_level
                    if(myid .eq. n2p(nm))then
                        nbl = nm - ((nm-1)/global_level)*(global_level-1)
                        idimc = nxt(nbl)
                        jdimc = nyt(nbl)
                        kdimc = nzt(nbl)
                        write(*,'("myid=",i3," blocks=",i3," idimc=",i3," jdimc=",i3," kdimc=",i3)') myid,nm,idimc,jdimc,kdimc
                        call mem_allocate_block(mesh%blocks(nm),idimc,jdimc,kdimc)
                        do ncg=1,global_level-1
                            idimc = idimc/2 + 1
                            jdimc = jdimc/2 + 1
                            kdimc = kdimc/2 + 1
                            call mem_allocate_block(mesh%blocks(nm+ncg),idimc,jdimc,kdimc)
                            write(*,'("myid=",i3," blocks=",i3," idimc=",i3," jdimc=",i3," kdimc=",i3)') myid,nm+ncg,idimc,jdimc,kdimc
                        end do
                    end if
                end do
            else
                do nm=1,nblocks,global_level
                    if(myid .eq. n2p(nm))then
                        idimc = nxt(nm)
                        jdimc = nyt(nm)
                        kdimc = nzt(nm)
                        call mem_allocate_block(mesh%blocks(nm),idimc,jdimc,kdimc)
                    end if
                end do
            end if
        end if
        !>
        !>
        !>
        immax = 0
        jmmax = 0
        kmmax = 0
        nbl   = 0
        !>
        !>
        do nm=1,nblocks,global_level
            nbl = nm - (nm/global_level)*(global_level-1)
            immax = max(immax,nxt(nbl))
            jmmax = max(jmmax,nyt(nbl))
            kmmax = max(kmmax,nzt(nbl))
        end do
        xyz_buf = 3*immax*jmmax*kmmax
        allocate(nobl_send_xyz(xyz_buf))
        call mpi_barrier(mycomm,ierr)
        !>
        !> send the level of multi-grid to other processors from host nodes
        !> send the coordinate to other processors
        do nbl=1,nblocks
            !>
            !>
            if(myid == myhost)then
                !>
                !>
                iblock => mesh%blocks(nbl)
                xyz_buf = 3*iblock%idim*iblock%jdim*iblock%kdim
                !>
                !>
                nd_dest = n2p(nbl)
                mytag   = nbl
                ii = 1
                do kkk =2,iblock%kdim+1
                    do jjj=2,iblock%jdim+1
                        do iii =2,iblock%idim+1
                            nobl_send_xyz(ii) = iblock%coordinate(iii,jjj,kkk)%x
                            ii = ii + 1
                            nobl_send_xyz(ii) = iblock%coordinate(iii,jjj,kkk)%y
                            ii = ii + 1
                            nobl_send_xyz(ii) = iblock%coordinate(iii,jjj,kkk)%z
                            ii = ii + 1

                        end do
                    end do
                end do
                if(xyz_buf .ne. ii-1)then
                    write(*,*) 'send the coordinate have some error,check the buffer size  '
                    stop
                end if
                call mpi_send(nobl_send_xyz,xyz_buf,mpi_double_precision,nd_dest,mytag,mycomm,ierr)
                !>
                !>
                !>
            endif
            !>recv form host nodes
            !>
            !>
            if(myid == n2p(nbl))then
                iblock => mesh%blocks(nbl)
                xyz_buf = 3*iblock%idim*iblock%jdim*iblock%kdim
                mytag = nbl
                jj = 1
                call mpi_recv(nobl_send_xyz,xyz_buf,mpi_double_precision,myhost,mytag,mycomm,istat,ierr)
                do kkk=2,iblock%kdim+1
                    do jjj = 2,iblock%jdim+1
                        do iii = 2,iblock%idim+1
                            iblock%coordinate(iii,jjj,kkk)%x = nobl_send_xyz(jj)
                            jj = jj + 1
                            iblock%coordinate(iii,jjj,kkk)%y = nobl_send_xyz(jj)
                            jj = jj + 1
                            iblock%coordinate(iii,jjj,kkk)%z = nobl_send_xyz(jj)
                            jj = jj + 1
                        end do
                    end do
                end do
                if(xyz_buf .ne. jj-1)then
                    write(*,*) 'recv the coordinate have some error,check the buffer size  '
                    stop
                end if
            endif
!
        end do
        call mpi_barrier(mycomm,ierr)
        deallocate(nobl_send_xyz)

#endif
        !>
        !>
#if defined PMPI
        call mpi_barrier(mycomm,ierr)
#endif

        !> compute the free stream conditions
        !>
        !>
        !> initial the flow condition
        call initial
        !>
        !>
        !> initial the cell variable
        !>
        !>
        do nbl = 1, nblocks
#if defined PMPI
            if(myid == n2p(nbl) .or. myid == myhost)then
#endif
                iblock => mesh%blocks(nbl)
                do k  = 0,iblock%kdim + 2
                    do j  = 0,iblock%jdim + 2
                        do i  = 0,iblock%idim + 2
                            !>
                            !>
                            !>
                            iblock%cell(i,j,k)%r           = r0
                            iblock%cell(i,j,k)%u           = u0
                            iblock%cell(i,j,k)%v           = v0
                            iblock%cell(i,j,k)%w           = w0
                            iblock%cell(i,j,k)%p           = p0
                            !>
                            !> when the cell at the no-slip wall bc the wall_blank = 2.0
                            !>
                            iblock%cell(i,j,k)%blank       = 1.0
                            iblock%cell(i,j,k)%wall_blank  = 1.0
                        end do
                    end do
                end do
#if defined PMPI
            endif
#endif
        end do
        !>
        !>
        !>define the viscous variable
        !>
        !>
        if(ivisc_i .ge. 1 .or. ivisc_j .ge. 1 .or. ivisc_k .ge. 1)then
            do nbl = 1, nblocks
#if defined PMPI
                if(myid == n2p(nbl) )then
#endif
                    iblock => mesh%blocks(nbl)
                    do k  = 0,iblock%kdim + 2
                        do j  = 0,iblock%jdim + 2
                            do i  = 0,iblock%idim + 2
                            !>
                            !>
                            iblock%cell(i,j,k)%viscous  = 0.0
                            !>
                            !>
                            end do
                        end do
                    end do
#if defined PMPI
                endif
#endif
            end do
        end if
        !>
        !>
        !> define the turbulent variable
        !>
        !>
        if(ivisc_i .ge. 2 .or. ivisc_j .ge. 2 .or. ivisc_k .ge. 2)then
            do nbl = 1, nblocks
#if defined PMPI
                if(myid == n2p(nbl))then
#endif
                    iblock => mesh%blocks(nbl)
                    do k  = 0,iblock%kdim + 2
                        do j  = 0,iblock%jdim + 2
                            do i  = 0,iblock%idim + 2
                            !>
                            !>
                            iblock%turbulent(i,j,k)%viscous     = 0.0
                            iblock%turbulent(i,j,k)%vorticity   = 0.0
                            !>
                            !>
                            end do
                        end do
                    end do
#if defined PMPI
                endif
#endif
            end do
        end if
        !>
        !>
        !**********************************************
    end subroutine input_grids_cgns
