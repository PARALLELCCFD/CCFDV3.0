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
    !> \brief
    !> purpose:interpolate the solution or the correction from
    !> a coarser meshs to a finer meshs
    !> and the methods from handbooks
    !> amd we used the area law method to prolong operator
    !>********************************************************
    !> the prolong also can used the straight injection
    !> we will update this method later!
    subroutine prolonged(nbl)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        !>
        !>
        implicit none
        !>

        integer :: idim,jdim,kdim
        integer :: idim_c,jdim_c,kdim_c
        integer :: idim_c1,jdim_c1,kdim_c1
        integer :: nbl,nbl1,iss,npoint
        integer :: i,j,k,n
        integer :: ii,jj,kk
        integer :: ii1,jj1,kk1
        type(blocks_type),pointer  :: iblock
        type(blocks_type),pointer  :: iblock_c
        type(overlap_type),pointer :: mesh
        !>
        !>
        !>
        real(kind = dprec) alpq,phiq,betq,tt1,tt2,logical_depending,qq
        !>
        !>
        real(kind = dprec),dimension(:,:,:,:),allocatable  :: delta_q
        real(kind = dprec),dimension(:,:,:,:),allocatable  :: tmp_q1
        real(kind = dprec),dimension(:,:,:,:),allocatable  :: tmp_q2
        !> idim,jdim,kdim finer meshs indices
        !> idim_c,jdim_c,kdim_c coarser meshs indices
        !> coarser blocks num nbl
        !> so finer meshs blocks num is nbl1=nbl-1
        !> interpolate solution to finer mesh (mode =0)
        !> interpolate correction to finer mesh (mode = 1)
        !> note mode=1 must only be used for rho,u,v,w,p(primitive variable)
        !> i(n-1,n)delta_q(n-1) to n mesh
        !> q(n) = q(n) + i(n-1,n)(q(n-1)-i(n,n-1)q(n))
        !> and used the delta_q to correct the q value of n meshes
        !>  by liuxz

        nbl1     =  nbl - 1
        mesh     => grids(imesh)
        iblock   => mesh%blocks(nbl1)
        iblock_c => mesh%blocks(nbl)
        !>
        !>
        !> the coordinate of the finer meshes
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        !>
        !>
        !> the coordiante of the coarser meshes
        idim_c = iblock_c%idim
        jdim_c = iblock_c%jdim
        kdim_c = iblock_c%kdim
!        write(*,'("nbl=",i3," nbl1=",i3,"  fine axis=",3i4," coarse axis=",3i3)') nbl,nbl1,idim,jdim,kdim,idim_c,jdim_c,kdim_c

        !>
        !>
        !>
        !> allocate the array for delta q
        allocate(delta_q(2:idim,2:jdim,2:kdim,1:7))
        allocate(tmp_q1(2:idim_c,2:jdim_c,2:kdim_c,1:7))
        allocate(tmp_q2(2:idim,2:jdim,2:kdim,1:7))

        !> inint the dela_q
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    do n=1,5
                      delta_q(i,j,k,n) = 0.e0
                    end do
                end do
            end do
        end do

        !> mode = 0 interpolate solution to finer mesh
        if(mode .eq. 0)then
            do k=2,kdim_c
                do j=2,jdim_c
                    do i=2,idim_c
                        tmp_q1(i,j,k,1) =  iblock_c%cell(i,j,k)%r
                        tmp_q1(i,j,k,2) =  iblock_c%cell(i,j,k)%u
                        tmp_q1(i,j,k,3) =  iblock_c%cell(i,j,k)%v
                        tmp_q1(i,j,k,4) =  iblock_c%cell(i,j,k)%w
                        tmp_q1(i,j,k,5) =  iblock_c%cell(i,j,k)%p
                        if(nvisc .ge. 3)then
                            !> used the full multi-grid method ,when the level at the
                            !> top of mesh sequence level ,need prolonged
                            !> the turbulent from coarse meshed to finer meshes
                            !>
                            !>
                            tmp_q1(i,j,k,6) = iblock_c%turbulent(i,j,k)%tur_save
                            tmp_q1(i,j,k,7) = iblock_c%turbulent(i,j,k)%viscous
                        end if
                    enddo
                enddo
            enddo
        !> mode =1 interpolate correction to finer mesh
        !> delta_q = i
        else !mode = 1
            do k=2,kdim_c
                do j=2,jdim_c
                    do i=2,idim_c
                        tmp_q1(i,j,k,1) = iblock_c%cell(i,j,k)%r - iblock_c%variable_mg(i,j,k)%r_mg
                        tmp_q1(i,j,k,2) = iblock_c%cell(i,j,k)%u - iblock_c%variable_mg(i,j,k)%u_mg
                        tmp_q1(i,j,k,3) = iblock_c%cell(i,j,k)%v - iblock_c%variable_mg(i,j,k)%v_mg
                        tmp_q1(i,j,k,4) = iblock_c%cell(i,j,k)%w - iblock_c%variable_mg(i,j,k)%w_mg
                        tmp_q1(i,j,k,5) = iblock_c%cell(i,j,k)%p - iblock_c%variable_mg(i,j,k)%p_mg
                        !>
                        !>
                    end do
                enddo
            enddo
            !> use the irs methods to smooth the rhs
            !> smooth corrections
            if(ics .ne. 0)then
               call post_smoothing(idim_c,jdim_c,kdim_c,tmp_q1)
            end if
            !>
            !>
        endif
        !> interpolate in j-direction
        !> jdim_c1 = jdim_c - 1
        !> kdim_c1 = kdim_c - 1
        !>
        do i=2,idim_c
            do k=2,kdim_c
                !>
                tmp_q2(i   ,2,k,1)    = tmp_q1(i,2     ,k,1)
                tmp_q2(i,jdim,k,1)    = tmp_q1(i,jdim_c,k,1)
                tmp_q2(i,   2,k,2)    = tmp_q1(i,     2,k,2)
                tmp_q2(i,jdim,k,2)    = tmp_q1(i,jdim_c,k,2)
                tmp_q2(i,   2,k,3)    = tmp_q1(i,2     ,k,3)
                tmp_q2(i,jdim,k,3)    = tmp_q1(i,jdim_c,k,3)
                tmp_q2(i,   2,k,4)    = tmp_q1(i,     2,k,4)
                tmp_q2(i,jdim,k,4)    = tmp_q1(i,jdim_c,k,4)
                tmp_q2(i,   2,k,5)    = tmp_q1(i,     2,k,5)
                tmp_q2(i,jdim,k,5)    = tmp_q1(i,jdim_c,k,5)
                if(nvisc .ge. 3 .and. mode ==0 )then
                    tmp_q2(i,   2,k,6)    = tmp_q1(i,     2,k,6)
                    tmp_q2(i,jdim,k,6)    = tmp_q1(i,jdim_c,k,6)
                    tmp_q2(i,   2,k,7)    = tmp_q1(i,     2,k,7)
                    tmp_q2(i,jdim,k,7)    = tmp_q1(i,jdim_c,k,7)
                end if
                jj = 2
                do j=2,jdim_c-1
                    jj = jj + 1
                    !> rho
                    tmp_q2(i,jj,k,1) = 0.75e0*tmp_q1(i,j,k,1) + 0.25e0*tmp_q1(i,j+1,k,1)
                    !> ru
                    tmp_q2(i,jj,k,2) = 0.75e0*tmp_q1(i,j,k,2) + 0.25e0*tmp_q1(i,j+1,k,2)
                    !> rv
                    tmp_q2(i,jj,k,3) = 0.75e0*tmp_q1(i,j,k,3) + 0.25e0*tmp_q1(i,j+1,k,3)
                    !> rw
                    tmp_q2(i,jj,k,4) = 0.75e0*tmp_q1(i,j,k,4) + 0.25e0*tmp_q1(i,j+1,k,4)
                    !> p
                    tmp_q2(i,jj,k,5) = 0.75e0*tmp_q1(i,j,k,5) + 0.25e0*tmp_q1(i,j+1,k,5)
                    !>
                    if(nvisc .ge. 3 .and. mode ==0)then
                        !>
                        !> tur_save
                        tmp_q2(i,jj,k,6) = 0.75e0*tmp_q1(i,j,k,6) + 0.25e0*tmp_q1(i,j+1,k,7)
                        !> turbulent viscous
                        tmp_q2(i,jj,k,7) = 0.75e0*tmp_q1(i,j,k,6) + 0.25e0*tmp_q1(i,j+1,k,7)
                    end if
                    !>
                    jj = jj + 1
                    !> rho
                    tmp_q2(i,jj,k,1) = 0.25e0*tmp_q1(i,j,k,1) + 0.75e0*tmp_q1(i,j+1,k,1)
                    !> ru
                    tmp_q2(i,jj,k,2) = 0.25e0*tmp_q1(i,j,k,2) + 0.75e0*tmp_q1(i,j+1,k,2)
                    !> rv
                    tmp_q2(i,jj,k,3) = 0.25e0*tmp_q1(i,j,k,3) + 0.75e0*tmp_q1(i,j+1,k,3)
                    !> rw
                    tmp_q2(i,jj,k,4) = 0.25e0*tmp_q1(i,j,k,4) + 0.75e0*tmp_q1(i,j+1,k,4)
                    !> p
                    tmp_q2(i,jj,k,5) = 0.25e0*tmp_q1(i,j,k,5) + 0.75e0*tmp_q1(i,j+1,k,5)
                    !>
                    if(nvisc .ge. 3 .and. mode ==0)then
                        !> turbulent save values
                        tmp_q2(i,jj,k,6) = 0.25e0*tmp_q1(i,j,k,6) + 0.75e0*tmp_q1(i,j+1,k,6)
                        !> turbulent viscous
                        tmp_q2(i,jj,k,7) = 0.25e0*tmp_q1(i,j,k,7) + 0.75e0*tmp_q1(i,j+1,k,7)
                    end if
                end do
            end do
        end do

        !> interpolate in k-direction
        !>
        ii = idim - idim_c + 1
        do i=2,idim_c
            ii = ii + 1
            do j=2,jdim
                !> delta rho
                delta_q(ii,j,2,1) = tmp_q2(i,j,2,1)
                !> delta u
                delta_q(ii,j,2,2) = tmp_q2(i,j,2,2)
                !> delta v
                delta_q(ii,j,2,3) = tmp_q2(i,j,2,3)
                !> delta w
                delta_q(ii,j,2,4) = tmp_q2(i,j,2,4)
                !> delta p
                delta_q(ii,j,2,5) = tmp_q2(i,j,2,5)
                !>
                !>
                if(nvisc .ge. 3 .and. mode ==0)then
                    delta_q(ii,j,2,6) = tmp_q2(i,j,2,6)
                    delta_q(ii,j,2,7) = tmp_q2(i,j,2,7)
                end if
            end do
            !>
            do jj1=2,jdim
                !> delta rho
                delta_q(ii,jj1,kdim,1) = tmp_q2(i,jj1,kdim_c,1)
                !> delta u
                delta_q(ii,jj1,kdim,2) = tmp_q2(i,jj1,kdim_c,2)
                !> delta v
                delta_q(ii,jj1,kdim,3) = tmp_q2(i,jj1,kdim_c,3)
                !> delta w
                delta_q(ii,jj1,kdim,4) = tmp_q2(i,jj1,kdim_c,4)
                !> delta p
                delta_q(ii,jj1,kdim,5) = tmp_q2(i,jj1,kdim_c,5)
                !>
                if(nvisc .ge. 3 .and. mode ==0)then
                    delta_q(ii,jj1,kdim,6) = tmp_q2(i,jj1,kdim_c,6)
                    delta_q(ii,jj1,kdim,7) = tmp_q2(i,jj1,kdim_c,7)
                end if
            end do
            !>
            kk = 2
            do k=2,kdim_c-1
                kk = kk + 1
                do jj1=2,jdim
                    !> delta rho
                    delta_q(ii,jj1,kk,1) = 0.75e0*tmp_q2(i,jj1,k,1) + 0.25e0*tmp_q2(i,jj1,k+1,1)
                    !> delta u
                    delta_q(ii,jj1,kk,2) = 0.75e0*tmp_q2(i,jj1,k,2) + 0.25e0*tmp_q2(i,jj1,k+1,2)
                    !> delta v
                    delta_q(ii,jj1,kk,3) = 0.75e0*tmp_q2(i,jj1,k,3) + 0.25e0*tmp_q2(i,jj1,k+1,3)
                    !> delta w
                    delta_q(ii,jj1,kk,4) = 0.75e0*tmp_q2(i,jj1,k,4) + 0.25e0*tmp_q2(i,jj1,k+1,4)
                    !> delta re
                    delta_q(ii,jj1,kk,5) = 0.75e0*tmp_q2(i,jj1,k,5) + 0.25e0*tmp_q2(i,jj1,k+1,5)
                    !>
                    !>
                    if(nvisc .ge. 3 .and. mode ==0)then
                        delta_q(ii,jj1,kk,6) = 0.75e0*tmp_q2(i,jj1,k,6) + 0.25e0*tmp_q2(i,jj1,k+1,6)
                        delta_q(ii,jj1,kk,7) = 0.75e0*tmp_q2(i,jj1,k,7) + 0.25e0*tmp_q2(i,jj1,k+1,7)
                    end if
                end do
                kk = kk + 1
                do jj1=2,jdim
                    !> delat rho
                    delta_q(ii,jj1,kk,1) = 0.25e0*tmp_q2(i,jj1,k,1) + 0.75e0*tmp_q2(i,jj1,k+1,1)
                    !> delat u
                    delta_q(ii,jj1,kk,2) = 0.25e0*tmp_q2(i,jj1,k,2) + 0.75e0*tmp_q2(i,jj1,k+1,2)
                    !> delat v
                    delta_q(ii,jj1,kk,3) = 0.25e0*tmp_q2(i,jj1,k,3) + 0.75e0*tmp_q2(i,jj1,k+1,3)
                    !> delat w
                    delta_q(ii,jj1,kk,4) = 0.25e0*tmp_q2(i,jj1,k,4) + 0.75e0*tmp_q2(i,jj1,k+1,4)
                    !> delat p
                    delta_q(ii,jj1,kk,5) = 0.25e0*tmp_q2(i,jj1,k,5) + 0.75e0*tmp_q2(i,jj1,k+1,5)
                    !>
                    !>
                    if(nvisc .ge. 3 .and. mode ==0)then
                        delta_q(ii,jj1,kk,6) = 0.25e0*tmp_q2(i,jj1,k,6) + 0.75e0*tmp_q2(i,jj1,k+1,6)
                        delta_q(ii,jj1,kk,7) = 0.25e0*tmp_q2(i,jj1,k,7) + 0.75e0*tmp_q2(i,jj1,k+1,7)
                    end if
                end do
            end do
            !>
        end do
        !> interpolate in i-direction
        !> three-dimension meshs
        if(idim .gt. 2)then
            iss    = idim - idim_c + 1 + 1
            do j=2,jdim
                do k=2,kdim
                    !> rho
                    delta_q(2,j,k,1) = delta_q(iss,j,k,1)
                    !> u
                    delta_q(2,j,k,2) = delta_q(iss,j,k,2)
                    !> v
                    delta_q(2,j,k,3) = delta_q(iss,j,k,3)
                    !> w
                    delta_q(2,j,k,4) = delta_q(iss,j,k,4)
                    !> p
                    delta_q(2,j,k,5) = delta_q(iss,j,k,5)
                    !>
                    if(nvisc .ge. 3 .and. mode ==0)then
                        delta_q(2,j,k,6) = delta_q(iss,j,k,6)
                        delta_q(2,j,k,7) = delta_q(iss,j,k,7)
                    end if
                    !> rho
                    delta_q(idim,j,k,1) = delta_q(idim,j,k,1)
                    !> u
                    delta_q(idim,j,k,2) = delta_q(idim,j,k,2)
                    !> v
                    delta_q(idim,j,k,3) = delta_q(idim,j,k,3)
                    !> w
                    delta_q(idim,j,k,4) = delta_q(idim,j,k,4)
                    !> p
                    delta_q(idim,j,k,5) = delta_q(idim,j,k,5)
                    !>
                    !>
                    if(nvisc .ge. 3 .and. mode ==0)then
                        delta_q(idim,j,k,6) = delta_q(idim,j,k,6)
                        delta_q(idim,j,k,7) = delta_q(idim,j,k,7)
                    end if
                end do
            end do
            !>
            ii = 2
            do i=iss,idim-1
                ii = ii + 1
                !>
                do j=2,jdim
                    do k=2,kdim
                        !> rho
                        delta_q(ii,j,k,1) = 0.75e0*delta_q(i,j,k,1) + 0.25e0*delta_q(i+1,j,k,1)
                        !> u
                        delta_q(ii,j,k,2) = 0.75e0*delta_q(i,j,k,2) + 0.25e0*delta_q(i+1,j,k,2)
                        !> v
                        delta_q(ii,j,k,3) = 0.75e0*delta_q(i,j,k,3) + 0.25e0*delta_q(i+1,j,k,3)
                        !> w
                        delta_q(ii,j,k,4) = 0.75e0*delta_q(i,j,k,4) + 0.25e0*delta_q(i+1,j,k,4)
                        !> p
                        delta_q(ii,j,k,5) = 0.75e0*delta_q(i,j,k,5) + 0.25e0*delta_q(i+1,j,k,5)
                        !>
                        !>
                        if(nvisc .ge. 3 .and. mode ==0)then
                            delta_q(ii,j,k,6) = 0.75e0*delta_q(i,j,k,6) + 0.25e0*delta_q(i+1,j,k,6)
                            delta_q(ii,j,k,7) = 0.75e0*delta_q(i,j,k,7) + 0.25e0*delta_q(i+1,j,k,7)
                        end if
                    end do
                end do
                !>
                ii = ii + 1
                !> j,k-faces
                do j=2,jdim
                    do k=2,kdim
                        !> rho
                        delta_q(ii,j,k,1) = 0.25e0*delta_q(i,j,k,1) + 0.75e0*delta_q(i+1,j,k,1)
                        !> u
                        delta_q(ii,j,k,2) = 0.25e0*delta_q(i,j,k,2) + 0.75e0*delta_q(i+1,j,k,2)
                        !> v
                        delta_q(ii,j,k,3) = 0.25e0*delta_q(i,j,k,3) + 0.75e0*delta_q(i+1,j,k,3)
                        !> w
                        delta_q(ii,j,k,4) = 0.25e0*delta_q(i,j,k,4) + 0.75e0*delta_q(i+1,j,k,4)
                        !> p
                        delta_q(ii,j,k,5) = 0.25e0*delta_q(i,j,k,5) + 0.75e0*delta_q(i+1,j,k,5)
                        !>
                        !>
                        if(nvisc .ge. 3 .and. mode ==0)then
                            delta_q(ii,j,k,6) = 0.25e0*delta_q(i,j,k,6) + 0.75e0*delta_q(i+1,j,k,6)
                            delta_q(ii,j,k,7) = 0.25e0*delta_q(i,j,k,7) + 0.75e0*delta_q(i+1,j,k,7)
                        end if
                    end do
                end do
            end do
        end if
        !>
        !> interpolate solution to finer mesh (mode =0)
        !>
        if(mode .eq. 0)then
            do k=2,kdim
                do j=2,jdim
                    do i=2,idim
                        !> update rho
                        iblock%cell(i,j,k)%r  = delta_q(i,j,k,1)
                        !> update ru
                        iblock%cell(i,j,k)%u  = delta_q(i,j,k,2)
                        !> update rv
                        iblock%cell(i,j,k)%v  = delta_q(i,j,k,3)
                        !> update rw
                        iblock%cell(i,j,k)%w  = delta_q(i,j,k,4)
                        !> update the pressure
                        iblock%cell(i,j,k)%p  = delta_q(i,j,k,5)
                        !>
                        !>
                        if(nvisc .ge. 3)then
                            iblock%turbulent(i,j,k)%tur_save = delta_q(i,j,k,6)
                            iblock%turbulent(i,j,k)%viscous  = delta_q(i,j,k,7)
                        end if
                    end do
                end do
            end do
        !>
        !> interpolate correction to finer mesh (mode = 1)
        else
            !> update density and pressure to ensure positivity
            !> "cut-off" point is determined by alpq
            !> minimum value of density is equal to 1/phiq
            alpq = -0.2
            phiq = 1.0/0.5
            betq = 1.0 + alpq*phiq
            do k=2,kdim
                do j=2,jdim
                    do i=2,idim
                        tt1 = delta_q(i,j,k,1)/iblock%cell(i,j,k)%r
                        tt2 = delta_q(i,j,k,1)/(betq+abs(tt1)*phiq)
                        !> output the density  is negative number
                        delta_q(i,j,k,1) = logical_depending(tt2,delta_q(i,j,k,1),(real(tt1) .lt. real(alpq)))
                        !>
#if defined DEBUG
                        if(real(tt1) .lt. real(alpq))then
                            open(unit=69,file='ccfd.multigrid',status='unknown')
                            write(69,'("At the prolonged sub,and iter=",i4," at block[",i5,"]@x-",i4," y-",i4," z-",i4,&
                            " density is negative number and the values= ",e20.12," and rho=",e20.12)') icyc,nbl1,&
                            i,j,k,delta_q(i,j,k,1),iblock%cell(i,j,k)%r
                        end if
#endif
                        tt1 = delta_q(i,j,k,5)/iblock%cell(i,j,k)%p
                        tt2 = delta_q(i,j,k,5)/(betq+abs(tt1)*phiq)
                        !> output the  pressure is negative number
                        delta_q(i,j,k,5) = logical_depending(tt2,delta_q(i,j,k,5),(real(tt1) .lt. real(alpq)))
#if defined DEBUG
                        if(real(tt1) .lt. real(alpq))then
                            open(unit=69,file='ccfd.multigrid',status='unknown')
                                                        write(69,'("At the prolonged sub,and iter=",i4," at block[",i5,"]@x-",i4," y-",i4," z-",i4,&
                            " pressure is negative number and the values= ",e20.12," and rho=",e20.12)') icyc,nbl1,&
                            i,j,k,delta_q(i,j,k,5),iblock%cell(i,j,k)%p
                        end if
#endif
                    end do
                end do
            end do
            !> update the primitive variables
            !> and used the blanks of each cell
            do k=2,kdim
                do j=2,jdim
                    do i=2,idim
                        !> update the rho variables
                        iblock%cell(i,j,k)%r   = iblock%cell(i,j,k)%r  + iblock%cell(i,j,k)%blank*delta_q(i,j,k,1)
                        !> update the u-direction variables
                        iblock%cell(i,j,k)%u   = iblock%cell(i,j,k)%u  + iblock%cell(i,j,k)%blank*delta_q(i,j,k,2)
                        !> update the v-direction variables
                        iblock%cell(i,j,k)%v   = iblock%cell(i,j,k)%v  + iblock%cell(i,j,k)%blank*delta_q(i,j,k,3)
                        !> update the w-direction variables
                        iblock%cell(i,j,k)%w   = iblock%cell(i,j,k)%w  + iblock%cell(i,j,k)%blank*delta_q(i,j,k,4)
                        !> update the pressure variables
                        iblock%cell(i,j,k)%p   = iblock%cell(i,j,k)%p  + iblock%cell(i,j,k)%blank*delta_q(i,j,k,5)
                        if(iblock%cell(i,j,k)%p .lt. 0.0)then
                            write(*,'("have negative values when prolonged values from coarse mesh to finer mesh, cycle=",i6," block=",i6," axis=",3i4," p=",e24.16)'),icyc,nbl,i,j,k,iblock%cell(i,j,k)%p
                        end if
                    end do
                end do
            end do
            !> end the mode of select
        endif
        !> end the prolonged operating and deallocate the array of used!
        deallocate(delta_q)
        deallocate(tmp_q1)
        deallocate(tmp_q2)
        return

    end subroutine prolonged

