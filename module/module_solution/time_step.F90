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
    subroutine time_step(nbl)
    !======================================================================
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use metric_module
        use turbulent_module
        use cell_module
        use nodes_paras
        !>
        !>
		implicit none
		!>
		!>
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !>
        !>
        integer :: nbl ,&
                   i,j,k,&
                   idim,jdim,kdim
        real(kind = dprec) :: vterm,&
                              uban,&
                              vban,&
                              wban,&
                              a,&
                              t14,&
                              t24,&
                              t34,&
                              tin,&
                              b1,&
                              tvi,&
                              b2
        !>
        !>
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        !> computing the time step for time advanced
        !>
        !>
        vterm = max(1.3333,gamma/prandtl)*rmre

        do k = 2,kdim
            do j = 2,jdim
                do i = 2,idim
                    !>
                    !>
                    uban  =   iblock%cell(i,j,k)%u  *0.5e0*(iblock%metric(i,j,k)%fi  + iblock%metric(i+1,j,k)%fi) + &
                              iblock%cell(i,j,k)%v  *0.5e0*(iblock%metric(i,j,k)%fj  + iblock%metric(i+1,j,k)%fj) + &
                              iblock%cell(i,j,k)%w  *0.5e0*(iblock%metric(i,j,k)%fk  + iblock%metric(i+1,j,k)%fk)
                    !>
                    !>
                    vban  =   iblock%cell(i,j,k)%u  *0.5e0*(iblock%metric(i,j,k)%gi  + iblock%metric(i,j+1,k)%gi) + &
                              iblock%cell(i,j,k)%v  *0.5e0*(iblock%metric(i,j,k)%gj  + iblock%metric(i,j+1,k)%gj) + &
                              iblock%cell(i,j,k)%w  *0.5e0*(iblock%metric(i,j,k)%gk  + iblock%metric(i,j+1,k)%gk)
                    !>
                    !>
                    wban  =   iblock%cell(i,j,k)%u  *0.5e0*(iblock%metric(i,j,k)%hi  + iblock%metric(i,j,k+1)%hi) + &
                              iblock%cell(i,j,k)%v  *0.5e0*(iblock%metric(i,j,k)%hj  + iblock%metric(i,j,k+1)%hj) + &
                              iblock%cell(i,j,k)%w  *0.5e0*(iblock%metric(i,j,k)%hk  + iblock%metric(i,j,k+1)%hk)
                    !>
                    !>
                    a  =   sqrt(gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r)

                    !>
                    t14  = 0.5e0*(iblock%metric(i,j,k)%ff + iblock%metric(i+1,j,k)%ff)
                    t24  = 0.5e0*(iblock%metric(i,j,k)%gg + iblock%metric(i,j+1,k)%gg)
                    t34  = 0.5e0*(iblock%metric(i,j,k)%hh + iblock%metric(i,j,k+1)%hh)
                    !>
                    !>
                    tin =(abs(uban)+a)*t14+(abs(vban)+a)*t24+(abs(wban)+a)*t34
                    tin = tin*iblock%metric(i,j,k)%jacobian
                    !> computing the time step for viscous
                    !>
                    if(nvisc .eq. 1)then
                        !>
                        !>
                        b1   = 2.0*vterm*iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r
                        tvi  = b1*(t14*t14+t24*t24+t34*t34)*iblock%metric(i,j,k)%jacobian*iblock%metric(i,j,k)%jacobian
                        iblock%cell(i,j,k)%dt = 1.0/(tin+tvi)
                        !>
                        !>
                    elseif(nvisc .ge. 2)then
                        !>
                        !>
                        b1   = 2.0*vterm*(iblock%cell(i,j,k)%viscous + iblock%turbulent(i,j,k)%viscous)/iblock%cell(i,j,k)%r
                        tvi  = b1*(t14*t14+t24*t24+t34*t34)*iblock%metric(i,j,k)%jacobian*iblock%metric(i,j,k)%jacobian
                        iblock%cell(i,j,k)%dt = 1.0/(tin+tvi)
!                        write(501,'(3i4,5e24.16)') i-1,j-1,k-1,cfl*iblock%cell(i,j,k)%dt,iblock%metric(i,j,k)%volume,&
!                        b1*t14*iblock%metric(i,j,k)%jacobian,b1*t24*iblock%metric(i,j,k)%jacobian,b1*t34*iblock%metric(i,j,k)%jacobian
                        !>
                        !>
                    else
                        !>
                        !>
                        iblock%cell(i,j,k)%dt = 1.0/tin
                    end if
                    !>
                    !>
                    !>check the dt value not equit 0
                    !> nan will the dt time = 0
                    if(iblock%cell(i,j,k)%dt <= 1.e-15) then
                        open(unit=69,file='ccfd.error',status='unknown')
                        write(69,'("at block=",i4," i-",i4," j-",i4," k-",i4," dt=",e20.12," tin=",e20.12)') nbl,i,j,k,iblock%cell(i,j,k)%dt,tin
                        write(69,'(4e20.12)') t14,t24,t34,iblock%metric(i,j,k)%jacobian
                    end if

                end do
            end do
        end do
!        !>
        !open(112,file='step.txt',status='unknown')
        do  i = 2,idim
            do  k = 2,kdim
                do  j = 2,jdim
                    iblock%cell(i,j,k)%dt = cfl*iblock%cell(i,j,k)%dt
!                    if(nbl==93 .and. icyc==2) write(50093,'(1e24.16)') iblock%cell(i,j,k)%dt!,iblock%metric(i,j,k)%volume,iblock%metric(i,j,k)%jacobian
                end do
            end do
        end do
        do i = 2,idim
            do j = 2,jdim
                iblock%cell(i,j,kdim+1)%dt =iblock%cell(i,j,kdim)%dt
                iblock%cell(i,j,     1)%dt =iblock%cell(i,j, 2)%dt
            end do
        end do
        !>
        do j = 2,jdim
            do k = 2,kdim
                iblock%cell( 1    ,j,k)%dt = iblock%cell(   2,j,k)%dt
                iblock%cell(idim+1,j,k)%dt = iblock%cell(idim,j,k)%dt
            end do
        end do
        !>
        !>
        do i = 2,idim
            do k = 2,kdim
                iblock%cell(i,jdim+1,k)%dt =iblock%cell(i,jdim,k)%dt
                iblock%cell(i,     1,k)%dt =iblock%cell(i,   2,k)%dt
            end do
        end do

        return
    end subroutine time_step

