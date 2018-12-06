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
    subroutine vorticity_magnitude(nbl)
        !> purpose:evaluate vorticity magnitude  for use in
        !> detemining the turbulent eddy viscosity
        !> defined the cell_vor_tur array for store the vorticity magnitude
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
		implicit none
        !>
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        integer :: nbl,&
                   i,j,k,&
                   idim,jdim,kdim
        !> cross-flow temporaries
        real(kind =dprec),dimension(:,:,:),allocatable::cross_t
        !> components of vorticity
        real(kind =dprec),dimension(:,:,:,:),allocatable::vor_c
        !>
        !>
        real(kind =dprec):: term,&
                            factor
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        !>
        allocate(cross_t(1:idim+1,1:jdim+1,1:3))
        allocate(vor_c(1:idim+1,1:jdim+1,1:kdim+1,1:3))
        !>
        !>
        !>***************************************************************************************************
        !>
        !>
        !> calculate the vorticity magnitude
        !> at i-direction
        !> cycle through interface
        !>***************************************************************************************************
        do k=2,kdim
            do j=2,jdim
                do i=3,idim
                    term = iblock%metric(i,j,k)%ff / (iblock%metric(i,j,k)%volume + iblock%metric(i-1,j,k)%volume)
                    !>
                    cross_t(i,2,1) = term*((iblock%cell(i,j,k)%w       - iblock%cell(i-1,j,k)%w ) * iblock%metric(i,j,k)%fj - &
                                           (iblock%cell(i,j,k)%v       - iblock%cell(i-1,j,k)%v ) * iblock%metric(i,j,k)%fk)
                    !>
                    cross_t(i,2,2) = term*((iblock%cell(i,j,k)%u       - iblock%cell(i-1,j,k)%u ) * iblock%metric(i,j,k)%fk - &
                                           (iblock%cell(i,j,k)%w       - iblock%cell(i-1,j,k)%w ) * iblock%metric(i,j,k)%fi)
                    !>
                    cross_t(i,2,3) = term*((iblock%cell(i,j,k)%v       - iblock%cell(i-1,j,k)%v ) * iblock%metric(i,j,k)%fi - &
                                           (iblock%cell(i,j,k)%u       - iblock%cell(i-1,j,k)%u ) * iblock%metric(i,j,k)%fj)
                end do
                !>contribution i0-face
                !>
                !>
                i = 2
                term    =   iblock%metric(i,j,k)%ff / (iblock%metric(i,j,k)%volume + iblock%metric(i-1,j,k)%volume)
!                write(*,'("term=",3i4,4e20.12)') i,j,k,term,iblock%metric(i,j,k)%volume,iblock%metric(i-1,j,k)%volume,iblock%cell(i-1,j,k)%wall_blank
                factor  =   iblock%cell(i-1,j,k)%wall_blank
                term    =   factor*term
                !>
                cross_t(i,2,1) = term*((iblock%cell(i,j,k)%w       - iblock%cell(i-1,j,k)%w ) * iblock%metric(i,j,k)%fj - &
                                       (iblock%cell(i,j,k)%v       - iblock%cell(i-1,j,k)%v ) * iblock%metric(i,j,k)%fk)
                !>
                cross_t(i,2,2) = term*((iblock%cell(i,j,k)%u       - iblock%cell(i-1,j,k)%u ) * iblock%metric(i,j,k)%fk - &
                                       (iblock%cell(i,j,k)%w       - iblock%cell(i-1,j,k)%w ) * iblock%metric(i,j,k)%fi)
                !>
                cross_t(i,2,3) = term*((iblock%cell(i,j,k)%v       - iblock%cell(i-1,j,k)%v ) * iblock%metric(i,j,k)%fi - &
                                       (iblock%cell(i,j,k)%u       - iblock%cell(i-1,j,k)%u ) * iblock%metric(i,j,k)%fj)
                !>contribution idim-face
                !>
                !>
                i = idim+1
                term    =   iblock%metric(i,j,k)%ff / (iblock%metric(i,j,k)%volume + iblock%metric(i-1,j,k)%volume)
                factor  =   iblock%cell(i,j,k)%wall_blank
                term    =   factor*term
                !>
                cross_t(i,2,1) = term*((iblock%cell(i,j,k)%w       - iblock%cell(i-1,j,k)%w  ) * iblock%metric(i,j,k)%fj - &
                                       (iblock%cell(i,j,k)%v       - iblock%cell(i-1,j,k)%v  ) * iblock%metric(i,j,k)%fk)
                !>
                cross_t(i,2,2) = term*((iblock%cell(i,j,k)%u       - iblock%cell(i-1,j,k)%u  ) * iblock%metric(i,j,k)%fk - &
                                       (iblock%cell(i,j,k)%w       - iblock%cell(i-1,j,k)%w  ) * iblock%metric(i,j,k)%fi)
                !>
                cross_t(i,2,3) = term*((iblock%cell(i,j,k)%v       - iblock%cell(i-1,j,k)%v  ) * iblock%metric(i,j,k)%fi - &
                                       (iblock%cell(i,j,k)%u       - iblock%cell(i-1,j,k)%u  ) * iblock%metric(i,j,k)%fj)
                !>
                do i=2,idim
                    vor_c(i,j,k,1) = cross_t(i,2,1) + cross_t(i+1,2,1)
                    vor_c(i,j,k,2) = cross_t(i,2,2) + cross_t(i+1,2,2)
                    vor_c(i,j,k,3) = cross_t(i,2,3) + cross_t(i+1,2,3)
                end do
            end do
        end do

        !>***************************************************************************************************
        !>
        !>
        !> calculate the vorticity magnitude
        !> at j-direction
        !> cycle through interface
        !>***************************************************************************************************
        do k=2,kdim
            !>
            !>
            do j=3,jdim
                do i=2,idim
                    term = iblock%metric(i,j,k)%gg / (iblock%metric(i,j,k)%volume + iblock%metric(i,j-1,k)%volume)
                    !>
                    cross_t(i,j,1) = term*((iblock%cell(i,j,k)%w       - iblock%cell(i,j-1,k)%w  ) * iblock%metric(i,j,k)%gj - &
                                           (iblock%cell(i,j,k)%v       - iblock%cell(i,j-1,k)%v  ) * iblock%metric(i,j,k)%gk)
                        !>
                    cross_t(i,j,2) = term*((iblock%cell(i,j,k)%u       - iblock%cell(i,j-1,k)%u  ) * iblock%metric(i,j,k)%gk - &
                                           (iblock%cell(i,j,k)%w       - iblock%cell(i,j-1,k)%w  ) * iblock%metric(i,j,k)%gi)
                        !>
                    cross_t(i,j,3) = term*((iblock%cell(i,j,k)%v       - iblock%cell(i,j-1,k)%v  ) * iblock%metric(i,j,k)%gi - &
                                           (iblock%cell(i,j,k)%u       - iblock%cell(i,j-1,k)%u  ) * iblock%metric(i,j,k)%gj)
                end do
            end do
            !>contribution for j0-face
            !>
            !>
            j = 2
            !>
            do i=2,idim
                term = iblock%metric(i,j,k)%gg / (iblock%metric(i,j,k)%volume + iblock%metric(i,j-1,k)%volume)
                factor= iblock%cell(i,j-1,k)%wall_blank
    !            write(*,'("term=",3i4,4e20.12)') i,j,k,term,iblock%metric(i,j,k)%volume,iblock%metric(i,j-1,k)%volume,iblock%cell(i,j-1,k)%wall_blank
                term  = factor*term
                !>
                cross_t(i,j,1) = term*((iblock%cell(i,j,k)%w      - iblock%cell(i,j-1,k)%w  ) * iblock%metric(i,j,k)%gj - &
                                       (iblock%cell(i,j,k)%v      - iblock%cell(i,j-1,k)%v  ) * iblock%metric(i,j,k)%gk)
                        !>
                cross_t(i,j,2) = term*((iblock%cell(i,j,k)%u      - iblock%cell(i,j-1,k)%u  ) * iblock%metric(i,j,k)%gk - &
                                       (iblock%cell(i,j,k)%w      - iblock%cell(i,j-1,k)%w  ) * iblock%metric(i,j,k)%gi)
                        !>
                cross_t(i,j,3) = term*((iblock%cell(i,j,k)%v      - iblock%cell(i,j-1,k)%v  ) * iblock%metric(i,j,k)%gi - &
                                       (iblock%cell(i,j,k)%u      - iblock%cell(i,j-1,k)%u  ) * iblock%metric(i,j,k)%gj)
            end do
            !> contribution jdim-face
            !>
            !>
            j = jdim+1
            do i=2,idim
                !>
                !>
                term = iblock%metric(i,j,k)%gg / (iblock%metric(i,j,k)%volume + iblock%metric(i,j-1,k)%volume)
                factor= iblock%cell(i,j,k)%wall_blank
                term  = factor*term
                !>
                cross_t(i,j,1) = term*((iblock%cell(i,j,k)%w       - iblock%cell(i,j-1,k)%w ) * iblock%metric(i,j,k)%gj - &
                                       (iblock%cell(i,j,k)%v       - iblock%cell(i,j-1,k)%v  ) * iblock%metric(i,j,k)%gk)
                        !>
                cross_t(i,j,2) = term*((iblock%cell(i,j,k)%u       - iblock%cell(i,j-1,k)%u ) * iblock%metric(i,j,k)%gk - &
                                       (iblock%cell(i,j,k)%w       - iblock%cell(i,j-1,k)%w ) * iblock%metric(i,j,k)%gi)
                        !>
                cross_t(i,j,3) = term*((iblock%cell(i,j,k)%v       - iblock%cell(i,j-1,k)%v ) * iblock%metric(i,j,k)%gi - &
                                       (iblock%cell(i,j,k)%u       - iblock%cell(i,j-1,k)%u  ) * iblock%metric(i,j,k)%gj)
            end do
            !>
            !>
            do j=2,jdim
                do i=2,idim
                    vor_c(i,j,k,1) = vor_c(i,j,k,1) + cross_t(i,j,1) + cross_t(i,j+1,1)
                    vor_c(i,j,k,2) = vor_c(i,j,k,2) + cross_t(i,j,2) + cross_t(i,j+1,2)
                    vor_c(i,j,k,3) = vor_c(i,j,k,3) + cross_t(i,j,3) + cross_t(i,j+1,3)
                end do
            end do
            !>
            !>
        end do
        !>***************************************************************************************************
        !>
        !>
        !> calculate the vorticity magnitude
        !> at k-direction
        !> cycle through interface
        !>***************************************************************************************************
        !> for the three-dimension
        !>
        do k=3,kdim
            do j=2,jdim
                do i=2,idim
                    term = iblock%metric(i,j,k)%hh / (iblock%metric(i,j,k)%volume + iblock%metric(i,j,k-1)%volume)
                    !>
                    cross_t(i,j,1) = term*((iblock%cell(i,j,k)%w      - iblock%cell(i,j,k-1)%w  ) * iblock%metric(i,j,k)%hj - &
                                            (iblock%cell(i,j,k)%v      - iblock%cell(i,j,k-1)%v  ) * iblock%metric(i,j,k)%hk )
                    !>
                    cross_t(i,j,2) = term*((iblock%cell(i,j,k)%u      - iblock%cell(i,j,k-1)%u  ) * iblock%metric(i,j,k)%hk - &
                                            (iblock%cell(i,j,k)%w      - iblock%cell(i,j,k-1)%w  ) * iblock%metric(i,j,k)%hi)
                    !>
                    cross_t(i,j,3) = term*((iblock%cell(i,j,k)%v      - iblock%cell(i,j,k-1)%v  ) * iblock%metric(i,j,k)%hi - &
                                            (iblock%cell(i,j,k)%u      - iblock%cell(i,j,k-1)%u  ) * iblock%metric(i,j,k)%hj)
                    !>
                    !>
                    !>
                    !>
                    vor_c(i,j,k-1,1) = vor_c(i,j,k-1,1) + cross_t(i,j,1)
                    vor_c(i,j,k-1,2) = vor_c(i,j,k-1,2) + cross_t(i,j,2)
                    vor_c(i,j,k-1,3) = vor_c(i,j,k-1,3) + cross_t(i,j,3)
                    !>
                    !>
                    vor_c(i,j,k,1) = vor_c(i,j,k,1) + cross_t(i,j,1)
                    vor_c(i,j,k,2) = vor_c(i,j,k,2) + cross_t(i,j,2)
                    vor_c(i,j,k,3) = vor_c(i,j,k,3) + cross_t(i,j,3)
                end do
            end do
        end do
        !>contribution k0-face
        !>
        !>
        k = 2
        do j=2,jdim
            do i=2,idim
                !>
                !>
                term   = iblock%metric(i,j,k)%hh / (iblock%metric(i,j,k)%volume + iblock%metric(i,j,k-1)%volume)
                factor = iblock%cell(i,j,k-1)%wall_blank
    !            write(*,'("term=",3i4,4e20.12)') i,j,k,term,iblock%metric(i,j,k)%volume,iblock%metric(i,j,k-1)%volume,iblock%cell(i,j,k-1)%wall_blank
                term   = factor*term
                !>
                cross_t(i,j,1) = term*((iblock%cell(i,j,k)%w       - iblock%cell(i,j,k-1)%w  ) * iblock%metric(i,j,k)%hj - &
                                       (iblock%cell(i,j,k)%v       - iblock%cell(i,j,k-1)%v  ) * iblock%metric(i,j,k)%hk )
                !>
                cross_t(i,j,2) = term*((iblock%cell(i,j,k)%u       - iblock%cell(i,j,k-1)%u  ) * iblock%metric(i,j,k)%hk - &
                                       (iblock%cell(i,j,k)%w       - iblock%cell(i,j,k-1)%w  ) * iblock%metric(i,j,k)%hi)
                !>
                cross_t(i,j,3) = term*((iblock%cell(i,j,k)%v       - iblock%cell(i,j,k-1)%v  ) * iblock%metric(i,j,k)%hi - &
                                       (iblock%cell(i,j,k)%u       - iblock%cell(i,j,k-1)%u  ) * iblock%metric(i,j,k)%hj)
                !>
                !>
                vor_c(i,j,2,1) = vor_c(i,j,2,1) + cross_t(i,j,1)
                vor_c(i,j,2,2) = vor_c(i,j,2,2) + cross_t(i,j,2)
                vor_c(i,j,2,3) = vor_c(i,j,2,3) + cross_t(i,j,3)
            end do
        end do
        !>contribution kdim-face
        !>
        !>
        k = kdim+1
        do j=2,jdim
            do i=2,idim
                !>
                !>
                term   = iblock%metric(i,j,k)%hh / (iblock%metric(i,j,k)%volume + iblock%metric(i,j,k-1)%volume)
    !            write(*,'("term=",3i4,4e20.12)') i,j,k,term,iblock%metric(i,j,k)%volume,iblock%metric(i,j,k-1)%volume,iblock%cell(i,j,k)%wall_blank
                factor = iblock%cell(i,j,k)%wall_blank
                term   = factor*term
                !>
                cross_t(i,j,1) = term*((iblock%cell(i,j,k)%w       - iblock%cell(i,j,k-1)%w  ) * iblock%metric(i,j,k)%hj - &
                                       (iblock%cell(i,j,k)%v       - iblock%cell(i,j,k-1)%v  ) * iblock%metric(i,j,k)%hk )
                !>
                cross_t(i,j,2) = term*((iblock%cell(i,j,k)%u       - iblock%cell(i,j,k-1)%u  ) * iblock%metric(i,j,k)%hk - &
                                       (iblock%cell(i,j,k)%w       - iblock%cell(i,j,k-1)%w  ) * iblock%metric(i,j,k)%hi)
                !>
                cross_t(i,j,3) = term*((iblock%cell(i,j,k)%v       - iblock%cell(i,j,k-1)%v  ) * iblock%metric(i,j,k)%hi - &
                                       (iblock%cell(i,j,k)%u       - iblock%cell(i,j,k-1)%u  ) * iblock%metric(i,j,k)%hj)
                !>
                !>
                vor_c(i,j,kdim,1) = vor_c(i,j,kdim,1) + cross_t(i,j,1)
                vor_c(i,j,kdim,2) = vor_c(i,j,kdim,2) + cross_t(i,j,2)
                vor_c(i,j,kdim,3) = vor_c(i,j,kdim,3) + cross_t(i,j,3)
            end do
        end do
        !>
        !>
        !>calculate the vorticity magnitude
        !>
        !>
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    iblock%turbulent(i,j,k)%vorticity = sqrt(vor_c(i,j,k,1)*vor_c(i,j,k,1) + vor_c(i,j,k,2)*vor_c(i,j,k,2) + vor_c(i,j,k,3)*vor_c(i,j,k,3))
                end do
            end do
        end do
        !>
        !>
        deallocate(cross_t)
        deallocate(vor_c)

        return
        !>
    end subroutine vorticity_magnitude
