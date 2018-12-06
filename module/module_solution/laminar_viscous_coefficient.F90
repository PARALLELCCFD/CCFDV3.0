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
    subroutine laminar_viscous_coefficient(nbl)
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use nodes_paras
        !>
        !>
		implicit none
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        real(kind = dprec) :: t,suth
        integer :: nbl
        integer :: i,j,k
        !>
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        !> sutherland's constant
        suth   = 198.6/tinf
        !>
        !> laminar flow
        !> calculate the laminar viscous coefficient
        !>
        do k = 1,iblock%kdim + 1
            do j = 1,iblock%jdim + 1
                do i = 1,iblock%idim + 1
                    t            = gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r
                    iblock%cell(i,j,k)%viscous = t*sqrt(t)*(1.0+suth)/(t+suth)
!                    write(180801,'(4i4,6e24.16)') icyc,i,j,k,iblock%cell(i,j,k)%r,iblock%cell(i,j,k)%u,iblock%cell(i,j,k)%v,iblock%cell(i,j,k)%w,iblock%cell(i,j,k)%p,iblock%cell(i,j,k)%viscous
                end do
            end do
        end do
        !>
        !>
        return
    end subroutine laminar_viscous_coefficient


