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
    subroutine turbulent_viscous_coefficient(nbl,level)
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
        integer :: nbl,level
        !>
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        !> turbulent viscosity - finest grids only
        !>
        !>
        if((.not. meshseqflags .and. level .ge. global_level) .or. (meshseqflags .and. level .ge. imeshseque))then
            !>
            !>
            if(ivisc_i .eq. 2 .or. ivisc_j .eq. 2 .or. ivisc_k .eq. 2)then
                !>
                !>
                !>
                !> computing the vorticity magnitude
                !> and the value store at the cell centre
                !> the array cell_vor_tur(idim-1,jdim-1,kdim-1)
                !> so,the ghost cell or boundary cell must computings
                !> at the turbulent  model subroutine
                !> the turbulent viscosity  -finest meshes only
                call vorticity_magnitude(nbl)
                !>
                !>
                !> calculate the turbulent viscous coeefficient
                !> on the wall,
                !> baldwin-lomax two-layer eddy-viscosity model
                !>
                call baldwin_lomax(nbl)
                !>
                !>
                !>
            elseif(ivisc_i .eq. 3 .or. ivisc_j .eq. 3 .or. ivisc_k .eq. 3)then
                !>
                !>
                !>
                !> computing the vorticity magnitude
                !> and the value store at the cell centre
                !> the array cell_vor_tur(idim-1,jdim-1,kdim-1)
                !> so,the ghost cell or boundary cell must computings
                !> at the turbulent  model subroutine
                !> the turbulent viscosity  -finest meshes only
                call vorticity_magnitude(nbl)
                !>
                !>
                !> calculate the turbulent viscous coeefficient
                !> on the wall,
                !> spalart-allmaras model
                !>
                call spalart_allmaras(nbl)
                !>
                !>
                !>
            elseif(ivisc_i .eq. 4 .or. ivisc_j .eq. 4 .or. ivisc_k .eq. 4)then
                !>
                !>
                !>
                !> other the turbulent model will be export
                !>
                !>
            end if
            !>
            !>
        end if
        !>
        !>
        return
    end subroutine turbulent_viscous_coefficient


