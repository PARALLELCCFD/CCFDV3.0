!   _______________________________________________________________________________
!    _______/\\\\\\\\\_______/\\\\\\\\\______/\\\\\\\\\\\\\\__/\\\\\\\\\\\\\\\______
!     _____/\\\////////_____/\\\////////_____\/\\\//////////__\/\/////////////________
!      ___/\\\/____________/\\\/______________\/\\\____________\/\\\________/\\\\_______
!       __/\\\____________/\\\_________________\/\\\\\\\\\\\\\__\/\\\_______\////________
!        _\/\\\___________\/\\\_________________\/\\\/////////___\/\\\________/\\\\_______
!         _\//\\\__________\//\\\________________\/\\\____________\/\\\_______\////________
!          __\///\\\_________\///\\\______________\/\\\____________\/\\\______/\\\\_________
!           ____\////\\\\\\\\\__\////\\\\\\\\\_____\/\\\____________\/\\\\\\\\\////__________
!            _______\/////////____  \/////////______\///_____________\/////////_______________
!            __________________________________________________________________________________
!----------------------------------------------------------------------------------------------
!>  subroutine input_files_cgns
!>  last edit 2016-05-24
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
	module bc_module
		!>
		!>
!		use global_parameter
		!>
		!>
		implicit none
		!>
		!> public types
		!>
		type ,public :: bc_types
            !>
			!>1-1 cut boundary conditons
			!>bc_type = cut
			!>
			!>  far boundary conditons
			!>bc_type = farfield
			!>
			!> sysmmetry  boundary conditons
			!>bc_type = symmetryplane
			!>
			!>inviscous boundary conditons
			!>bc_type = wallinviscid
			!>
			!>no-slip wall  boundary conditons
			!>bc_type = wallviscid
			character*24 :: bc_type
!			integer :: bc_type
			!>
			!>
			integer :: block_index
			integer :: iblock_index
			!>
			!>
			integer :: istart
			integer :: istart_cut
			!>
			!>
			integer :: iend
			integer :: iend_cut
			!>
			!>
			integer :: jstart
			integer :: jstart_cut
			!>
			!>
			integer :: jend
			integer :: jend_cut
			!>
			!>
			integer :: kstart
			integer :: kstart_cut
			!>
			!>
			integer :: kend
			integer :: kend_cut
			!>
			!>
			integer :: irank
			integer :: jrank
			integer :: krank
			!>
			integer :: norm_index
			integer :: direction
			!>


			!>
		end type bc_types
		!>
!		type(bc_types),pointer :: iibc
		!>
!		type(bc_type),pointer :: bc
		!>
		!>
		!>
		!>
	end module bc_module
