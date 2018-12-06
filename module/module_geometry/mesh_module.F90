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
!>  module  mesh_overlap_module
!>  the overlap grid and the blocks at each overlap meshes
!>
!----------------------------------------------------------------------------------------------
	module mesh_overlap_module
		!>
		!>
		use blocks_module
		use bc_module
		!>
		!>
		implicit none
		!>
		!>
		public
		!>
		type  :: overlap_type
			!>
			type(blocks_type),dimension(:),pointer :: blocks!=>null()
			!> store the informations of the boundary conditions
			!> about type,start,end....
			type(bc_types),   dimension(:),pointer :: bcs!=>null()
			!>
		end type overlap_type
		!>
		!>
		type(overlap_type),dimension(:),pointer :: grids!=> null()
		!>
!		!>
!		interface new
!            module procedure mem_allocate_mesh_blocks
!            module procedure mem_allocate_mesh_bc
!        end interface
        !>
        !>
!        interface delete
!            module procedure demem_allocate_block
!        end interface
		!>
		!>
        contains
        !>
        !>
!        subroutine mem_allocate_mesh_blocks(num,iblock)
!            implicit none
!            !>
!            type(blocks_type),dimension(:),pointer :: iblock
!            integer ::num
!            write(*,*) 'in the mesh module=',num
!            !>
!            allocate(grids(imesh)%blocks(num))
!            iblock => grids(imesh)%blocks
!            if(allocated(grids(imesh)%blocks(num))) write(*,*)'allocate the blocks array'
!            mesh%blocks(1)%idim = 1
!            write(*,*) mesh%blocks(1)%idim
!            write(*,*) 'in the mesh module='
!
!        end subroutine
!        subroutine mem_allocate_mesh_bc(mesh,num)
!            implicit none
!            !>
!            type(overlap_type) :: mesh
!            integer ::num
!            !>
!            allocate(mesh%bcs(num))
!
!        end subroutine
		!>



	end module mesh_overlap_module
