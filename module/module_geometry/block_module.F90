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
	module blocks_module
		use global_parameter
		use cell_module
		use bc_module
		use turbulent_module
		use coordinate_module
		use metric_module
		use variable_module
		use variable_mg_module
		!>
		!>
		implicit none
		!>
		!>
		!>

		!>
		type ,public ::blocks_type
			!>
			integer :: idim=0
			integer :: jdim=0
			integer :: kdim=0
			!>
			!>
			!> store the flow variable
			!> like density,pressure,velocity of three directions

			type(cell_type),pointer,dimension(:,:,:) :: cell!=> null()
			!>
			!> store the coordinate of the three directions
			!>
			type(coordinate_type),pointer,dimension(:,:,:) :: coordinate!=> null()
			!>
			!>
			!>store the metric about nine vector
			type(metric_type),pointer,dimension(:,:,:) :: metric!=> null()
			!>
			!> the module for the variable
			!>
			!> res_1,...,res_5 of the right hand term
            type(variable_types),pointer,dimension(:,:,:) :: variable!=> null()
            !>
            !>
            type(variable_mg_types),pointer,dimension(:,:,:) :: variable_mg!=> null()
			!>
			!> if the visc ==2 used the b-l model
			!> if the visc ==3 used the s-a model
            type(turbulent_type),pointer,dimension(:,:,:) :: turbulent!=> null()
			!>
		end type blocks_type
		!>
!		!>
!		interface new
!            module procedure mem_allocate_block
!        end interface
!        !>
!        !>
!        interface delete
!            module procedure demem_allocate_block
!        end interface
		!>
		!>
!		type(nblocks_type),pointer,dimension(:) :: nblocks
        contains
		!>
        !>
        !>
        !>
        subroutine mem_allocate_block(nbl,i,j,k)
            implicit none
            type(blocks_type) :: nbl
            integer :: i,j,k
            integer :: ghost_i,ghost_j,ghost_k
            integer :: error
            error = 0
            !>
            !>
            ghost_i = i + 2
            ghost_j = j + 2
            ghost_k = k + 2
            !>
            !>
            !>
            nbl%idim = i
            nbl%jdim = j
            nbl%kdim = k
            !>
            !>
            allocate(nbl%cell(0:ghost_i,0:ghost_j,0:ghost_k),stat = error)
            !>
            !>
            if(error .ne. 0 )then
                write(*,*) 'when allocate the cell array ,allocate th bytes too large in mem_allocate_cell  '
                stop
            end if
            !>
            !>
            if(grid_dim .eq. 3)then
                allocate(nbl%coordinate(0:ghost_i,0:ghost_j,0:ghost_k),stat = error)
                if(error .ne. 0 )then
                    write(*,*) 'when allocate the coordinate array ,allocate th bytes too large in mem_allocate_cell  '
                    stop
                end if
            else
                allocate(nbl%coordinate(1:2,0:ghost_j,0:ghost_k),stat = error)
                if(error .ne. 0 )then
                    write(*,*) 'when allocate the coordinate array ,allocate th bytes too large in mem_allocate_cell  '
                    stop
                end if
            end if
            !>
            !>
            allocate(nbl%metric(0:ghost_i+1,0:ghost_j+1,0:ghost_k+1),stat = error)
            if(error .ne. 0 )then
                write(*,*) 'when allocate the metric array ,allocate th bytes too large in mem_allocate_cell  '
                stop
            end if
            !>
            !> store the variable for  multigrid methods
            !>

            allocate(nbl%variable(0:ghost_i,0:ghost_j,0:ghost_k),stat= error)
            if(error .ne. 0 )then
                write(*,*) 'when allocate the variable array ,allocate th bytes too large in mem_allocate_cell'
                stop
            end if
            !>
            !>
            if(mgflags .ne. 0)then
                allocate(nbl%variable_mg(0:ghost_i,0:ghost_j,0:ghost_k),stat= error)
                if(error .ne. 0 )then
                    write(*,*) 'when allocate the variable array ,allocate th bytes too large in mem_allocate_cell'
                    stop
                end if
            end if
            !>
            !>
            if(ivisc_i >= 2 .or. ivisc_j >= 2 .or. ivisc_k >=2 )then
                allocate(nbl%turbulent(0:ghost_i,0:ghost_j,0:ghost_k),stat = error)
                if(error .ne. 0 )then
                    write(*,*) 'when allocate the turbulent array ,allocate th bytes too large in mem_allocate_cell  '
                    stop
                end if
            end if
            !>

        end subroutine mem_allocate_block
        !>
        !>
        subroutine demem_allocate_block(nbl)
            implicit none
            !>
            !>
            type(blocks_type) :: nbl
            !>
            !>
            if(associated(nbl%cell))         deallocate(nbl%cell)
            if(associated(nbl%coordinate))   deallocate(nbl%coordinate)
            if(associated(nbl%metric))       deallocate(nbl%metric)
            if(associated(nbl%variable))     deallocate(nbl%variable)
            if(mgflags .ne. 0 )then
                if(associated(nbl%variable_mg))  deallocate(nbl%variable_mg)
            end if
            if(ivisc_i >= 2 .or. ivisc_j >= 2 .or. ivisc_k >=2 )then
                if(associated(nbl%turbulent))    deallocate(nbl%turbulent)
            end if
            !>
            !>

        end subroutine
	end module blocks_module
