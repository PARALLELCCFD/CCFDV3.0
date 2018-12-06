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
    subroutine bc_cut_array_allocate()
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use bc_module
        use nodes_paras
        use nodes_var_bc
        use nodes_var
		!>
		!>
		implicit none
#if defined PMPI
        include "mpif.h"
#endif
#if defined PMPI
        type(overlap_type),pointer :: mesh
        integer :: i_size,j_size,k_size
        integer :: mm,dblock,sblock,ixyz,i
        !*************************************************************
        !>
        !>
        i_size1to1 = 0
        mesh => grids(imesh)
        do mm=1,num_bc
            i_size = 1
            j_size = 1
            k_size = 1

            if(mesh%bcs(mm)%bc_type .eq. 'cut')then
                !>
                !>
                sblock  = n2p(mesh%bcs(mm)%block_index)
                dblock  = n2p(mesh%bcs(mm)%iblock_index)
                !>
                !>
                if( myid .eq. sblock .and.  myid .ne. dblock )then
                    !>
                    !>
                    i_size = i_size + abs(mesh%bcs(mm)%iend -mesh%bcs(mm)%istart)
                    j_size = j_size + abs(mesh%bcs(mm)%jend -mesh%bcs(mm)%jstart)
                    k_size = k_size + abs(mesh%bcs(mm)%kend -mesh%bcs(mm)%kstart)
                    i_size1to1 = i_size1to1 + i_size*j_size*k_size
                end if
                !>
                !>
            end if
        end do
        i_size1to1 = i_size1to1*2*16
        !write(*,'("size of buffer=",i20," myid =",i6)') i_size1to1,myid
        call allocate_nodes_var_bc(num_bc,i_size1to1,global_level)
        !*************************************************************
#endif
        return
    end subroutine bc_cut_array_allocate
  !**********************************************
  !**********************************************
  !**********************************************
