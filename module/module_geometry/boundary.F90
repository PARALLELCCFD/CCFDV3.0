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
    subroutine boundary(nbl,level)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use bc_module
        use nodes_var
        use nodes_paras
        use nodes_var_bc
        use nodes_mg
		!>
		!>
		implicit none

        type(bc_types),pointer :: ibc
        integer :: nbl,&
                   num_of_bc,&
                   level,&
                   iblock,&
                   nd_src
        type(overlap_type),pointer :: mesh
        !       boundary conditions
        !
        !------------------------------------------------------------
        !       boundary conditions
        !------------------------------------------------------------
        mesh => grids(imesh)
        do num_of_bc = 1,num_bc
            !>
            !>
            ibc => mesh%bcs(num_of_bc)
            iblock = ibc%block_index
            !>
            !>
            !>
            if(level .eq. level_mg(iblock) .and.  iblock .eq. nbl)then

#if defined PMPI
                nd_src = n2p(iblock)
                !>
                !>
                if(myid .eq. nd_src)then
#endif
                    if(iblock .ne. nbl) cycle
                    !>
!                    write(*,'("in the boundary sub ,and the num of bc=",i4)')num_of_bc
                    if(ibc%bc_type .eq. 'cut')then
                        !>
#if defined PMPI
                        !> null
#else
                        !>
                        !> the cut boundary conditions
!                        write(*,*)'cut=',icyc,num_of_bc
                        call bc_cut(num_of_bc)
!                        write(*,*) 'bc_cut=',nbl,num_of_bc
#endif
                    else if(ibc%bc_type .eq. 'farfield')then
                        !>
                        !>for the far flow conditions
                        !>
!                        write(*,*)'farfield=',icyc,num_of_bc
                        call bc_far(num_of_bc,nbl)
!                        write(*,*) 'farfield=',nbl,num_of_bc
                    else if(ibc%bc_type .eq. 'symmetryplane')then
                        !>
                        !> for the symmetry boundary conditions
                        !>
!                        write(*,*)'symmetryplane=',icyc,num_of_bc
                        call bc_symmetry(num_of_bc,nbl)
!                        write(*,*) 'symmetryplane=',nbl,num_of_bc
                    else if(ibc%bc_type .eq. 'wallinviscid')then
                        !>
                        !> for viscous wall boundary conditions
                        !> for the n-s
!                        write(*,*)'wallinviscid=',icyc,num_of_bc
                        call bc_wall(num_of_bc,nbl)
!                        write(*,*) 'wallinviscid=',nbl,num_of_bc
                    else if(ibc%bc_type .eq. 'wallviscid')then
                        !>
                        !> for viscous wall boundary conditions
                        !> for the n-s
!                        write(*,*)'wallviscid=',icyc,num_of_bc
                        call bc_wall(num_of_bc,nbl)
!                        write(*,*) 'wallviscid=',nbl,num_of_bc
                    else
                    end if
                    !>
                    !>
                    !>
#if defined PMPI
                end if
#endif
            end if

        end do

        return
    end subroutine boundary
