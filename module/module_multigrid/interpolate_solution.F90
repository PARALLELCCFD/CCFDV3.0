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
    subroutine interpolate_solution()
    !***********************************************************************
    !     purpose:interpolate the solutin from a coarser mesh to finer mesh.
    !     when used the full methods and computing at low level and interpolate
    !     the solution to finer mesh and the solution will used as a inii solution
    !     at finer level
    !***********************************************************************
        !>
        !>
        use global_parameter
        use nodes_paras
        use nodes_var
        use nodes_mg
		!>
		implicit none
        !>
        integer :: i,j,k
        integer :: n,level,nbl
        !> interpolate solution to finer mesh - mesh sequence
        !>
        if(imeshseque .le. meshseque .and. imeshseque .gt. 1)then
            if(myid .eq. myhost)then
                write(*,*)'*******beginning prolonged the solutions to finer level *****'
            end if
            !>
            !>
            do nbl = 1,nblocks
                if(level_mg(nbl) .eq. imeshseque-1)then

#if defined PMPI
                    if(myid .eq. n2p(nbl))then
#endif
                        call prolonged(nbl)
#if defined PMPI
                    endif
#endif
                end if
            end do
            !>
            !>
        end if
        !>
        !>
        return
    end subroutine interpolate_solution
