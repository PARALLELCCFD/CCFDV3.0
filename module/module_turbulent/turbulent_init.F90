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
!>  purpose :: initialize the turbulent initial conditions
!>  last edit 2018-04-16
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
    subroutine turbulent_init()
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use nodes_var
        use nodes_paras
		implicit none
#if defined PMPI
        include "mpif.h"
#endif
        type(overlap_type),pointer  ::  mesh
        type(blocks_type),pointer   ::  iblock
        integer                     ::  i,j,k
        integer                     ::  nbl
        real(kind =dprec)           ::  chi
        !>
        !>
        chi  =  1.341946
        mesh => grids(imesh)
        !>
        !>
        do nbl=1,nblocks
#if defined PMPI
            if( n2p(nbl) .eq. myid)then
#endif
                iblock => mesh%blocks(nbl)
                if(ivisc_i .gt. 2 .or. ivisc_j .gt. 2 .or.  ivisc_k .gt. 2)then
                    !>
                    !>
                    !> fill in edge values for turbulent array for safety before passing the data
                    do k=2,iblock%kdim
                        do j=2,iblock%jdim
                            do i=2,iblock%idim
                                iblock%turbulent(i,j,k)%tur_save = chi
                                iblock%turbulent(i,j,k)%viscous  = chi*(chi**3/(chi**3 + 357.911))
                            end do
                        end do
                    end do
                else
                    do k=2,iblock%kdim
                        do j=2,iblock%jdim
                            do i=2,iblock%idim
                                iblock%turbulent(i,j,k)%viscous = 0.0
                            end do
                        end do
                    end do
                end if
#if defined PMPI
            end if
#endif
        end do
        return

	end subroutine turbulent_init
