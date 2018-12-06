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
    subroutine restart_output()
		!>
		!>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
        !>
		implicit none
        integer :: nbl,i,j,k

        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh

        !>
        !>
        mesh => grids(imesh)

        open (unit=63,file = 'ccfd.restart', form='formatted',status='unknown')
        !>
        !>
        write(63,'(3i6)') nblocks,icyc-1,global_level
        do nbl =1,nblocks,global_level
            iblock => mesh%blocks(nbl)
            write(63,'(3i6)') iblock%idim,iblock%jdim,iblock%kdim
        end do
        !>
        !>
        write(63,'(3e24.12)') mach,alpha,reue
        !>
        !>
        do nbl =1,nblocks,global_level
            iblock => mesh%blocks(nbl)
            write(63,'(5e24.12)') (((iblock%cell(i,j,k)%r,i=2,iblock%idim),j=2,iblock%jdim),k=1,iblock%kdim)
            write(63,'(5e24.12)') (((iblock%cell(i,j,k)%u,i=2,iblock%idim),j=2,iblock%jdim),k=1,iblock%kdim)
            write(63,'(5e24.12)') (((iblock%cell(i,j,k)%v,i=2,iblock%idim),j=2,iblock%jdim),k=1,iblock%kdim)
            write(63,'(5e24.12)') (((iblock%cell(i,j,k)%w,i=2,iblock%idim),j=2,iblock%jdim),k=1,iblock%kdim)
            write(63,'(5e24.12)') (((iblock%cell(i,j,k)%p,i=2,iblock%idim),j=2,iblock%jdim),k=1,iblock%kdim)

        end do

        close(63)


        return
    end subroutine restart_output


