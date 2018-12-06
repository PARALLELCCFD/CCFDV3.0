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
    subroutine boundary_check
        !>
        !>
        !> check the boundary and set some informations
        !>

        use global_parameter
        use mesh_overlap_module
        use bc_module
        use blocks_module
        use nodes_paras
        use nodes_var
        use nodes_var_bc
        implicit none
#if defined PMPI
        include "mpif.h"
#endif
        integer::i,j,k,ii,&
                 nbl,inbl,&
                 is,ie,&
                 js,je,&
                 ks,ke
        integer::block_index,&
                 iblock_index,&
                 indexp,&
                 ierr ! mpi
        integer :: iid ,jjd,kkd,ipoint,bc_type
        character(len=20) iface
        integer,dimension(:),allocatable   :: bc_send
        !>
        !>
        type(bc_types),pointer :: ibc
        type(overlap_type),pointer :: mesh
        !>
        !       boundary conditions  check
        !
        !------------------------------------------------------------
        !       check boundary conditions  informations
        !       before comouting the mesh
        !------------------------------------------------------------
        !------------------------------------------------------------
        mesh => grids(imesh)

#if defined PMPI
        if(myid .eq. myhost)then
#endif

            do ii = 1,num_bc
                !>
                !>
                ibc => mesh%bcs(ii)
                is  = ibc%istart
                ie  = ibc%iend
                js  = ibc%jstart
                je  = ibc%jend
                ks  = ibc%kstart
                ke  = ibc%kend
                nbl = ibc%block_index
                inbl= ibc%iblock_index
                !>
                !> the bc not at the local block
                !>
                if(is .eq. ie)then
                    iface = 'i'
                    indexp = is
                else if(js .eq. je)then
                    iface = 'j'
                    indexp = js
                else if(ks .eq. ke)then
                    iface = 'k'
                    indexp = ks
                else
                    write(*,*) 'ccfdv3.0 :: have some wrong of the grid ,the bc maybe not a face!'
                    write(*,'("bc=",i4," type=",a24," range=",6i6)') ii,ibc%bc_type,is,ie,js,je,ks,ke
                    stop
                end if
                    !>
                    !>
                if(iface  .eq. 'i')then
                    if(indexp .le. 2)then
                        !>
                        ibc%direction = -1
                        !>
                        !> the inviscous wall or no-slip viscous wall
                        !>
                        if(ibc%bc_type .eq. 'wallinviscid' .or. ibc%bc_type .eq. 'wallviscid')then
                            ibc%norm_index  = 1
                        else
                            ibc%norm_index  = -1
                        end if
                    else
                        !>
                        ibc%direction = 1
                        !>
                        !> the inviscous wall or no-slip viscous wall
                        !>
                        if(ibc%bc_type .eq. 'wallinviscid' .or. ibc%bc_type .eq. 'wallviscid')then
                            ibc%norm_index  = -1
                        else
                            ibc%norm_index  = 1
                        end if
                    end if
                else if(iface  .eq. 'j')then
                    if(indexp .le. 2)then
                        !>
                        ibc%direction = -1
                        !>
                        !> the inviscous wall or no-slip viscous wall
                        !>
                        if(ibc%bc_type .eq. 'wallinviscid' .or. ibc%bc_type .eq. 'wallviscid')then
                            ibc%norm_index  = 1
                        else
                            ibc%norm_index  = -1
                        end if
                    else
                        !>
                        ibc%direction = 1
                        !>
                        !> the inviscous wall or no-slip viscous wall
                        !>
                        if(ibc%bc_type .eq. 'wallinviscid' .or. ibc%bc_type .eq. 'wallviscid')then
                            ibc%norm_index  = -1
                        else
                            ibc%norm_index  = 1
                        end if
                    end if
                else if(iface  .eq. 'k')then
                    if(indexp .le. 2)then
                        !>
                        ibc%direction = -1
                        !>
                        !> the inviscous wall or no-slip viscous wall
                        !>
                        if(ibc%bc_type .eq. 'wallinviscid' .or. ibc%bc_type .eq. 'wallviscid')then
                            ibc%norm_index  = 1
                        else
                            ibc%norm_index  = -1
                        end if
                    else
                        !>
                        ibc%direction = 1
                        !>
                        !> the inviscous wall or no-slip viscous wall
                        !>
                        if(ibc%bc_type .eq. 'wallinviscid' .or. ibc%bc_type .eq. 'wallviscid')then
                            ibc%norm_index  = -1
                        else
                            ibc%norm_index  = 1
                        end if
                    end if
                else
                    write(*,*) 'ccfdv3.0 :: have some error at the bc face ,check the bc '
                    stop
                end if
                !>
                !> i face
                !>
                !>reset the dimension for the bc
                !>
                !>
                !> i face
                !>
                !>reset the dimension for the bc
                !>
                if(is .lt. ie )then
                    ibc%istart = ibc%istart  + 1
                else if(is .eq. ie)then
                    if(indexp .le. 1)then
                        ibc%istart = ibc%istart +1
                        ibc%iend   = ibc%iend +1
                    end if
                else
                    ibc%iend   = ibc%iend +1
                end if
                !>
                !>reset the dimension for the 1to1 bc
                !>
                if(ibc%istart_cut .lt. ibc%iend_cut )then
                    ibc%istart_cut = ibc%istart_cut  + 1
                else if(ibc%istart_cut .eq. ibc%iend_cut)then
                    if(ibc%istart_cut .le. 1)then
                        ibc%istart_cut = ibc%istart_cut +1
                        ibc%iend_cut   = ibc%iend_cut +1
                    end if
                else
                    ibc%iend_cut   = ibc%iend_cut +1
                end if
                !>
                !>
                !> j face
                !>reset the dimension for the bc
                !>
                if(js .lt. je )then
                    ibc%jstart = ibc%jstart  + 1
                else if(js .eq. je)then
                    if(indexp .le. 1)then
                        ibc%jstart = ibc%jstart +1
                        ibc%jend   = ibc%jend +1
                    end if
                else
                    ibc%jend   = ibc%jend +1
                end if
                !>
                !>reset the dimension for the 1to1 bc
                !>
                if(ibc%jstart_cut .lt. ibc%jend_cut )then
                    ibc%jstart_cut = ibc%jstart_cut  + 1
                else if(ibc%jstart_cut .eq. ibc%jend_cut)then
                    if(ibc%jstart_cut .le. 1)then
                        ibc%jstart_cut = ibc%jstart_cut +1
                        ibc%jend_cut   = ibc%jend_cut +1
                    end if
                else
                    ibc%jend_cut   = ibc%jend_cut +1
                end if
                !>
                !>
                !> k face
                !>reset the dimension for the bc
                !>
                if(ks .lt. ke )then
                    ibc%kstart = ibc%kstart  + 1
                else if(ks .eq. ke)then
                    if(indexp .le. 1)then
                        ibc%kstart = ibc%kstart +1
                        ibc%kend   = ibc%kend +1
                    end if
                else
                    ibc%kend   = ibc%kend +1
                end if
                !>
                !>reset the dimension for the 1to1 bc
                !>
                if(ibc%kstart_cut .lt. ibc%kend_cut )then
                    ibc%kstart_cut = ibc%kstart_cut  + 1
                else if(ibc%kstart_cut .eq. ibc%kend_cut)then
                    if(ibc%kstart_cut .le. 1)then
                        ibc%kstart_cut = ibc%kstart_cut +1
                        ibc%kend_cut   = ibc%kend_cut +1
                    end if
                else
                    ibc%kend_cut   = ibc%kend_cut +1
                end if
                write(*,'("num_bc=",7i6)') ii,ibc%istart,ibc%iend,ibc%jstart,ibc%jend,ibc%kstart,ibc%kend
                !>
                !> end the num_bc cycle
            end do
# if defined PMPI
        end if
        !> bcast the bc information from host processor to other processors
        !> shared the bc informations
        !>

        allocate(bc_send(20*num_bc))
        !>
        !>
!
        if(myid .eq. myhost)then
            ipoint  = 1
            do ii=1,num_bc
                if(mesh%bcs(ii)%bc_type .eq.  'cut')then
                    bc_type = 315
                else if(mesh%bcs(ii)%bc_type .eq.  'farfield')then
                    bc_type = 316
                else if(mesh%bcs(ii)%bc_type .eq.  'symmetryplane')then
                    bc_type = 317
                else if(mesh%bcs(ii)%bc_type .eq.  'wallinviscid')then
                    bc_type = 318
                else if(mesh%bcs(ii)%bc_type .eq.  'wallviscid')then
                    bc_type = 319
                end if
                bc_send(ipoint)= bc_type
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%block_index
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%iblock_index
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%istart
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%istart_cut
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%jstart
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%jstart_cut
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%kstart
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%kstart_cut
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%iend
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%iend_cut
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%jend
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%jend_cut
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%kend
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%kend_cut
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%irank
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%jrank
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%krank
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%norm_index
                ipoint = ipoint + 1
                bc_send(ipoint)= mesh%bcs(ii)%direction
                ipoint = ipoint + 1
!                write(100+myid,'("bc",20i4,2x,a24)') ii,mesh%bcs(ii)%block_index,mesh%bcs(ii)%iblock_index,&
!                mesh%bcs(ii)%istart,mesh%bcs(ii)%istart_cut,mesh%bcs(ii)%iend,mesh%bcs(ii)%iend_cut,&
!                mesh%bcs(ii)%jstart,mesh%bcs(ii)%jstart_cut,mesh%bcs(ii)%jend,mesh%bcs(ii)%jend_cut,&
!                mesh%bcs(ii)%kstart,mesh%bcs(ii)%kstart_cut,mesh%bcs(ii)%kend,mesh%bcs(ii)%kend_cut,&
!                mesh%bcs(ii)%irank,mesh%bcs(ii)%jrank,mesh%bcs(ii)%krank,mesh%bcs(ii)%norm_index,mesh%bcs(ii)%direction,&
!                mesh%bcs(ii)%bc_type
            end do
        end if
        !>
        !>
        call mpi_bcast(bc_send,20*num_bc,mpi_integer,myhost,mycomm,ierr)
        !>
        !>
        call mpi_barrier(mycomm,ierr)
		!>
		!>
		if(myid .ne. myhost )then
            ipoint  = 1
            do ii = 1,num_bc
                bc_type = bc_send(ipoint)
                ipoint  = ipoint  + 1
                if(bc_type .eq.  315)then
                    mesh%bcs(ii)%bc_type = 'cut'
                else if(bc_type .eq.  316)then
                    mesh%bcs(ii)%bc_type = 'farfield'
                else if(bc_type .eq.  317)then
                    mesh%bcs(ii)%bc_type = 'symmetryplane'
                else if(bc_type .eq.  318)then
                    mesh%bcs(ii)%bc_type = 'wallinviscid'
                else if(bc_type .eq.  319)then
                    mesh%bcs(ii)%bc_type = 'wallviscid'
                end if
                mesh%bcs(ii)%block_index = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%iblock_index = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%istart = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%istart_cut = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%jstart = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%jstart_cut = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%kstart = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%kstart_cut = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%iend = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%iend_cut = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%jend = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%jend_cut = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%kend = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%kend_cut = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%irank = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%jrank = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%krank = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%norm_index = bc_send(ipoint)
                ipoint  = ipoint  + 1
                mesh%bcs(ii)%direction = bc_send(ipoint)
                ipoint  = ipoint  + 1
!                write(100+myid,'("bc",20i4,2x,a24)') ii,mesh%bcs(ii)%block_index,mesh%bcs(ii)%iblock_index,&
!                mesh%bcs(ii)%istart,mesh%bcs(ii)%istart_cut,mesh%bcs(ii)%iend,mesh%bcs(ii)%iend_cut,&
!                mesh%bcs(ii)%jstart,mesh%bcs(ii)%jstart_cut,mesh%bcs(ii)%jend,mesh%bcs(ii)%jend_cut,&
!                mesh%bcs(ii)%kstart,mesh%bcs(ii)%kstart_cut,mesh%bcs(ii)%kend,mesh%bcs(ii)%kend_cut,&
!                mesh%bcs(ii)%irank,mesh%bcs(ii)%jrank,mesh%bcs(ii)%krank,mesh%bcs(ii)%norm_index,mesh%bcs(ii)%direction,&
!                mesh%bcs(ii)%bc_type
            end do

        end if
        !>
        !>
        deallocate(bc_send)
#endif
        !>
        !>
        !>
        !>
        return

    end subroutine boundary_check
  !
