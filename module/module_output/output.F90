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
!>  subroutine output
!>  last edit 2018-05-24
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
    subroutine output()
		!>
		!>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
        use bc_module
        !>
        !>
		implicit none
        integer :: nbl,i,j,k,&
                   idim,jdim,kdim,&
                   ii,sig,&
                   mi,mj,mk

        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        type(bc_types),pointer     :: ibc

        !>
        !>
        mesh => grids(imesh)
        !>
        !>
        !>output the restart files  to the output files
        !>
        !>
        !>
        call restart_output()
        !>
        !couputing the dimension of x,y,z
        !>
        !>
        !>
        do nbl = 1,nblocks,global_level
            !>
            iblock => mesh%blocks(nbl)
            idim  = iblock%idim
            jdim  = iblock%jdim
            kdim  = iblock%kdim
            !>
            open(54,file='ccfd_sys.plt',status='unknown')
            open(55,file='ccfd_wall.plt',status='unknown')
            !> output the wall bc
            !>
            !>
            !>
            do ii=1,num_bc
                ibc => grids(imesh)%bcs(ii)
                if(ibc%bc_type .eq. 'symmetryplane' )then
                    !>
                    if(ibc%block_index .eq. nbl)then
                        write(54,'(1x,a100)') 'variables="x","y","z","density","u-velocity","v-velocity","w-velocity","pressure","blank" '
                        !>
                        !>
                        !>
                        mi = ibc%iend -ibc%istart + 1
                        mj = ibc%jend -ibc%jstart + 1
                        mk = ibc%kend -ibc%kstart + 1
                        sig = 0
                        if(mi == 1)then
                            if (ibc%istart .gt. 2)then
                                sig = 1
                            end if
                            write(54,3001) nbl, mk+1,mj+1,mi
                            if(ialph .eq. 0)then
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%x,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%y,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%z,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                            else
                                write(54,'(3e15.5)') ((( iblock%coordinate(i,j,k)%x,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                                write(54,'(3e15.5)') ((( iblock%coordinate(i,j,k)%z,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                                write(54,'(3e15.5)') (((-iblock%coordinate(i,j,k)%y,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                            endif
                            !>
                            !>
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%r,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%u,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%v,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%w,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%p,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%blank,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                        else if(mj == 1)then
                            if (ibc%jstart .gt. 2)then
                                sig = 1
                            end if
                            write(54,3001) nbl, mk+1,mi+1,mj
                            if(ialph .eq. 0)then
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%x,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%y,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%z,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                            else
                                write(54,'(3e15.5)') ((( iblock%coordinate(i,j,k)%x,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                                write(54,'(3e15.5)') ((( iblock%coordinate(i,j,k)%z,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                                write(54,'(3e15.5)') (((-iblock%coordinate(i,j,k)%y,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                            endif

                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%r,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%u,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%v,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%w,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%p,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%blank,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                        else if(mk == 1)then
                            if (ibc%kstart .gt. 2)then
                                sig = 1
                            end if
                            write(54,3001) nbl, mi+1,mj+1,mk
                            if(ialph .eq. 0)then
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%x,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%y,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                                write(54,'(3e15.5)') (((iblock%coordinate(i,j,k)%z,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                            else
                                write(54,'(3e15.5)') ((( iblock%coordinate(i,j,k)%x,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                                write(54,'(3e15.5)') ((( iblock%coordinate(i,j,k)%z,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                                write(54,'(3e15.5)') (((-iblock%coordinate(i,j,k)%y,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                            endif


                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%r,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%u,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%v,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%w,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%p,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(54,'(3e15.5)') (((iblock%cell(i,j,k)%blank,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                        end if
                    end if
                end if
                !>
                !>
                if( ibc%bc_type .eq. 'wallinviscid' .or. ibc%bc_type .eq. 'wallviscid')then
                    !>
                    if(ibc%block_index .eq. nbl)then
                        write(55,'(1x,a100)') 'variables="x","y","z","density","u-velocity","v-velocity","w-velocity","pressure","blank" '
                        !>
                        !>
                        !>

                        mi = ibc%iend -ibc%istart + 1
                        mj = ibc%jend -ibc%jstart + 1
                        mk = ibc%kend -ibc%kstart + 1

                        sig = 0
                        if(mi == 1)then
                            if (ibc%istart .gt. 2)then
                                sig = 1
                            end if
                            write(55,1001) nbl, mk+1,mj+1,mi
                            if(ialph .eq. 0)then
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%x,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%y,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%z,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                            else
                                write(55,'(3e15.5)') ((( iblock%coordinate(i,j,k)%x,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                                write(55,'(3e15.5)') ((( iblock%coordinate(i,j,k)%z,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                                write(55,'(3e15.5)') (((-iblock%coordinate(i,j,k)%y,k=ibc%kstart,ibc%kend+1),j=ibc%jstart,ibc%jend+1),i=ibc%istart+sig,ibc%iend+sig)
                            endif

                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%r,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%u,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%v,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%w,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%p,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%blank,k=ibc%kstart,ibc%kend),j=ibc%jstart,ibc%jend),i=ibc%istart,ibc%iend)
                        else if(mj == 1)then
                            if (ibc%jstart .gt. 2)then
                                sig = 1
                            end if
                            write(55,1001) nbl, mk+1,mi+1,mj
                            if(ialph .eq. 0)then
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%x,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%y,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%z,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                            else
                                write(55,'(3e15.5)') ((( iblock%coordinate(i,j,k)%x,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                                write(55,'(3e15.5)') ((( iblock%coordinate(i,j,k)%z,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                                write(55,'(3e15.5)') (((-iblock%coordinate(i,j,k)%y,k=ibc%kstart,ibc%kend+1),i=ibc%istart,ibc%iend+1),j=ibc%jstart+sig,ibc%jend+sig)
                            endif
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%r,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%u,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%v,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%w,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%p,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%blank,k=ibc%kstart,ibc%kend),i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend)
                        else if(mk == 1)then
                            if (ibc%kstart .gt. 2)then
                                sig = 1
                            end if
                            write(55,1001) nbl, mi+1,mj+1,mk
                            if(ialph .eq. 0)then
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%x,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%y,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                                write(55,'(3e15.5)') (((iblock%coordinate(i,j,k)%z,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                            else
                                write(55,'(3e15.5)') ((( iblock%coordinate(i,j,k)%x,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                                write(55,'(3e15.5)') ((( iblock%coordinate(i,j,k)%z,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                                write(55,'(3e15.5)') (((-iblock%coordinate(i,j,k)%y,i=ibc%istart,ibc%iend+1),j=ibc%jstart,ibc%jend+1),k=ibc%kstart+sig,ibc%kend+sig)
                            endif
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%r,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%u,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%v,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%w,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%p,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                            write(55,'(3e15.5)') (((iblock%cell(i,j,k)%blank,i=ibc%istart,ibc%iend),j=ibc%jstart,ibc%jend),k=ibc%kstart,ibc%kend)
                        end if
                    end if
                end if
                    !>
                    !>
            end do !>do ii=1,num_bc
        end do !> nbl=1,nblocks
        !>
        !>
1001    format('ZONE T="WALL_BLK',i6,' " ',',I=',I4,',J=',I4,',K=',I4, ',DATAPACKING=BLOCK, VARLOCATION=([4-9]=CELLCENTERED) ')
3001    format('ZONE T="SYS_BLK',i6,' " ',',I=',I4,',J=',I4,',K=',I4, ',DATAPACKING=BLOCK, VARLOCATION=([4-9]=CELLCENTERED) ')
        close(54)
        close(55)
        !>output the global filed result at  finest meshs
        !>
        !>
        open(56,file='ccfd_global.plt',status='unknown')
        do nbl=1,nblocks,global_level
            iblock => mesh%blocks(nbl)
            write(56,*) 'variables="x","y","z","density","u-velocity","v-velocity","w-velocity","pressure" '
            write(56,1002) nbl,iblock%idim,iblock%jdim,iblock%kdim
            write(56,'(5e14.6)') ((( iblock%coordinate(i,j,k)%x , i=2,iblock%idim + 1),j=2,iblock%jdim + 1),k=2,iblock%kdim + 1)
            write(56,'(5e14.6)') ((( iblock%coordinate(i,j,k)%y , i=2,iblock%idim + 1),j=2,iblock%jdim + 1),k=2,iblock%kdim + 1)
            write(56,'(5e14.6)') ((( iblock%coordinate(i,j,k)%z , i=2,iblock%idim + 1),j=2,iblock%jdim + 1),k=2,iblock%kdim + 1)
            !>
            !>
            write(56,'(5e14.6)') (((iblock%cell(i,j,k)%r,i=2,iblock%idim),j=2,iblock%jdim),k=2,iblock%kdim )
            write(56,'(5e14.6)') (((iblock%cell(i,j,k)%u,i=2,iblock%idim),j=2,iblock%jdim),k=2,iblock%kdim )
            write(56,'(5e14.6)') (((iblock%cell(i,j,k)%v,i=2,iblock%idim),j=2,iblock%jdim),k=2,iblock%kdim )
            write(56,'(5e14.6)') (((iblock%cell(i,j,k)%w,i=2,iblock%idim),j=2,iblock%jdim),k=2,iblock%kdim )
            write(56,'(5e14.6)') (((iblock%cell(i,j,k)%p,i=2,iblock%idim),j=2,iblock%jdim),k=2,iblock%kdim )
        end do
1002    format('zone t="blk',i6,' " ',',i=',i5,',j=',i5,',k=',i5, ',datapacking=block, varlocation=([4-8]=cellcentered) ')
        close(56)


        return
    end subroutine output


