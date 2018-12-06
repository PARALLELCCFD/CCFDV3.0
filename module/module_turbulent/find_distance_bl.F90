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
    subroutine find_distance_bl()
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use metric_module
        use coordinate_module
        use nodes_paras
        use nodes_var
        use nodes_var_bc
        use nodes_mg
		implicit none
#if defined PMPI
        include "mpif.h"
#endif
        type(overlap_type),pointer :: mesh
        real(kind= dprec),dimension(:,:,:,:),allocatable :: work_for_tmp
        type(blocks_type),pointer  :: iblock
        integer::i,j,k
        integer::idim,jdim,kdim
        integer::nbl
        !>
        !>
        mesh => grids(imesh)
        !> evaluate directed distance from wall face to inner cell
        !> for use in evaluating the baldwin-lomax turbulence model
        !> **********baldwin-lomax turbulence model****************
        if(ivisc_i .eq. 2 .or. ivisc_j .eq. 2 .or. ivisc_k .eq. 2)then
        !>
        do nbl=1 ,nblocks
            iblock => mesh%blocks(nbl)
#if defined PMPI
            if(myid .eq. n2p(nbl))then
#endif

                !>
                !>
                idim  =   iblock%idim
                jdim  =   iblock%jdim
                kdim  =   iblock%kdim
                !>
                !>
                allocate(work_for_tmp(idim,jdim,kdim,4))
                !> cell-center locations
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            !>x -dimension
                            work_for_tmp(i,j,k,1)= 0.125*(iblock%coordinate(i  ,j  ,k  )%x + &
                                                          iblock%coordinate(i+1,j  ,k  )%x + &
                                                          iblock%coordinate(i  ,j+1,k  )%x + &
                                                          iblock%coordinate(i  ,j  ,k+1)%x + &
                                                          iblock%coordinate(i+1,j+1,k  )%x + &
                                                          iblock%coordinate(i  ,j+1,k+1)%x + &
                                                          iblock%coordinate(i+1,j  ,k+1)%x + &
                                                          iblock%coordinate(i+1,j+1,k+1)%x   )
                            !>
                            !>
                            !>y -dimension
                            work_for_tmp(i,j,k,2)= 0.125*(iblock%coordinate(i  ,j  ,k  )%y + &
                                                          iblock%coordinate(i+1,j  ,k  )%y + &
                                                          iblock%coordinate(i  ,j+1,k  )%y + &
                                                          iblock%coordinate(i  ,j  ,k+1)%y + &
                                                          iblock%coordinate(i+1,j+1,k  )%y + &
                                                          iblock%coordinate(i  ,j+1,k+1)%y + &
                                                          iblock%coordinate(i+1,j  ,k+1)%y + &
                                                          iblock%coordinate(i+1,j+1,k+1)%y   )
                            !>z -dimension
                            work_for_tmp(i,j,k,3)= 0.125*(iblock%coordinate(i  ,j  ,k  )%z + &
                                                          iblock%coordinate(i+1,j  ,k  )%z + &
                                                          iblock%coordinate(i  ,j+1,k  )%z + &
                                                          iblock%coordinate(i  ,j  ,k+1)%z + &
                                                          iblock%coordinate(i+1,j+1,k  )%z + &
                                                          iblock%coordinate(i  ,j+1,k+1)%z + &
                                                          iblock%coordinate(i+1,j  ,k+1)%z + &
                                                          iblock%coordinate(i+1,j+1,k+1)%z   )
                        end do
                    end do
                end do
                !>
                !>**********************************************************************************************
                !> the min-distance on i0 wall face
                !> face-center locations
                !>
                i = 2
                do k=2,kdim
                    do j=2,jdim
                        work_for_tmp(i,j,k,4) = 0.25e0*((  iblock%coordinate(i,j,k)%x      + &
                                                           iblock%coordinate(i,j+1,k)%x    + &
                                                           iblock%coordinate(i,j,k+1)%x    + &
                                                           iblock%coordinate(i,j+1,k+1)%x    &
                                                         )*iblock%metric(i,j,k)%fi + &
                                                         ( iblock%coordinate(i,j,k)%y     + &
                                                           iblock%coordinate(i,j+1,k)%y   + &
                                                           iblock%coordinate(i,j,k+1)%y   + &
                                                           iblock%coordinate(i,j+1,k+1)%y   &
                                                         )*iblock%metric(i,j,k)%fj + &
                                                         ( iblock%coordinate(i,j,k)%z     + &
                                                           iblock%coordinate(i,j+1,k)%z   + &
                                                           iblock%coordinate(i,j,k+1)%z   + &
                                                           iblock%coordinate(i,j+1,k+1)%z &
                                                          )*iblock%metric(i,j,k)%fk)
                    end do
                end do
                !>
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            iblock%turbulent(i,j,k)%distance_tur_1  = &
                                                            (work_for_tmp(i,j,k,1)*iblock%metric(2,j,k)%fi + &
                                                             work_for_tmp(i,j,k,2)*iblock%metric(2,j,k)%fj + &
                                                             work_for_tmp(i,j,k,3)*iblock%metric(2,j,k)%fk) - &
                                                             work_for_tmp(2,j,k,4)
                        end do
                    end do
                end do
                !>
                !>**********************************************************************************************
                !> the min-distance on idim wall face
                !>
                i = idim
                !> face-center locations
                !>
                do k=2,kdim
                    do j=2,jdim
                        work_for_tmp(i,j,k,4) = 0.25e0*((iblock%coordinate(i,j  ,k  )%x   + &
                                                         iblock%coordinate(i,j+1,k  )%x   + &
                                                         iblock%coordinate(i,j  ,k+1)%x   + &
                                                         iblock%coordinate(i,j+1,k+1)%x     &
                                                         )*iblock%metric(i,j,k)%fi        + &
                                                        (iblock%coordinate(i,j  ,k  )%y   + &
                                                         iblock%coordinate(i,j+1,k  )%y   + &
                                                         iblock%coordinate(i,j  ,k+1)%y   + &
                                                         iblock%coordinate(i,j+1,k+1)%y     &
                                                         )*iblock%metric(i,j,k)%fj        + &
                                                        (iblock%coordinate(i,j  ,k  )%z   + &
                                                         iblock%coordinate(i,j+1,k  )%z   + &
                                                         iblock%coordinate(i,j  ,k+1)%z   + &
                                                         iblock%coordinate(i,j+1,k+1)%z     &
                                                         )*iblock%metric(i,j,k)%fk )
                    end do
                end do
                !>
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            iblock%turbulent(i,j,k)%distance_tur_2 = work_for_tmp(idim,j,k,4) - &
                                                            (work_for_tmp(i,j,k,1)*iblock%metric(idim,j,k)%fi   + &
                                                             work_for_tmp(i,j,k,2)*iblock%metric(idim,j,k)%fj   + &
                                                             work_for_tmp(i,j,k,3)*iblock%metric(idim,j,k)%fk )
                        end do
                    end do
                end do
                !>
                !>**********************************************************************************************
                !> the min-distance on j0 wall face
                j = 2
                !> face-center locations
                !>
                do k=2,kdim
                    do i=2,idim
                        work_for_tmp(i,j,k,4) = 0.25e0*((iblock%coordinate(i  ,j,k  )%x   + &
                                                         iblock%coordinate(i+1,j,k  )%x   + &
                                                         iblock%coordinate(i  ,j,k+1)%x   + &
                                                         iblock%coordinate(i+1,j,k+1)%x    &
                                                         )*iblock%metric(i,j,k)%gi  + &
                                                        (iblock%coordinate(i  ,j,k  )%y   + &
                                                         iblock%coordinate(i+1,j,k  )%y   + &
                                                         iblock%coordinate(i  ,j,k+1)%y   + &
                                                         iblock%coordinate(i+1,j,k+1)%y   &
                                                         )*iblock%metric(i,j,k)%gj  + &
                                                        (iblock%coordinate(i  ,j,k  )%z   + &
                                                         iblock%coordinate(i+1,j,k  )%z   + &
                                                         iblock%coordinate(i  ,j,k+1)%z   + &
                                                         iblock%coordinate(i+1,j,k+1)%z   )* &
                                                         iblock%metric(i,j,k)%gk )
                    end do
                end do
                !>
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            iblock%turbulent(i,j,k)%distance_tur_3 =(work_for_tmp(i,j,k,1)*iblock%metric(i,2,k)%gi + &
                                                                           work_for_tmp(i,j,k,2)*iblock%metric(i,2,k)%gj + &
                                                                           work_for_tmp(i,j,k,3)*iblock%metric(i,2,k)%gk ) - &
                                                                           work_for_tmp(i,2,k,4)
                        end do

                    end do
                end do
                !>
                !>**********************************************************************************************
                !> the min-distance on jdim wall face
                j = jdim
                !> face-center locations
                !>
                do k=2,kdim
                    do i=2,idim
                        work_for_tmp(i,j,k,4) = 0.25e0*((iblock%coordinate(i  ,j,k  )%x   + &
                                                         iblock%coordinate(i+1,j,k  )%x   + &
                                                         iblock%coordinate(i  ,j,k+1)%x   + &
                                                         iblock%coordinate(i+1,j,k+1)%x   )*&
                                                         iblock%metric(i,j,k)%gi  + &
                                                        (iblock%coordinate(i  ,j,k  )%y   + &
                                                         iblock%coordinate(i+1,j,k  )%y   + &
                                                         iblock%coordinate(i  ,j,k+1)%y   + &
                                                         iblock%coordinate(i+1,j,k+1)%y   )*&
                                                         iblock%metric(i,j,k)%gj  + &
                                                        (iblock%coordinate(i  ,j,k  )%z   + &
                                                         iblock%coordinate(i+1,j,k  )%z   + &
                                                         iblock%coordinate(i  ,j,k+1)%z   + &
                                                         iblock%coordinate(i+1,j,k+1)%z   )*&
                                                         iblock%metric(i,j,k)%gk  )
                    end do
                end do
                !>
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            iblock%turbulent(i,j,k)%distance_tur_4  = work_for_tmp(i,jdim,k,4) - &
                                              (work_for_tmp(i,j,k,1)*iblock%metric(i,jdim,k)%gi + &
                                               work_for_tmp(i,j,k,2)*iblock%metric(i,jdim,k)%gj + &
                                               work_for_tmp(i,j,k,3)*iblock%metric(i,jdim,k)%gk)
                        end do

                    end do
                end do
                !>
                !>**********************************************************************************************
                !> the min-distance on k0 wall face
                k = 2
                !> face-center locations
                !>
                do j=2,jdim
                    do i=2,idim
                        work_for_tmp(i,j,k,4) = 0.25e0*((iblock%coordinate(i  ,j  ,k)%x  + &
                                                         iblock%coordinate(i+1,j  ,k)%x  + &
                                                         iblock%coordinate(i  ,j+1,k)%x  + &
                                                         iblock%coordinate(i+1,j+1,k)%x &
                                                         )*iblock%metric(i,j,k)%hi  +   &
                                                        (iblock%coordinate(i  ,j  ,k)%y  + &
                                                         iblock%coordinate(i+1,j  ,k)%y  + &
                                                         iblock%coordinate(i  ,j+1,k)%y  + &
                                                         iblock%coordinate(i+1,j+1,k)%y &
                                                         )*iblock%metric(i,j,k)%hj  +   &
                                                        (iblock%coordinate(i  ,j  ,k)%z  + &
                                                         iblock%coordinate(i+1,j  ,k)%z  + &
                                                         iblock%coordinate(i  ,j+1,k)%z  + &
                                                         iblock%coordinate(i+1,j+1,k)%z &
                                                         )*iblock%metric(i,j,k)%hk )
                    end do
                end do
                !>
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            iblock%turbulent(i,j,k)%distance_tur_5  =( work_for_tmp(i,j,k,1)*iblock%metric(i,j,2)%hi + &
                                                                       work_for_tmp(i,j,k,2)*iblock%metric(i,j,2)%hj + &
                                                                       work_for_tmp(i,j,k,3)*iblock%metric(i,j,2)%hk ) - &
                                                                       work_for_tmp(i,j,2,4)
                        end do
                    end do
                end do
                !>
                !>**********************************************************************************************
                !> the min-distance on kdim wall face
                k = kdim
                !> face-center locations
                !>
                do j=2,jdim
                    do i=2,idim
                        work_for_tmp(i,j,k,4) = 0.25e0*((iblock%coordinate(i  ,j  ,k)%x   + &
                                                         iblock%coordinate(i+1,j  ,k)%x   + &
                                                         iblock%coordinate(i  ,j+1,k)%x   + &
                                                         iblock%coordinate(i+1,j+1,k)%x   )* &
                                                         iblock%metric(i,j,k)%hi + &
                                                        (iblock%coordinate(i  ,j  ,k)%y   + &
                                                         iblock%coordinate(i+1,j  ,k)%y   + &
                                                         iblock%coordinate(i  ,j+1,k)%y   + &
                                                         iblock%coordinate(i+1,j+1,k)%y   )* &
                                                         iblock%metric(i,j,k)%hj + &
                                                        (iblock%coordinate(i  ,j  ,k)%z   + &
                                                         iblock%coordinate(i+1,j  ,k)%z   + &
                                                         iblock%coordinate(i  ,j+1,k)%z   + &
                                                         iblock%coordinate(i+1,j+1,k)%z   )* &
                                                         iblock%metric(i,j,k)%hk )
                    end do
                end do
                !>
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            iblock%turbulent(i,j,k)%distance_tur_6  =  work_for_tmp(i,j,kdim,4) - &
                                                                     ( work_for_tmp(i,j,k ,1)*iblock%metric(i,j,kdim)%hi + &
                                                                       work_for_tmp(i,j,k ,2)*iblock%metric(i,j,kdim)%hj + &
                                                                       work_for_tmp(i,j,k ,3)*iblock%metric(i,j,kdim)%hk)
                        end do

                    end do
                end do
                !>
                deallocate(work_for_tmp)
#if defined PMPI
            end if
#endif
        end do
        !>
        end if

        return

	end subroutine
