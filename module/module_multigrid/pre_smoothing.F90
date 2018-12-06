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
    subroutine pre_smoothing(nbl)
        !>
        !>
        use  global_parameter
        use  blocks_module
        use  mesh_overlap_module
		!>
		!>
		implicit none
        real(kind = dprec) smoopi,smoopj,smoopk,t
        real(kind = dprec),dimension(:,:),allocatable::d
        !>
        !>
        integer i,j,k,nbl
        integer idim,jdim,kdim
        !>
        !>
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !> \brief
        !> the methods of smoothing is implicit residual smooting
        !> and form mozeyao's treatise the smooth operator in each level must less than three times,
        !> so the operator can divided into pre_smoothing and post_smoothing,
        !> but the pre_smoothing  and  post_smoothing have same operating,
        !> and we write the subroutine smooting to use irs to operating the rhs.

        !> smoopi residual smoothing coefficient for i  directions
        !> smoopj residual smoothing coefficient for j  directions
        !> smoopk residual smoothing coefficient for k  directions

        !> default values:0.3 0.3 0.3
        !> if the value is zero then not used the method of smooting
        !> also can compute this value by implicit residual smooting method
        !> we will coding  the irs method at ccfdv3.1 version
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        smoopi = smooth_r_coe!default 0.3
        smoopj = smooth_r_coe!default 0.3
        smoopk = smooth_r_coe!default 0.3
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        allocate(d(jdim+1,idim+kdim+2))
        !> smooting in i-direction
        !
        !
        !> computing the no-two-dimension case
        if(idim .gt. 1 .and. abs(smoopi) .gt. 0.0 )then
            do k=2,kdim
                do j=2,jdim
                    t = 1.0/(1.0+smoopi+smoopi)
                    d(j,2) = t*smoopi
                    iblock%variable(2,j,k)%res_1  = t*iblock%variable(2,j,k)%res_1
                    iblock%variable(2,j,k)%res_2  = t*iblock%variable(2,j,k)%res_2
                    iblock%variable(2,j,k)%res_3  = t*iblock%variable(2,j,k)%res_3
                    iblock%variable(2,j,k)%res_4  = t*iblock%variable(2,j,k)%res_4
                    iblock%variable(2,j,k)%res_5  = t*iblock%variable(2,j,k)%res_5
                end do
                do i=3,idim
                    do j=2,jdim
                        t           = 1.0/(1.0+smoopi+smoopi -smoopi*d(j,i-1))
                        d(j,i)      = t*smoopi
                        iblock%variable(i,j,k)%res_1 = t*(iblock%variable(i,j,k)%res_1  +smoopi*iblock%variable(i-1,j,k)%res_1)
                        iblock%variable(i,j,k)%res_2 = t*(iblock%variable(i,j,k)%res_2  +smoopi*iblock%variable(i-1,j,k)%res_2)
                        iblock%variable(i,j,k)%res_3 = t*(iblock%variable(i,j,k)%res_3  +smoopi*iblock%variable(i-1,j,k)%res_3)
                        iblock%variable(i,j,k)%res_4 = t*(iblock%variable(i,j,k)%res_4  +smoopi*iblock%variable(i-1,j,k)%res_4)
                        iblock%variable(i,j,k)%res_5 = t*(iblock%variable(i,j,k)%res_5  +smoopi*iblock%variable(i-1,j,k)%res_5)

                    end do
                end do
                do i=idim-1,2,-1
                    do j=2,jdim
                        iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1  + d(j,i)*iblock%variable(i+1,j,k)%res_1
                        iblock%variable(i,j,k)%res_2 = iblock%variable(i,j,k)%res_2  + d(j,i)*iblock%variable(i+1,j,k)%res_2
                        iblock%variable(i,j,k)%res_3 = iblock%variable(i,j,k)%res_3  + d(j,i)*iblock%variable(i+1,j,k)%res_3
                        iblock%variable(i,j,k)%res_4 = iblock%variable(i,j,k)%res_4  + d(j,i)*iblock%variable(i+1,j,k)%res_4
                        iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_5  + d(j,i)*iblock%variable(i+1,j,k)%res_5
                    end do
                end do
            end do
        end if
        !> \brief
        ! smooting in j-direction
        !
        !
        if(jdim .gt. 1 .and. abs(smoopj) .gt. 0.0 )then
            do i=2,idim
                do k=2,kdim
                    t = 1.0/(1.0+smoopj+smoopj)
                    d(2,k) = t*smoopj
                    iblock%variable(i,2,k)%res_1 = t*iblock%variable(i,2,k)%res_1
                    iblock%variable(i,2,k)%res_2 = t*iblock%variable(i,2,k)%res_2
                    iblock%variable(i,2,k)%res_3 = t*iblock%variable(i,2,k)%res_3
                    iblock%variable(i,2,k)%res_4 = t*iblock%variable(i,2,k)%res_4
                    iblock%variable(i,2,k)%res_5 = t*iblock%variable(i,2,k)%res_5
                end do
                do j=3,jdim
                    do k=2,kdim
                        t           = 1.0/(1.0+smoopj+smoopj -smoopj*d(j-1,k))
                        d(j,k)      = t*smoopj
                        iblock%variable(i,j,k)%res_1 = t*( iblock%variable(i,j,k)%res_1 + smoopj*iblock%variable(i,j-1,k)%res_1 )
                        iblock%variable(i,j,k)%res_2 = t*( iblock%variable(i,j,k)%res_2 + smoopj*iblock%variable(i,j-1,k)%res_2 )
                        iblock%variable(i,j,k)%res_3 = t*( iblock%variable(i,j,k)%res_3 + smoopj*iblock%variable(i,j-1,k)%res_3 )
                        iblock%variable(i,j,k)%res_4 = t*( iblock%variable(i,j,k)%res_4 + smoopj*iblock%variable(i,j-1,k)%res_4 )
                        iblock%variable(i,j,k)%res_5 = t*( iblock%variable(i,j,k)%res_5 + smoopj*iblock%variable(i,j-1,k)%res_5 )
                    end do
                end do
                do j=jdim-1,2,-1
                    do k=2,kdim
                        iblock%variable(i,j,k)%res_1  = iblock%variable(i,j,k)%res_1   +d(j,k)*iblock%variable(i,j+1,k)%res_1
                        iblock%variable(i,j,k)%res_2  = iblock%variable(i,j,k)%res_2   +d(j,k)*iblock%variable(i,j+1,k)%res_2
                        iblock%variable(i,j,k)%res_3  = iblock%variable(i,j,k)%res_3   +d(j,k)*iblock%variable(i,j+1,k)%res_3
                        iblock%variable(i,j,k)%res_4  = iblock%variable(i,j,k)%res_4   +d(j,k)*iblock%variable(i,j+1,k)%res_4
                        iblock%variable(i,j,k)%res_5  = iblock%variable(i,j,k)%res_5   +d(j,k)*iblock%variable(i,j+1,k)%res_5
                    end do
                end do
            end do
        end if
        !> \brief
        ! smooting in k-direction
        !
        !
        if( kdim .gt. 1 .and. abs(smoopk) .gt. 0.0 )then
            do i=2,idim
                do j=2,jdim 
                    t = 1.0/(1.0+smoopk+smoopk)
                    d(j,2) = t*smoopk
                    iblock%variable(i,j,2)%res_1  = t * iblock%variable(i,j,2)%res_1
                    iblock%variable(i,j,2)%res_2  = t * iblock%variable(i,j,2)%res_2
                    iblock%variable(i,j,2)%res_3  = t * iblock%variable(i,j,2)%res_3
                    iblock%variable(i,j,2)%res_4  = t * iblock%variable(i,j,2)%res_4
                    iblock%variable(i,j,2)%res_5  = t * iblock%variable(i,j,2)%res_5
                end do
                do k=3,kdim
                    do j=2,jdim
                        t           = 1.0/(1.0+ smoopk + smoopk -smoopk*d(j,k-1))
                        d(j,k)      = t*smoopk
                        iblock%variable(i,j,k)%res_1  = t * ( iblock%variable(i,j,k)%res_1   + smoopk*iblock%variable(i,j,k-1)%res_1)
                        iblock%variable(i,j,k)%res_2  = t * ( iblock%variable(i,j,k)%res_2   + smoopk*iblock%variable(i,j,k-1)%res_2)
                        iblock%variable(i,j,k)%res_3  = t * ( iblock%variable(i,j,k)%res_3   + smoopk*iblock%variable(i,j,k-1)%res_3)
                        iblock%variable(i,j,k)%res_4  = t * ( iblock%variable(i,j,k)%res_4   + smoopk*iblock%variable(i,j,k-1)%res_4)
                        iblock%variable(i,j,k)%res_5  = t * ( iblock%variable(i,j,k)%res_5   + smoopk*iblock%variable(i,j,k-1)%res_5)
                    end do
                end do
                do k=kdim-1,2,-1
                    do j=2,jdim
                        iblock%variable(i,j,k)%res_1  = iblock%variable(i,j,k)%res_1   + d(j,k)*iblock%variable(i,j,k+1)%res_1
                        iblock%variable(i,j,k)%res_2  = iblock%variable(i,j,k)%res_2   + d(j,k)*iblock%variable(i,j,k+1)%res_2
                        iblock%variable(i,j,k)%res_3  = iblock%variable(i,j,k)%res_3   + d(j,k)*iblock%variable(i,j,k+1)%res_3
                        iblock%variable(i,j,k)%res_4  = iblock%variable(i,j,k)%res_4   + d(j,k)*iblock%variable(i,j,k+1)%res_4
                        iblock%variable(i,j,k)%res_5  = iblock%variable(i,j,k)%res_5   + d(j,k)*iblock%variable(i,j,k+1)%res_5

                    end do
                end do
            end do
        end if
        deallocate(d)
        return
    end subroutine pre_smoothing

