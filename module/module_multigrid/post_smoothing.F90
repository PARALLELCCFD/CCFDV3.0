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
    subroutine post_smoothing(idim,jdim,kdim,deltaq)
		!>
		!>
        use global_parameter
        !>
        !>
		implicit none
        !>
        real(kind = dprec) smoopi,smoopj,smoopk,t
        real(kind = dprec),dimension(:,:),allocatable :: d
        real(kind = dprec):: deltaq(2:idim,2:jdim,2:kdim,1:7)
        integer i,j,k
        integer idim,jdim,kdim
        allocate(d(jdim+1,idim+kdim+2))
        !> \brief
        !> the methods of smoothing is implicit residual smooting at prolongated procedure
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

        smoopi = smooth_c_coe!0.3
        smoopj = smooth_c_coe!0.3
        smoopk = smooth_c_coe!0.3

        !allocate(d(jdim,idim+kdim))
        !> smooting in i-direction
        !
        !
        !> computing the three-dimension
        if(idim .gt. 1 .and. abs(smoopi) .gt. 0.0 )then
            do k=2,kdim
                do j=2,jdim
                    t = 1./(1.0+smoopi+smoopi)
                    d(j,2) = t*smoopi
                    deltaq(2,j,k,1) = t*deltaq(2,j,k,1)
                    deltaq(2,j,k,2) = t*deltaq(2,j,k,2)
                    deltaq(2,j,k,3) = t*deltaq(2,j,k,3)
                    deltaq(2,j,k,4) = t*deltaq(2,j,k,4)
                    deltaq(2,j,k,5) = t*deltaq(2,j,k,5)
                    !>
                    !>
                end do
                do i=3,idim
                    do j=2,jdim
                        t           = 1.0/(1.0+smoopi+smoopi -smoopi*d(j,i-1))
                        d(j,i)      = t*smoopi
                        deltaq(i,j,k,1) = t*(deltaq(i,j,k,1)  +smoopi*deltaq(i-1,j,k,1))
                        deltaq(i,j,k,2) = t*(deltaq(i,j,k,2)  +smoopi*deltaq(i-1,j,k,2))
                        deltaq(i,j,k,3) = t*(deltaq(i,j,k,3)  +smoopi*deltaq(i-1,j,k,3))
                        deltaq(i,j,k,4) = t*(deltaq(i,j,k,4)  +smoopi*deltaq(i-1,j,k,4))
                        deltaq(i,j,k,5) = t*(deltaq(i,j,k,5)  +smoopi*deltaq(i-1,j,k,5))
                    end do
                end do
                do i=idim-1,2,-1
                    do j=2,jdim
                        deltaq(i,j,k,1) = deltaq(i,j,k,1)  +d(j,i)*deltaq(i+1,j,k,1)
                        deltaq(i,j,k,2) = deltaq(i,j,k,2)  +d(j,i)*deltaq(i+1,j,k,2)
                        deltaq(i,j,k,3) = deltaq(i,j,k,3)  +d(j,i)*deltaq(i+1,j,k,3)
                        deltaq(i,j,k,4) = deltaq(i,j,k,4)  +d(j,i)*deltaq(i+1,j,k,4)
                        deltaq(i,j,k,5) = deltaq(i,j,k,5)  +d(j,i)*deltaq(i+1,j,k,5)
                        !>
                        !>
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
                    deltaq(i,2,k,1) = t*deltaq(i,2,k,1)
                    deltaq(i,2,k,2) = t*deltaq(i,2,k,2)
                    deltaq(i,2,k,3) = t*deltaq(i,2,k,3)
                    deltaq(i,2,k,4) = t*deltaq(i,2,k,4)
                    deltaq(i,2,k,5) = t*deltaq(i,2,k,5)
                    !>
                    !>
                end do
                do j=3,jdim
                    do k=2,kdim
                        t           = 1.0/(1.0+smoopj+smoopj -smoopj*d(j-1,k))
                        d(j,k)      = t*smoopj
                        deltaq(i,j,k,1) = t*(deltaq(i,j,k,1)  +smoopj*deltaq(i,j-1,k,1))
                        deltaq(i,j,k,2) = t*(deltaq(i,j,k,2)  +smoopj*deltaq(i,j-1,k,2))
                        deltaq(i,j,k,3) = t*(deltaq(i,j,k,3)  +smoopj*deltaq(i,j-1,k,3))
                        deltaq(i,j,k,4) = t*(deltaq(i,j,k,4)  +smoopj*deltaq(i,j-1,k,4))
                        deltaq(i,j,k,5) = t*(deltaq(i,j,k,5)  +smoopj*deltaq(i,j-1,k,5))
                        !>
                        !>
                    end do
                end do
                do j=jdim-1,2,-1
                    do k=2,kdim
                        deltaq(i,j,k,1) = deltaq(i,j,k,1)  +d(j,k)*deltaq(i,j+1,k,1)
                        deltaq(i,j,k,2) = deltaq(i,j,k,2)  +d(j,k)*deltaq(i,j+1,k,2)
                        deltaq(i,j,k,3) = deltaq(i,j,k,3)  +d(j,k)*deltaq(i,j+1,k,3)
                        deltaq(i,j,k,4) = deltaq(i,j,k,4)  +d(j,k)*deltaq(i,j+1,k,4)
                        deltaq(i,j,k,5) = deltaq(i,j,k,5)  +d(j,k)*deltaq(i,j+1,k,5)
                        !>
                        !>
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
                    deltaq(i,j,2,1) = t * deltaq(i,j,2,1)
                    deltaq(i,j,2,2) = t * deltaq(i,j,2,2)
                    deltaq(i,j,2,3) = t * deltaq(i,j,2,3)
                    deltaq(i,j,2,4) = t * deltaq(i,j,2,4)
                    deltaq(i,j,2,5) = t * deltaq(i,j,2,5)
                    !>
                    !>
                end do
                do k=3,kdim
                    do j=2,jdim
                        t           = 1.0/(1.0+smoopk+smoopk -smoopk*d(j,k-1))
                        d(j,k)      = t*smoopk
                        deltaq(i,j,k,1) = t * (deltaq(i,j,k,1)  + smoopk*deltaq(i,j,k-1,1))
                        deltaq(i,j,k,2) = t * (deltaq(i,j,k,2)  + smoopk*deltaq(i,j,k-1,2))
                        deltaq(i,j,k,3) = t * (deltaq(i,j,k,3)  + smoopk*deltaq(i,j,k-1,3))
                        deltaq(i,j,k,4) = t * (deltaq(i,j,k,4)  + smoopk*deltaq(i,j,k-1,4))
                        deltaq(i,j,k,5) = t * (deltaq(i,j,k,5)  + smoopk*deltaq(i,j,k-1,5))
                        !>
                        !>
                    end do
                end do
                do k=kdim-1,2,-1
                    do j=2,jdim
                        deltaq(i,j,k,1) = deltaq(i,j,k,1)  + d(j,k)*deltaq(i,j,k+1,1)
                        deltaq(i,j,k,2) = deltaq(i,j,k,2)  + d(j,k)*deltaq(i,j,k+1,2)
                        deltaq(i,j,k,3) = deltaq(i,j,k,3)  + d(j,k)*deltaq(i,j,k+1,3)
                        deltaq(i,j,k,4) = deltaq(i,j,k,4)  + d(j,k)*deltaq(i,j,k+1,4)
                        deltaq(i,j,k,5) = deltaq(i,j,k,5)  + d(j,k)*deltaq(i,j,k+1,5)
                        !>
                        !>
                    end do
                end do
            end do
        end if
        deallocate(d)
        !>
        !>
        !>
        return
    end subroutine post_smoothing
