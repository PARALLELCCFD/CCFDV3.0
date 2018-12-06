!>***********************************************************************
!>     purpose:   Shell sort methods
!>***********************************************************************
    subroutine shell_sort(n,a,indx)
        !>
        use global_parameter
        implicit none
        !>
        real(kind = dprec) :: a(1:n)
        integer :: indx(1:n)
        !>
        !>
        integer :: i,j,gap,inx,n
        !>
        !>
        gap = n / 2
        !>
        !>
!        write(*,'("shells=",i8)') n
        do while(gap .ge. 1)
            do i=gap+1,n
                j=i-gap
#if 1
                do while(j .gt. 0 )
                    if(a(indx(j)) .gt. a(indx(j+gap)))then
                        inx         = indx(j)
                        indx(j)     = indx(j+gap)
                        indx(j+gap) = inx
                        j = j - gap
                    else
                        exit
                    end if
                end do
#else
                100 continue
                if(j .le. 0 ) cycle
                if(a(indx(j)) .le. a(indx(j+gap))) cycle
                inx         = indx(j)
                indx(j)     = indx(j+gap)
                indx(j+gap) = inx
                j = j - gap
                goto 100
#endif
                !>
            end do
            !>
            gap = gap / 2
        end do
!        do i=1,n
!            write(180720,'("shell_sort",i6,e24.16)') indx(i),a(indx(i))
!        end do
        !>
		!>
        return
    end subroutine
