!>***********************************************************************
!>     purpose:   Heap sort methods
!>***********************************************************************
    subroutine heap_sort(n,a,indx)
        !>
        use global_parameter
        implicit none
        !>
        real(kind = dprec) :: a(1:n)
        real(kind = dprec) :: temp
        integer :: indx(1:n)
        integer :: n,l1,r1,i,j,k,ip
        !>
        !>

        if(n .lt. 2) return
        !>
        !>
        l1 = n/2 + 1
        r1 = n
        10   continue
        if(l1 .gt. 1)then
            l1    = l1-1
            ip    = indx(l1)
            temp  = a(ip)
        else
            ip    = indx(r1)
            temp     = a(ip)
            indx(r1) = indx(1)
            r1       = r1 - 1
            if(r1 .eq. 1)then
                indx(1) = ip
                return
            endif
        endif
        i = l1
        j = l1+l1
        20  continue
        if(j .le. r1)then
            if(j .lt. r1)then
                if(real(a(indx(j))) .lt. real(a(indx(j+1)))) j=j+1
            endif
            if(real(temp) .lt. real(a(indx(j))))then
                indx(i) = indx(j)
                i       = j
                j       = 2*j
            else
                j       = r1 + 1
            endif
            go to 20
        endif
        indx(i) = ip
        go to 10

		!>
        return
    end subroutine
