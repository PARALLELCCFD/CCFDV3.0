!>***********************************************************************
!>     purpose:   To find the closest distance from a field point to the
!>                actual surface (i.e. not simply the closest discrete surface
!>                point), using local triangulation of the surface.
!>***********************************************************************
    subroutine triangles(dist,x,y,z,ax1,ay1,az1,s1,ax2,ay2,az2,s2,ax3,ay3,az3,s3)
        !>
        use global_parameter
        implicit none
        !>
        real(kind=dprec) dist,dist2,x,y,z
        real(kind=dprec) ax1,ay1,az1,s1
        real(kind=dprec) ax2,ay2,az2,s2
        real(kind=dprec) ax3,ay3,az3,s3
        real(kind=dprec) vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3
        real(kind=dprec) ab,ac,bc,aa,bb,cc
        real(kind=dprec) dd,s,t
        real(kind=dprec) hx,hy,hz
        !>
        !>
        if(s1 .ne. 0.0 .or. s2 .ne. 0.0 .or. s3 .ne. 0.0 ) return
        !>
        !> vector a
        vx1 = ax2 - ax1
        vy1 = ay2 - ay1
        vz1 = az2 - az1
        !> vector b
        vx2 = ax3 - ax1
        vy2 = ay3 - ay1
        vz2 = az3 - az1
		!> vector c
		vx3 = x   - ax1
		vy3 = y   - ay1
		vz3 = z   - az1
		!>
		aa = vx1**2 + vy1**2 + vz1**2
		bb = vx2**2 + vy2**2 + vz2**2
		ab = vx1*vx2 + vy1*vy2 + vz1*vz2
		!>
		!>
		dd = ab**2 - aa*bb
		!>
		!>
		if(dd .ne. 0.0)then
            ac = vx1*vx3 + vy1*vy3 + vz1*vz3
            bc = vx2*vx3 + vy2*vy3 + vz2*vz3
            s  = (ab*bc - bb*ac )/dd
            t  = (ab*ac - aa*bc )/dd
            if(s .ge. 0 .and. t .ge. 0 .and. (t+s) .le. 1)then
                hx = vx3 - (s*vx1 + t*vx2)
                hy = vy3 - (s*vy1 + t*vy2)
                hz = vz3 - (s*vz1 + t*vz2)
                dist = min(dist,sqrt(hx**2 + hy**2 + hz**2))
                !>
                return
            end if
        end if
        !>
        !>
        dist2 = dist**2
        !>
        !>
        if(aa .ne. 0)then
            ac = vx1*vx3 + vy1*vy3 + vz1*vz3
            t  = ac / aa
            if(t .ge. 0 .and. t .le. 1)then
                hx = vx3- t*vx1
                hy = vy3- t*vy1
                hz = vz3- t*vz1
                dist2 = min(dist2,(hx**2 + hy**2 + hz**2))
                !>
            end if
        end if
        !>
        !>
        if( bb .ne. 0 )then
            bc = vx2*vx3 + vy2*vy3 + vz2*vz3
            t  = bc / bb
            if(t .ge. 0 .and. t .le. 1)then
                hx = vx3- t*vx2
                hy = vy3- t*vy2
                hz = vz3- t*vz2
                dist2 = min(dist2,(hx**2 + hy**2 + hz**2))
                !>
            end if
        end if
        !> vector c
        vx3 = x - ax2
        vy3 = y - ay2
        vz3 = z - az2
        !> vector a
        vx1 = ax3 - ax2
        vy1 = ay3 - ay2
        vz1 = az3 - az2
        !>
        aa = vx1**2 + vy1**2 + vz1**2
        if(aa .ne. 0)then
            ac = vx1*vx3 + vy1*vy3 + vz1*vz3
            t = ac / aa
            if(t .ge. 0 .and. t .le. 1)then
                hx = vx3- t*vx1
                hy = vy3- t*vy1
                hz = vz3- t*vz1
                dist2 = min(dist2,(hx**2 + hy**2 + hz**2))
                !>
            end if
        end if
        !>
        !>
        if(dist2 .lt. (dist*dist)) dist = sqrt(dist2)
        !>
        !>
        return
    end subroutine
