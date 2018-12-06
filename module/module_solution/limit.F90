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
    subroutine limit(point,max_ijk,rr,pp,du,dv,dw,dp,dr,l1,l2,l3,l4,l5,r1,r2,r3,r4,r5)
		!>
		!>
        use global_parameter
        use nodes_paras
        !>
		implicit none
        !>
        integer::i,j,k,&
                 point,max_ijk
        real(kind=dprec)::du(0:max_ijk),dv(0:max_ijk),dw(0:max_ijk),dp(0:max_ijk),dr(0:max_ijk),&
                          rr(1:max_ijk),pp(1:max_ijk),&
                          dp1(0:max_ijk),dr1(0:max_ijk),&
                          dp2(0:max_ijk),dr2(0:max_ijk),&
                          du1(0:max_ijk),du2(0:max_ijk),dv1(0:max_ijk),dv2(0:max_ijk),dw1(0:max_ijk),dw2(0:max_ijk),&
                          l1(max_ijk),l2(max_ijk),l3(max_ijk),l4(max_ijk),l5(max_ijk),&
                          r1(max_ijk),r2(max_ijk),r3(max_ijk),r4(max_ijk),r5(max_ijk)
        real(kind=dprec)::eps,delx,aravg1,apavg1,aravg2,apavg2,logical_depending
        real(kind=dprec)::x1,x2,x3,x4,x5,lim_minmod,xterm
        real(kind=dprec)::t1,t2,t3,t4,phi,x6,x7




        eps  = 1.e-6
        !>kap value
        !>xkap= -1 upwind schemes
        !>xkap= 0 fromm schemes
        !>xkap= 1/3 at 2-d problem ,the schemes have thrid,and at 3-d problem ,the schemes just have second
        !>xkap= 1
        !xkap = -1.0!1./3. !
        phi  = (3.0-xkap)/(1.0-xkap+eps)
        !>
        if( xkap .lt. -1.0) then
            !>first order
            !>
            do i=1,point
                l1(i)  =0.
                l2(i)  =0.
                l3(i)  =0.
                l4(i)  =0.
                l5(i)  =0.
                r1(i)  =0.
                r2(i)  =0.
                r3(i)  =0.
                r4(i)  =0.
                r5(i)  =0.
            end do
            return
        end if
        !> higher order
        !>
        !> gradients by average values for pressure and density
        do i=2,point-1
            !>in pratice,the gradients of density are biased by
            !> an average value in order to improve the robustness
            !> of the calculation an higher mach number and in the
            !>early transient stages of a solution
            dr1(i)    = dr(i)  / (rr(i)  - 0.5* dr(i)  )
            dr2(i+1)  = dr(i+1)/ (rr(i)  + 0.5* dr(i+1))
            !>in pratice,the gradients of pressure are biased by
            !> an average value in order to improve the robustness
            !> of the calculation an higher mach number and in the
            !>early transient stages of a solution
            dp1(i)    = dp(i)  / (pp(i)  - 0.5* dp(i)  )
            dp2(i+1)  = dp(i+1)/ (pp(i)  + 0.5* dp(i+1))

            !> i-dir u
            du1(i)   = du(i)
            du2(i+1) = du(i+1)
            !> j-dir v
            dv1(i)   = dv(i)
            dv2(i+1) = dv(i+1)
            !> k-dir w
            dw1(i)   = dw(i)
            dw2(i+1) = dw(i+1)
        end do
        !>
        !> left boundary values
        dr1(1)  = dr(0)  / (rr(1)  - 0.5* dr(0))
        dr2(2)  = dr(1)  / (rr(1)  + 0.5* dr(1))
        !>
        dp1(1)  = dp(0)  / (pp(1)  - 0.5* dp(0))
        dp2(2)  = dp(1)  / (pp(1)  + 0.5* dp(1))
        !>
        du1(1)  = du(0)
        du2(2)  = du(1)
        !>
        dv1(1)  = dv(0)
        dv2(2)  = dv(1)
        !>
        dw1(1)  = dw(0)
        dw2(2)  = dw(1)
        !>
        !>right boundary values
        dr1(point)    = dr(point+1)  / (rr(point)  - 0.5* dr(point+1))
        dr2(point+1)  = dr(point+2)  / (rr(point)  + 0.5* dr(point+2))
        !>
        dp1(point)    = dp(point+1)  / (pp(point)  - 0.5* dp(point+1))
        dp2(point+1)  = dp(point+2)  / (pp(point)  + 0.5* dp(point+2))
        !>
        du1(point)    = du(point+1)
        du2(point+1)  = du(point+2)
        !>
        dv1(point)    = dv(point+1)
        dv2(point+1)  = dv(point+2)
        !>
        dw1(point)    = dw(point+1)
        dw2(point+1)  = dw(point+2)
        !> od - second order
        !>
        if(limiter .eq. -1 )then
    !        xterm = 0.25*xkap
            do i=1,point
                !>density of left
                r1(i) = logical_depending(0.5*dr1(i),dr2(i+1),(abs(dr1(i) ).lt. abs(dr2(i+1))) )
                l1(i) = r1(i)
                !>u of left
                r2(i) = logical_depending(0.5*du1(i),du2(i+1),(abs(du1(i) ).lt. abs(du2(i+1))) )
                l2(i) = r2(i)
                !>v of left
                r3(i) = logical_depending(0.5*dv1(i),dv2(i+1),(abs(dv1(i) ).lt. abs(dv2(i+1))) )
                l3(i) = r3(i)
                !>w of left
                r4(i) = logical_depending(0.5*dw1(i),du2(i+1),(abs(dw1(i) ).lt. abs(dw2(i+1))) )
                l4(i) = r4(i)
                !>pressuer of left
                r5(i) = logical_depending(0.5*dp1(i),dp2(i+1),(abs(dp1(i) ).lt. abs(dp2(i+1))) )
                l5(i) = r5(i)
            end do
        end if
        !>
        !> without the limiter
        !> kappa scheme and xkap =-1.0
        !> full upwind scheme
        !>
        if(limiter .eq. 0 )then
    !        xterm = 0.25*xkap
            do i=1,point
                !>density of left
                r1(i) = 0.25*((1.0-xkap)*dr1(i)   + (1.0+xkap)*dr2(i+1))
                !>density of right
                l1(i) = 0.25*((1.0-xkap)*dr2(i+1) + (1.0+xkap)*dr1(i))
                !>u of left
                r2(i) = 0.25*((1.0-xkap)*du1(i)    + (1.0+xkap)*du2(i+1))
                !>u of right
                l2(i) = 0.25*((1.0-xkap)*du2(i+1)  + (1.0+xkap)*du1(i))
                !>v of left
                r3(i) = 0.25*((1.0-xkap)*dv1(i)    + (1.0+xkap)*dv2(i+1))
                !>v of right
                l3(i) = 0.25*((1.0-xkap)*dv2(i+1)  + (1.0+xkap)*dv1(i))
                !>w of left
                r4(i) = 0.25*((1.0-xkap)*dw1(i)    + (1.0+xkap)*dw2(i+1))
                !>w of right
                l4(i) = 0.25*((1.0-xkap)*dw2(i+1)  + (1.0+xkap)*dw1(i))
                !>pressuer of left
                r5(i) = 0.25*((1.0-xkap)*dp1(i)    + (1.0+xkap)*dp2(i+1))
                !>pressuer of right
                l5(i) = 0.25*((1.0-xkap)*dp2(i+1)  + (1.0+xkap)*dp1(i))
            end do
        end if
        !>
        !>
        !> van albada's scheme limiter
        !>
        !>
        if(limiter .eq. 1)then
            do i=1,point
                !>
                !>density of left
                xterm = (2.0*dr2(i+1)*dr1(i)+eps)/(dr2(i+1)*dr2(i+1)+dr1(i)*dr1(i)+eps)
                !>
                r1(i) = 0.25*xterm*((1.-xkap*xterm)*dr1(i)   + (1.+xkap*xterm)*dr2(i+1))
                !>
                !>density of right
                l1(i) = 0.25*xterm*((1.-xkap*xterm)*dr2(i+1) + (1.+xkap*xterm)*dr1(i))
                !>
                !>u of left
                xterm = (2.0*du2(i+1)*du1(i)+eps)/(du2(i+1)*du2(i+1)+du1(i)*du1(i)+eps)
                !>
                r2(i) = 0.25*xterm*((1.-xkap*xterm)*du1(i)    + (1.+xkap*xterm)*du2(i+1))
                !>
                !>u of right
                l2(i) = 0.25*xterm*((1.-xkap*xterm)*du2(i+1)  + (1.+xkap*xterm)*du1(i))
                !>
                !>v of left
                xterm = (2.0*dv2(i+1)*dv1(i)+eps)/(dv2(i+1)*dv2(i+1)+dv1(i)*dv1(i)+eps)
                !>
                r3(i) = 0.25*xterm*((1.-xkap*xterm)*dv1(i)    + (1.+xkap*xterm)*dv2(i+1))
                !>
                !>v of right
                l3(i) = 0.25*xterm*((1.-xkap*xterm)*dv2(i+1)  + (1.+xkap*xterm)*dv1(i))
                !>
                !>w of left
                xterm = (2.0*dw2(i+1)*dw1(i)+eps)/(dw2(i+1)*dw2(i+1)+dw1(i)*dw1(i)+eps)
                !>
                r4(i) = 0.25*xterm*((1.-xkap*xterm)*dw1(i)    + (1.+xkap*xterm)*dw2(i+1))
                !>
                !>w of right
                l4(i) = 0.25*xterm*((1.-xkap*xterm)*dw2(i+1)  + (1.+xkap*xterm)*dw1(i))
                !>
                !>pressuer of left
                xterm = (2.0*dp2(i+1)*dp1(i)+eps)/(dp2(i+1)*dp2(i+1)+dp1(i)*dp1(i)+eps)
                !>
                r5(i) = 0.25*xterm*((1.-xkap*xterm)*dp1(i)   + (1.+xkap*xterm)*dp2(i+1))
                !>
                !>pressuer of right
                l5(i) = 0.25*xterm*((1.-xkap*xterm)*dp2(i+1) + (1.+xkap*xterm)*dp1(i))
            end do
        end if
        !>
        !>
        !>min-mod's scheme limiter
        !>
        !>
        !>

        if(limiter .eq. 2) then

            !>min-mod limiter
            do i=1,point
#if 1
                !>rho computing
                x4       = dr2(i+1)*dr1(i)
                if(real(x4) .lt. 0.0)then
                    x6 = 0.0
                    x7 = 0.0
                else
                    x6 = dr1(i)
                    x7 = dr2(i+1)
                end if
!                x6       = logical_depending(0.0,dr1(i),(real(x4) .lt. 0.0))
!                x7       = logical_depending(0.0,dr2(i+1),(real(x4) .lt. 0.0))
                x4       = abs(phi*x6)
                x3       = abs(phi*x7)
                x5       = abs(x7)
                if(real(x4) .lt. real(x5))then
                    x7 = phi*x6
                end if
!                x7       = logical_depending(phi*x6,x7,(real(x4) .lt. real(x5)))
                x5       = abs(x6)
                if(real(x3) .lt. real(x5))then
                    x6 = phi*x7
                end if
!                x6       = logical_depending(phi*x7,x6,(real(x3) .lt. real(x5)))
                x3       = 0.25*(x6 + x7)
                x4       = xkap*0.25*(x7 - x6)
                l1(i)    = x3-x4
                r1(i)    = x3+x4
                !> u computing
                x4       = du2(i+1)*du1(i)
                if(real(x4) .lt. 0.0)then
                    x6 = 0.0
                    x7 = 0.0
                else
                    x6 = du1(i)
                    x7 = du2(i+1)
                end if
!                x6       = logical_depending(0.0,du1(i),(real(x4) .lt. 0.0))
!                x7       = logical_depending(0.0,du2(i+1),(real(x4) .lt. 0.0))
                x4       = abs(phi*x6)
                x3       = abs(phi*x7)
                x5       = abs(x7)
                x7       = logical_depending(phi*x6,x7,(real(x4) .lt. real(x5)))
                x5       = abs(x6)
                x6       = logical_depending(phi*x7,x6,(real(x3) .lt. real(x5)))
                x3       = 0.25*(x6 + x7)
                x4       = xkap*0.25*(x7 - x6)
                l2(i)    = x3-x4
                r2(i)    = x3+x4
                !> v computing
                x4       = dv2(i+1)*dv1(i)
                if(real(x4) .lt. 0.0)then
                    x6 = 0.0
                    x7 = 0.0
                else
                    x6 = dv1(i)
                    x7 = dv2(i+1)
                end if
!                x6       = logical_depending(0.0,dv1(i),(real(x4) .lt. 0.0))
!                x7       = logical_depending(0.0,dv2(i+1),(real(x4) .lt. 0.0))
                x4       = abs(phi*x6)
                x3       = abs(phi*x7)
                x5       = abs(x7)
                if(real(x4) .lt. real(x5))then
                    x7 = phi*x6
                end if
!                x7       = logical_depending(phi*x6,x7,(real(x4) .lt. real(x5)))
                x5       = abs(x6)
                if(real(x3) .lt. real(x5))then
                    x6 = phi*x7
                end if
!                x6       = logical_depending(phi*x7,x6,(real(x3) .lt. real(x5)))
                x3       = 0.25*(x6 + x7)
                x4       = xkap*0.25*(x7 - x6)
                l3(i)    = x3-x4
                r3(i)    = x3+x4
                !> w computing
                x4       = dw2(i+1)*dw1(i)
                if(real(x4) .lt. 0.0)then
                    x6 = 0.0
                    x7 = 0.0
                else
                    x6 = dw1(i)
                    x7 = dw2(i+1)
                end if
!                x6       = logical_depending(0.0,dw1(i),(real(x4) .lt. 0.0))
!                x7       = logical_depending(0.0,dw2(i+1),(real(x4) .lt. 0.0))
                x4       = abs(phi*x6)
                x3       = abs(phi*x7)
                x5       = abs(x7)
                if(real(x4) .lt. real(x5))then
                    x7 = phi*x6
                end if
!                x7       = logical_depending(phi*x6,x7,(real(x4) .lt. real(x5)))
                x5       = abs(x6)
                if(real(x3) .lt. real(x5))then
                    x6 = phi*x7
                end if
!                x6       = logical_depending(phi*x7,x6,(real(x3) .lt. real(x5)))
                x3       = 0.25*(x6 + x7)
                x4       = xkap*0.25*(x7 - x6)
                l4(i)    = x3-x4
                r4(i)    = x3+x4
                !>pressure computing
                x4       = dp2(i+1)*dp1(i)
                if(real(x4) .lt. 0.0)then
                    x6 = 0.0
                    x7 = 0.0
                else
                    x6 = dp1(i)
                    x7 = dp2(i+1)
                end if
!                x6       = logical_depending(0.0,dp1(i),(real(x4) .lt. 0.0))
!                x7       = logical_depending(0.0,dp2(i+1),(real(x4) .lt. 0.0))
                x4       = abs(phi*x6)
                x3       = abs(phi*x7)
                x5       = abs(x7)
                if(real(x4) .lt. real(x5))then
                    x7 = phi*x6
                end if
!                x7       = logical_depending(phi*x6,x7,(real(x4) .lt. real(x5)))
                x5       = abs(x6)
                if(real(x3) .lt. real(x5))then
                    x6 = phi*x7
                end if
!                x6       = logical_depending(phi*x7,x6,(real(x3) .lt. real(x5)))
                x3       = 0.25*(x6 + x7)
                x4       = xkap*0.25*(x7 - x6)
                l5(i)    = x3-x4
                r5(i)    = x3+x4
#else
                !>
                !>
                !>rho
                !>left dleta rho in the cneter cell i+1/2
                r1(i) = 0.25*((1.- xkap)*lim_minmod(dr1(i),dr2(i+1),phi)+(1.+ xkap)*lim_minmod(dr2(i+1),dr1(i),phi))
                !>right dleta rho in the cneter cell i+1/2
                l1(i) = 0.25*((1.- xkap)*lim_minmod(dr2(i+1),dr1(i),phi)+(1.+ xkap)*lim_minmod(dr1(i),dr2(i+1),phi))
                !>
                !>u
                !>left dleta u in the cneter cell i+1/2
                r2(i) = 0.25*((1.- xkap)*lim_minmod(du1(i),du2(i+1),phi)+(1.+ xkap)*lim_minmod(du2(i+1),du1(i),phi))
                !>right dleta u in the cneter cell i+1/2
                l2(i) = 0.25*((1.- xkap)*lim_minmod(du2(i+1),du1(i),phi)+(1.+ xkap)*lim_minmod(du1(i),du2(i+1),phi))
                !>v
                !>left dleta v in the cneter cell i+1/2
                r3(i) = 0.25*((1.- xkap)*lim_minmod(dv1(i),dv2(i+1),phi)+(1.+ xkap)*lim_minmod(dv2(i+1),dv1(i),phi))
                !>right dleta v in the cneter cell i+1/2
                l3(i) = 0.25*((1.- xkap)*lim_minmod(dv2(i+1),dv1(i),phi)+(1.+ xkap)*lim_minmod(dv1(i),dv2(i+1),phi))
                !>
                !>w
                !>left dleta w in the cneter cell i+1/2
                r4(i) = 0.25*((1.- xkap)*lim_minmod(dw1(i),dw2(i+1),phi)+(1.+ xkap)*lim_minmod(dw2(i+1),dw1(i),phi))
                !>right dleta w in the cneter cell i+1/2
                l4(i) = 0.25*((1.- xkap)*lim_minmod(dw2(i+1),dw1(i),phi)+(1.+ xkap)*lim_minmod(dw1(i),dw2(i+1),phi))
                !>
                !>pressure
                !>left dleta pressure in the cneter cell i+1/2
                r5(i) = 0.25*((1.- xkap)*lim_minmod(dp1(i),dp2(i+1),phi)+(1.+ xkap)*lim_minmod(dp2(i+1),dp1(i),phi))
                !>right dleta pressure in the cneter cell i+1/2
                l5(i) = 0.25*((1.- xkap)*lim_minmod(dp2(i+1),dp1(i),phi)+(1.+ xkap)*lim_minmod(dp2(i),dp1(i+1),phi))
                !>
#endif
            end do
        end if
        !>
        !>
        !>spekreijse vankat's scheme limiter
        !>and turn xkap = 1/3.
        !>
        !>
        if(limiter .eq. 3) then
            !>
            delx = 10.0/float(point-2)
            eps  = delx**3
            !>
            do i=1,point
                !>
                !> density computing
                t1    = dr1(i)**2
                t2    = dr2(i+1)**2
                t3    = dr1(i)*dr2(i+1)
                t4    = dr1(i)+dr2(i+1)
                xterm = 0.5*(t3+eps)/(2.0*(t1+t2)-t3+3.0*eps)
                r1(i) = (dr2(i+1)+t4)*xterm
                l1(i) = (dr1(i)  +t4)*xterm
                !>
                !> u computing
                t1    = du1(i)**2
                t2    = du2(i+1)**2
                t3    = du1(i)*du2(i+1)
                t4    = du1(i)+du2(i+1)
                xterm = 0.5*(t3+eps)/(2.0*(t1+t2)-t3+3.0*eps)
                r2(i) = (du2(i+1)+t4)*xterm
                l2(i) = (du1(i)  +t4)*xterm
                !>
                !> v computing
                t1    = dv1(i)**2
                t2    = dv2(i+1)**2
                t3    = dv1(i)*dv2(i+1)
                t4    = dv1(i)+dv2(i+1)
                xterm = 0.5*(t3+eps)/(2.0*(t1+t2)-t3+3.0*eps)
                r3(i) = (dv2(i+1)+t4)*xterm
                l3(i) = (dv1(i)  +t4)*xterm
                !>
                !> w computing
                t1    = dw1(i)**2
                t2    = dw2(i+1)**2
                t3    = dw1(i)*dw2(i+1)
                t4    = dw1(i)+dw2(i+1)
                xterm = 0.5*(t3+eps)/(2.0*(t1+t2)-t3+3.0*eps)
                r4(i) = (dw2(i+1)+t4)*xterm
                l4(i) = (dw1(i)  +t4)*xterm
                !>
                !> pressure computing
                t1    = dp1(i)**2
                t2    = dp2(i+1)**2
                t3    = dp1(i)*dp2(i+1)
                t4    = dp1(i)+dp2(i+1)
                xterm = 0.5*(t3+eps)/(2.0*(t1+t2)-t3+3.0*eps)
                r5(i) = (dp2(i+1)+t4)*xterm
                l5(i) = (dp1(i)  +t4)*xterm
    !            write(222,'(10e18.11)')l1(i),r1(i),l2(i),r2(i),l3(i),r3(i),l4(i),r4(i),l5(i),r5(i)
                !>
                !>
            end do
        end if
        !>
        !> gradients by average values for pressure and density
        do i=1,point
            !>gradients by average values-density
            l1(i) =  rr(i)*l1(i)
            r1(i) =  rr(i)*r1(i)
            !>gradients by average values-pressure
            l5(i) =  pp(i)*l5(i)
            r5(i) =  pp(i)*r5(i)
        end do

!        return
    end  subroutine limit


