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
    subroutine roe(ax,ay,az,aa,rl,ul,vl,wl,pl,rr,ur,vr,wr,pr,f1,f2,f3,f4,f5)
		!>
		!>
        use global_parameter
		implicit none
		!>
		!>
        real(kind =dprec):: allspeedroe,m1,m2,m3,lmach
        real(kind =dprec):: epsaa,epsbb,epscc
        real(kind =dprec):: gx,gy,gz
        real(kind =dprec):: zero,gagm
        real(kind =dprec):: ax,ay,az,aa,&
                            rl,ul,vl,wl,pl,hl,&
                            rr,ur,vr,wr,pr,hr,&
                            f1,f2,f3,f4,f5,&
                            ra,ua,va,wa,pa,ha,&
                            c1,c2,c3,c4,c5,c6,c7,c8,c9,&
                            dun,uan,&
                            eig1,eig2,eig3,&
                            ca2,ca3,ca,&
                            prpl,eigra,&
                            unr,unl,t1,t2,t3,t4,t5,t6,t7,t8,&
                            t9,t10,t11,t12,t13,t14,t15,t16,t17,&
                            t18,t19,t20

#if 0
        zero   = 1.e-15
        !>roe-averaged
        !>
        gx = ax*aa
        gy = ay*aa
        gz = az*aa
        !>
        gagm = gamma / gm1
        c1   = sqrt(rr/rl)
        c2   = 1./(1.+c1)
        hl   = gagm*pl/rl+0.5*(ul*ul+vl*vl+wl*wl)
        hr   = gagm*pr/rr+0.5*(ur*ur+vr*vr+wr*wr)
        ra =sqrt(rr*rl)
        ua =(ul+ur*c1)*c2
        va =(vl+vr*c1)*c2
        wa =(wl+wr*c1)*c2
        ha =(hl+hr*c1)*c2
        ca =sqrt(gm1*(ha-0.5*(ua*ua+va*va+wa*wa)))

        dun =ax*(ur-ul)+ay*(vr-vl)+az*(wr-wl)
        uan =ax*ua     +ay*va     +az*wa
        !>entropy fix for eigenvalue (eig1,eig2 and eig3)
        !> roe schemes for all speed methods
        !>
#if 0
        lmach = (uban/aa)/ca
        m1   = 1 - lmach*lmach
        m2   = 1 + lmach*lmach
        m3   = lmach*sqrt(4+m1*m1)/m2
        allspeedroe = min(1.e0,m3)
        ca = allspeedroe * ca
        !> roe schemes for all speed methods
        !>
        !>
        if(xmach .lt. 0.3e0)then
            m1   = 1 - xmach*xmach
            m2   = 1 + xmach*xmach
            m3   = xmach*sqrt(4+m1*m1)/m2
            allspeedroe = min(1.e0,m3)
            ca = allspeedroe * ca
        end if
#endif
        !>

        eig1  =abs(uan)
        eig2  =abs(uan+ca)
        eig3  =abs(uan-ca)
        !> used the limiter eigenvalues a la-harten and gnoffo
        !> implicit tvd schemes for hyperbolic conservation laws in curvilinear coordinates
        !> the shock-caputring schemes
        !> 0.01  < epsa_r < 0.4
        !> 0.1   < epsaa  < 0.3     approximate range for blunt bodies
        !> 0.005 < epsaa  < 0.05    approximate range for slender bodies
        if(epsa_r .gt. 0)then
            epsaa = epsa_r*(ca + abs(ua) + abs(va) + abs(wa))
            epsbb = 0.25/max(epsaa,zero)
            epscc = 2.00*epsaa
            if(real(eig1) .lt. real(epscc)) eig1= eig1*eig1*epsbb+epsaa
            if(real(eig2) .lt. real(epscc)) eig2= eig2*eig2*epsbb+epsaa
            if(real(eig3) .lt. real(epscc)) eig3= eig3*eig3*epsbb+epsaa
        end if

        !> calculation of the |a|(qr-ql) terms

        ca2    = ca*ca
        ca3    = ra*ca*dun
        prpl   = pr-pl
        eigra  = eig1*ra

        c1 =eig1*(rr-rl-prpl/ca2)
        c2 =0.5*eig2*(prpl+ca3)/ca2
        c3 =0.5*eig3*(prpl-ca3)/ca2
        c4 =c1+c2+c3
        c5 =ca*(c2-c3)
        c6 =eigra*(ur-ul-ax*dun)
        c7 =eigra*(vr-vl-ay*dun)
        c8 =eigra*(wr-wl-az*dun)
        c9 =ha*c4+uan*c5+ua*c6+va*c7+wa*c8-c1*ca2*gm
        !> calculation the flux  f(i-1/2) = 0.5*( fr+fl -|a|(qr-ql))
        unr =gx*ur+gy*vr+gz*wr
        unl =gx*ul+gy*vl+gz*wl

        rr  =rr*unr
        rl  =rl*unl
        pr  =pr+pl
!        write(120,'("n6=",5e24.16)') UNR,UNL,RR,RL,PR
        f1 =0.5*(rr   +rl               - aa*c4)
        f2 =0.5*(rr*ur+rl*ul+gx*pr      - aa*(ua*c4+ax*c5+c6))
        f3 =0.5*(rr*vr+rl*vl+gy*pr      - aa*(va*c4+ay*c5+c7))
        f4 =0.5*(rr*wr+rl*wl+gz*pr      - aa*(wa*c4+az*c5+c8))
        f5 =0.5*(rr*hr+rl*hl            - aa*c9)
#else
        t1 = rr-rl
        t2 = ur-ul
        t3 = vr-vl
        t4 = wr-wl
        !>pressure and enthalpy
        t16  = 1.e0/rr
        t5   = pr
        pr   = (gamma/gm1)*pr*t16 + 0.5e0*(ur*ur+vr*vr+wr*wr)
        !>
        t15     = 1.e0/rl
        t19     = pl
        pl      = (gamma/gm1)*pl*t15+0.5e0*(ul*ul+vl*vl+wl*wl)
        !>unsplit contributions  f(r)+f(l)
        t18 = ax*ur+ay*vr+az*wr+0.0!at
        t17 = ax*ul+ay*vl+az*wl+0.0!at
        t6  = t18*rr
        t7  = t17*rl
        f1  = t6+t7
        f2  = t6*ur+t7*ul
        f3  = t6*vr+t7*vl
        f4  = t6*wr+t7*wl
        f5  = t6*pr+t7*pl
        t8  = t5+t19
        f2  = f2+ax*t8
        f3  = f3+ay*t8
        f4  = f4+az*t8
        f5  = f5!-at*t8
        !>roe averaged variables
        t6 = rr*t15
        t7 = sqrt(t6)
        t6 = 1.e0/(1.e0+t7)
        t8 = t7*t6
        !>average density
        rr = rl*t7
        !>u,v,w,h average
        t9  = ul*t6+ur*t8
        t10 = vl*t6+vr*t8
        t11 = wl*t6+wr*t8
        t12 = pl*t6+pr*t8
        !>extract sound speed
        t6  = (t9*t9+t10*t10+t11*t11)*0.5e0
        t7  = gm1*(t12-t6)
        t8  = sqrt(t7)
        !>
        t13 = t9*ax+t10*ay+t11*az
!        if(idir ==2) write(*,'(7e24.16,3e30.20)')t13,t9,t10,t11,ax,ay,az,t9*ax,t10*ay,t11*az
        rl  = rr*(t18-t17)
        ul  = (t5-t19)/t7
        !>
        t18 = t13!+at(i)
        t18 = abs(t18)
        t19 = t13+t8!+at
        t19 = abs(t19)
        t17 = t13-t8!+at
        t17 = abs(t17)
        !>limit eigenvalues
        if(epsa_r .gt. 0)then
            epsaa = epsa_r*(t8 + abs(t9) + abs(t10) + abs(t11))
            epsbb = 0.25/max(epsaa,zero)
            epscc = 2.00*epsaa
            if(real(t18) .lt. real(epscc)) t18= t18*t18*epsbb+epsaa
            if(real(t17) .lt. real(epscc)) t17= t17*t17*epsbb+epsaa
            if(real(t19) .lt. real(epscc)) t19= t19*t19*epsbb+epsaa
        end if
        !>
        t14 = t18*(t1-ul)
        t15 = 0.5e0*(ul+rl/t8)
        t16 = (ul-t15)*t17
        t15 = t15*t19
        !>
        ur = t18*(t2*rr-ax*rl)
        vr = t18*(t3*rr-ay*rl)
        wr = t18*(t4*rr-az*rl)
        pr = t9*ur+t10*vr+t11*wr
        !>
        rl = t14+t15+t16
        ul = t8*(t15-t16)
        !>
!        if(idir ==2) write(*,'(9e24.16)') f1,f5,t13,ax,ay,az,t17,t18,t19
        f1 = f1-rl
        f2 = f2-rl*t9 -ax*ul-ur
        f3 = f3-rl*t10-ay*ul-vr
        f4 = f4-rl*t11-az*ul-wr
        f5 = f5-rl*t12-t13*ul-pr+t7*t14/gm1
        !>
        t7 = 0.5e0*aa
        f1 = t7*f1
        f2 = t7*f2
        f3 = t7*f3
        f4 = t7*f4
        f5 = t7*f5
#endif
        return
    end subroutine roe
