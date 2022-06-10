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

        return
    end subroutine roe
