!>***********************************************************************
!>     purpose:  return a value of a or b
!>     depending the logical expression c
!>     when c is ture return a or return b
!>     by liuxz 2017/3/8
!>***********************************************************************
      function logical_depending(a,b,c)
        !>
        logical c
        real*8::a,b,logical_depending
		if(c)then
            logical_depending = a
        else
            logical_depending = b
		endif
		!>
		return

      end
!>***********************************************************************
      !> purpose:defined  the limiter functions
      !> used min-mod limiter
      !> input a and b and return the lim_minmod
!>***********************************************************************
      function lim_minmod(x,y,b)
        !>
        real*8::x,y,lim_minmod,b
!        lim_minmod=max(0.0,min(abs(a),b*sign(1.0,a)))*sign(1.0,a)
        lim_minmod=max(0.0,min(x*sign(1.0,y),b*y*sign(1.0,x)))*sign(1.0,x)
        !>test the minmod method
        !>
!        write(*,'("minmod=",e20.14)') max(0.0,min(0.15341*sign(1.0,0.17543),b*0.17543*sign(1.0,0.15341)))*sign(1.0,0.15341)
!        write(*,'("minmod=",e20.14)') max(0.0,min(0.17543*sign(1.0,0.15341),b*0.15341*sign(1.0,0.17543)))*sign(1.0,0.17543)
        !>
        !>
!        lim_minmod=max(0.0,min(x*sign(y),b*y*sign(x)))*sign(x)
		!>
		return

      end
