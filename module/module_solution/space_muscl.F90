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
    subroutine muscl(nbl)
		!>
		!>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
        use variable_module
        use bc_module
        use metric_module
        use nodes_paras
        !>
        !>
		implicit none
		!>
		!>
        integer :: nbl,&
                   i,j,k,&
                   idim,jdim,kdim,&
                   max_ijk,&
                   error,nwall,&
                   ii,jj,kk,&
                   ipw,jpw,kpw,&
                   iii,jjj,kkk,&
                   id1,jd1,kd1,&
                   minr1,maxr1,minr2,maxr2
        real(kind = dprec),dimension(:),allocatable :: du,dv,dw,dp,dr
        real(kind = dprec),dimension(:),allocatable :: fl1,fl2,fl3,fl4,fl5
        real(kind = dprec),dimension(:),allocatable :: fr1,fr2,fr3,fr4,fr5
        real(kind = dprec),dimension(:),allocatable :: sl1,sl2,sl3,sl4,sl5
        real(kind = dprec),dimension(:),allocatable :: sr1,sr2,sr3,sr4,sr5
        real(kind = dprec),dimension(:),allocatable :: f1,f2,f3,f4,f5
        real(kind = dprec),dimension(:),allocatable :: frr,fpp


        real(kind = dprec) :: rl,ul,vl,wl,pl,&
                              rr,ur,vr,wr,pr
        !>
        type(blocks_type),pointer  :: iblock
        type(bc_types),pointer     :: ibc
        type(overlap_type),pointer :: mesh


        !>
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        !>
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        !>
        !>
        max_ijk = 0
        max_ijk = max(idim,jdim)
        max_ijk = max(kdim,max_ijk)+3
        !>
        allocate(du(0:max_ijk))
        allocate(dv(0:max_ijk))
        allocate(dw(0:max_ijk))
        allocate(dr(0:max_ijk))
        allocate(dp(0:max_ijk))
        !>
        allocate(fl1(max_ijk))
        allocate(fl2(max_ijk))
        allocate(fl3(max_ijk))
        allocate(fl4(max_ijk))
        allocate(fl5(max_ijk))
        allocate(fr1(max_ijk))
        allocate(fr2(max_ijk))
        allocate(fr3(max_ijk))
        allocate(fr4(max_ijk))
        allocate(fr5(max_ijk))
        !>
        allocate(frr(max_ijk))
        allocate(fpp(max_ijk))
        !>
        allocate(sl1(max_ijk))
        allocate(sl2(max_ijk))
        allocate(sl3(max_ijk))
        allocate(sl4(max_ijk))
        allocate(sl5(max_ijk))
        allocate(sr1(max_ijk))
        allocate(sr2(max_ijk))
        allocate(sr3(max_ijk))
        allocate(sr4(max_ijk))
        allocate(sr5(max_ijk))
        allocate(f1(max_ijk))
        allocate(f2(max_ijk))
        allocate(f3(max_ijk))
        allocate(f4(max_ijk))
        allocate(f5(max_ijk))
        !>
        !>
        !>solution of fluxs in three-directions and composition of the fluxs
        !--------------------------------------------------------------------
        !> flux in i-directions
        !--------------------------------------------------------------------
        do k=2,kdim
            do j=2,jdim
                !>
                !>
                !> operations for limiter
                do i=3,idim
                    dr(i)  =( iblock%cell(i,j,k)%r - iblock%cell(i-1,j,k)%r )!*iblock%cell(i,j,k)%wall_blank
                    du(i)  =( iblock%cell(i,j,k)%u - iblock%cell(i-1,j,k)%u )!*iblock%cell(i,j,k)%wall_blank
                    dv(i)  =( iblock%cell(i,j,k)%v - iblock%cell(i-1,j,k)%v )!*iblock%cell(i,j,k)%wall_blank
                    dw(i)  =( iblock%cell(i,j,k)%w - iblock%cell(i-1,j,k)%w )!*iblock%cell(i,j,k)%wall_blank
                    dp(i)  =( iblock%cell(i,j,k)%p - iblock%cell(i-1,j,k)%p )!*iblock%cell(i,j,k)%wall_blank
                    !> store the pressure and density for gradients by average
                    frr(i)  =iblock%cell(i,j,k)%r
                    fpp(i)  =iblock%cell(i,j,k)%p
                end do
                !> left boundary values
                if(iblock%cell(1,j,k)%wall_blank .eq. 2)then
                    dr(2)  = iblock%cell(0,j,k)%r
                    du(2)  = iblock%cell(0,j,k)%u
                    dv(2)  = iblock%cell(0,j,k)%v
                    dw(2)  = iblock%cell(0,j,k)%w
                    dp(2)  = iblock%cell(0,j,k)%p
                else
                    dr(2)  = iblock%cell(2,j,k)%r - iblock%cell(1,j,k)%r
                    du(2)  = iblock%cell(2,j,k)%u - iblock%cell(1,j,k)%u
                    dv(2)  = iblock%cell(2,j,k)%v - iblock%cell(1,j,k)%v
                    dw(2)  = iblock%cell(2,j,k)%w - iblock%cell(1,j,k)%w
                    dp(2)  = iblock%cell(2,j,k)%p - iblock%cell(1,j,k)%p
                end if
                !> right boundary values
                if(iblock%cell(idim+1,j,k)%wall_blank .eq. 2)then
                    dr(idim+1)  = iblock%cell(idim+2,j,k)%r
                    du(idim+1)  = iblock%cell(idim+2,j,k)%u
                    dv(idim+1)  = iblock%cell(idim+2,j,k)%v
                    dw(idim+1)  = iblock%cell(idim+2,j,k)%w
                    dp(idim+1)  = iblock%cell(idim+2,j,k)%p
                else
                    dr(idim+1)  = iblock%cell(idim+1,j,k)%r - iblock%cell(idim,j,k)%r
                    du(idim+1)  = iblock%cell(idim+1,j,k)%u - iblock%cell(idim,j,k)%u
                    dv(idim+1)  = iblock%cell(idim+1,j,k)%v - iblock%cell(idim,j,k)%v
                    dw(idim+1)  = iblock%cell(idim+1,j,k)%w - iblock%cell(idim,j,k)%w
                    dp(idim+1)  = iblock%cell(idim+1,j,k)%p - iblock%cell(idim,j,k)%p
                end if
                frr(1)      =   iblock%cell(1,j,k)%r
                fpp(1)      =   iblock%cell(1,j,k)%p
                frr(2)      =   iblock%cell(2,j,k)%r
                fpp(2)      =   iblock%cell(2,j,k)%p
                frr(idim)   =   iblock%cell(idim,j,k)%r
                fpp(idim)   =   iblock%cell(idim,j,k)%p
                frr(idim+1) =   iblock%cell(idim+1,j,k)%r
                fpp(idim+1) =   iblock%cell(idim+1,j,k)%p
                !>
                !>
                !>left boundary
                !>
                dr(0) = iblock%cell(1,j,k)%r - iblock%cell(0,j,k)%r
                du(0) = iblock%cell(1,j,k)%u - iblock%cell(0,j,k)%u
                dv(0) = iblock%cell(1,j,k)%v - iblock%cell(0,j,k)%v
                dw(0) = iblock%cell(1,j,k)%w - iblock%cell(0,j,k)%w
                dp(0) = iblock%cell(1,j,k)%p - iblock%cell(0,j,k)%p
                !>
                !>
                dr(1) = iblock%cell(2,j,k)%r - iblock%cell(1,j,k)%r
                du(1) = iblock%cell(2,j,k)%u - iblock%cell(1,j,k)%u
                dv(1) = iblock%cell(2,j,k)%v - iblock%cell(1,j,k)%v
                dw(1) = iblock%cell(2,j,k)%w - iblock%cell(1,j,k)%w
                dp(1) = iblock%cell(2,j,k)%p - iblock%cell(1,j,k)%p
                !>right boundary
                !>
                dr(idim+2) = iblock%cell(idim+1,j,k)%r - iblock%cell(idim,j,k)%r
                du(idim+2) = iblock%cell(idim+1,j,k)%u - iblock%cell(idim,j,k)%u
                dv(idim+2) = iblock%cell(idim+1,j,k)%v - iblock%cell(idim,j,k)%v
                dw(idim+2) = iblock%cell(idim+1,j,k)%w - iblock%cell(idim,j,k)%w
                dp(idim+2) = iblock%cell(idim+1,j,k)%p - iblock%cell(idim,j,k)%p
                !>
                !>
                dr(idim+3) = iblock%cell(idim+2,j,k)%r - iblock%cell(idim+1,j,k)%r
                du(idim+3) = iblock%cell(idim+2,j,k)%u - iblock%cell(idim+1,j,k)%u
                dv(idim+3) = iblock%cell(idim+2,j,k)%v - iblock%cell(idim+1,j,k)%v
                dw(idim+3) = iblock%cell(idim+2,j,k)%w - iblock%cell(idim+1,j,k)%w
                dp(idim+3) = iblock%cell(idim+2,j,k)%p - iblock%cell(idim+1,j,k)%p
                !>
                !>
                call limit(idim+1,max_ijk,frr,fpp,du,dv,dw,dp,dr,fl1,fl2,fl3,fl4,fl5,fr1,fr2,fr3,fr4,fr5)

                !>  upwind step

                do i  = 2,idim+1

                    !> interpolation of the primitive variables q=(r,u,v,w,p)
                    !> i0 face
                    if(iblock%cell(i-1,j,k)%wall_blank .eq. 2.0 .and. i-1 .lt. 2)then
                        !>
                        rl = iblock%cell(i-1,j,k)%r
                        ul = iblock%cell(i-1,j,k)%u
                        vl = iblock%cell(i-1,j,k)%v
                        wl = iblock%cell(i-1,j,k)%w
                        pl = iblock%cell(i-1,j,k)%p
                        !>right values
                        rr =iblock%cell( i-1,j,k)%r! - fl1(i)
                        ur =iblock%cell( i-1,j,k)%u! - fl2(i)
                        vr =iblock%cell( i-1,j,k)%v! - fl3(i)
                        wr =iblock%cell( i-1,j,k)%w! - fl4(i)
                        pr =iblock%cell( i-1,j,k)%p! - fl5(i)
                    elseif(iblock%cell(i,j,k)%wall_blank .eq. 2.0  .and.  i .gt. 2)then
                        !> idim face
                        rl =iblock%cell(i,j,k)%r !+ fr1(i-1)
                        ul =iblock%cell(i,j,k)%u !+ fr2(i-1)
                        vl =iblock%cell(i,j,k)%v !+ fr3(i-1)
                        wl =iblock%cell(i,j,k)%w !+ fr4(i-1)
                        pl =iblock%cell(i,j,k)%p !+ fr5(i-1)
                        !>right values
                        !>
                        rr =iblock%cell(i,j,k)%r
                        ur =iblock%cell(i,j,k)%u
                        vr =iblock%cell(i,j,k)%v
                        wr =iblock%cell(i,j,k)%w
                        pr =iblock%cell(i,j,k)%p
                    else
                        !> inner grid
                        rl =iblock%cell(i-1,j,k)%r + fr1(i-1)
                        ul =iblock%cell(i-1,j,k)%u + fr2(i-1)
                        vl =iblock%cell(i-1,j,k)%v + fr3(i-1)
                        wl =iblock%cell(i-1,j,k)%w + fr4(i-1)
                        pl =iblock%cell(i-1,j,k)%p + fr5(i-1)

                        rr =iblock%cell( i,j,k)%r - fl1(i)
                        ur =iblock%cell( i,j,k)%u - fl2(i)
                        vr =iblock%cell( i,j,k)%v - fl3(i)
                        wr =iblock%cell( i,j,k)%w - fl4(i)
                        pr =iblock%cell( i,j,k)%p - fl5(i)
                    end if
!                    write(1411,'(3i4,10e20.12)')i,j,k,fl1(i-1),fl2(i-1),fl3(i-1),fl4(i-1),fl5(i-1),&
!                                                fr1(i),fr2(i),fr3(i),fr4(i),fr5(i)
                    !>
                    !>
                    !>modify the pressure and density
                    !>
                    if(pr<=0.)then
                        write(*,'("nbl=",i5," i-dir=",3i5," pr=",3e20.12)') nbl,i,j,k,pr,iblock%cell( i,j,k)%p,fr5(i)
                        pr = p0
                    end if
                    if(pl<=0.)then
                        write(*,'("nbl=",i5," i-dir=",3i5," pl=",3e20.12)') nbl,i-1,j,k,pl,iblock%cell( i-1,j,k)%p,fl5(i-1)
                        pl = p0
                    end if

                    if(rr<=0.)then
                        write(*,'("nbl=",i5," i-dir=",3i5," rr=",3e20.12)') nbl,i,j,k,rr,iblock%cell( i,j,k)%r,fr1(i)
                        rr = r0
                    end if
                    if(rl<=0.)then
                        write(*,'("nbl=",i5," i-dir=",3i5," rl=",3e20.12)') nbl,i-1,j,k,rl,iblock%cell(i-1,j,k)%r ,fl1(i-1)
                        rl = r0
                    end if
                    !>
                    !> upwind solution
                    !>
                    !>
!                    write(1411,'(3i4,10e20.12)')i,j,k,rl,ul,vl,wl,pl,rr,ur,vr,wr,pr
!                    write(14111,'(3i4,4e20.12)')i,j,k,iblock%metric(i,j,k)%fi,iblock%metric(i,j,k)%fj,iblock%metric(i,j,k)%fk,iblock%metric(i,j,k)%ff
!                    write(131,'("axis=",3i4,7e24.16)')i,j,k,dr(i),du(i),dv(i),dw(i),dp(i),frr(i),fpp(i)
!                    write(1411,'(3i4,5e24.16)') i,j,k,RL,UL,VL,WL,PL
!                    write(1412,'(3i4,5e24.16)') i,j,k,RR,UR,VR,WR,PR
                    call roe(iblock%metric(i,j,k)%fi,iblock%metric(i,j,k)%fj,iblock%metric(i,j,k)%fk,iblock%metric(i,j,k)%ff,&
                             rl,ul,vl,wl,pl,&
                             rr,ur,vr,wr,pr,&
                             f1(i),f2(i),f3(i),f4(i),f5(i))
!                    write(141,'(3i4,5e24.16)')i,j,k,f1(i),f2(i),f3(i),f4(i),f5(i)
                end do
                !>
                !>
                !> composition of the fluxes
                do i=2,idim
                    !>
                    !>
                    iblock%variable(i,j,k)%res_1 =    f1(i+1) -   f1(i)
                    iblock%variable(i,j,k)%res_2 =    f2(i+1) -   f2(i)
                    iblock%variable(i,j,k)%res_3 =    f3(i+1) -   f3(i)
                    iblock%variable(i,j,k)%res_4 =    f4(i+1) -   f4(i)
                    iblock%variable(i,j,k)%res_5 =    f5(i+1) -   f5(i)
!                    if(nbl == 1) write(14221,'(3i4,10e24.16)')i,j,k,f1(i+1),f1(i),f2(i+1),f2(i),f3(i+1),f3(i),&
!                    f4(i+1),f4(i),f5(i+1),f5(i)
!                    if(nbl == 3) write(142213,'(3i4,5e24.16)')i,j,k,f1(i+1)-f1(i),f2(i+1)-f2(i),f3(i+1)-f3(i),&
!                    f4(i+1)-f4(i),f5(i+1)-f5(i)
                end do
                !>
                !>
            end do !>jdim
        end do  ! >kdim
        !--------------------------------------------------------------------
        !>  flux in j-directions
        !--------------------------------------------------------------------
        do k=2,kdim
            do i=2,idim
                !>
                !>
                !> operations for limiter
                do j=3,jdim
                    dr(j)  =( iblock%cell(i,j,k)%r - iblock%cell(i,j-1,k)%r )
                    du(j)  =( iblock%cell(i,j,k)%u - iblock%cell(i,j-1,k)%u )
                    dv(j)  =( iblock%cell(i,j,k)%v - iblock%cell(i,j-1,k)%v )
                    dw(j)  =( iblock%cell(i,j,k)%w - iblock%cell(i,j-1,k)%w )
                    dp(j)  =( iblock%cell(i,j,k)%p - iblock%cell(i,j-1,k)%p )
                    !> store the pressure and density for gradients by average
                    frr(j)  =iblock%cell(i,j,k)%r
                    fpp(j)  =iblock%cell(i,j,k)%p
                end do
                !> left boundary values
                if(iblock%cell(i,1,k)%wall_blank .eq. 2)then
                    dr(2)  = iblock%cell(i,0,k)%r
                    du(2)  = iblock%cell(i,0,k)%u
                    dv(2)  = iblock%cell(i,0,k)%v
                    dw(2)  = iblock%cell(i,0,k)%w
                    dp(2)  = iblock%cell(i,0,k)%p
                else
                    dr(2)  = iblock%cell(i,2,k)%r - iblock%cell(i,1,k)%r
                    du(2)  = iblock%cell(i,2,k)%u - iblock%cell(i,1,k)%u
                    dv(2)  = iblock%cell(i,2,k)%v - iblock%cell(i,1,k)%v
                    dw(2)  = iblock%cell(i,2,k)%w - iblock%cell(i,1,k)%w
                    dp(2)  = iblock%cell(i,2,k)%p - iblock%cell(i,1,k)%p
                end if
                !> right boundary values
                if(iblock%cell(i,jdim+1,k)%wall_blank .eq. 2)then
                    dr(jdim+1)  = iblock%cell(i,jdim+2,k)%r
                    du(jdim+1)  = iblock%cell(i,jdim+2,k)%u
                    dv(jdim+1)  = iblock%cell(i,jdim+2,k)%v
                    dw(jdim+1)  = iblock%cell(i,jdim+2,k)%w
                    dp(jdim+1)  = iblock%cell(i,jdim+2,k)%p
                else
                    dr(jdim+1)  = iblock%cell(i,jdim+1,k)%r - iblock%cell(i,jdim,k)%r
                    du(jdim+1)  = iblock%cell(i,jdim+1,k)%u - iblock%cell(i,jdim,k)%u
                    dv(jdim+1)  = iblock%cell(i,jdim+1,k)%v - iblock%cell(i,jdim,k)%v
                    dw(jdim+1)  = iblock%cell(i,jdim+1,k)%w - iblock%cell(i,jdim,k)%w
                    dp(jdim+1)  = iblock%cell(i,jdim+1,k)%p - iblock%cell(i,jdim,k)%p
                end if
                frr(1)      =   iblock%cell(i,1,k)%r
                fpp(1)      =   iblock%cell(i,1,k)%p
                frr(2)      =   iblock%cell(i,2,k)%r
                fpp(2)      =   iblock%cell(i,2,k)%p
                frr(jdim)   =   iblock%cell(i,jdim,k)%r
                fpp(jdim)   =   iblock%cell(i,jdim,k)%p
                frr(jdim+1) =   iblock%cell(i,jdim+1,k)%r
                fpp(jdim+1) =   iblock%cell(i,jdim+1,k)%p
                !>
                !>
                !>left boundary
                !>
                dr(0) = iblock%cell(i,1,k)%r - iblock%cell(i,0,k)%r
                du(0) = iblock%cell(i,1,k)%u - iblock%cell(i,0,k)%u
                dv(0) = iblock%cell(i,1,k)%v - iblock%cell(i,0,k)%v
                dw(0) = iblock%cell(i,1,k)%w - iblock%cell(i,0,k)%w
                dp(0) = iblock%cell(i,1,k)%p - iblock%cell(i,0,k)%p
                !>
                !>
                dr(1) = iblock%cell(i,2,k)%r - iblock%cell(i,1,k)%r
                du(1) = iblock%cell(i,2,k)%u - iblock%cell(i,1,k)%u
                dv(1) = iblock%cell(i,2,k)%v - iblock%cell(i,1,k)%v
                dw(1) = iblock%cell(i,2,k)%w - iblock%cell(i,1,k)%w
                dp(1) = iblock%cell(i,2,k)%p - iblock%cell(i,1,k)%p
                !>right boundary
                !>
                dr(jdim+2) = iblock%cell(i,jdim+1,k)%r - iblock%cell(i,jdim,k)%r
                du(jdim+2) = iblock%cell(i,jdim+1,k)%u - iblock%cell(i,jdim,k)%u
                dv(jdim+2) = iblock%cell(i,jdim+1,k)%v - iblock%cell(i,jdim,k)%v
                dw(jdim+2) = iblock%cell(i,jdim+1,k)%w - iblock%cell(i,jdim,k)%w
                dp(jdim+2) = iblock%cell(i,jdim+1,k)%p - iblock%cell(i,jdim,k)%p
                !>
                !>
                dr(jdim+3) = iblock%cell(i,jdim+2,k)%r - iblock%cell(i,jdim+1,k)%r
                du(jdim+3) = iblock%cell(i,jdim+2,k)%u - iblock%cell(i,jdim+1,k)%u
                dv(jdim+3) = iblock%cell(i,jdim+2,k)%v - iblock%cell(i,jdim+1,k)%v
                dw(jdim+3) = iblock%cell(i,jdim+2,k)%w - iblock%cell(i,jdim+1,k)%w
                dp(jdim+3) = iblock%cell(i,jdim+2,k)%p - iblock%cell(i,jdim+1,k)%p
                !>
                !>
                call limit(jdim+1,max_ijk,frr,fpp,du,dv,dw,dp,dr,fl1,fl2,fl3,fl4,fl5,fr1,fr2,fr3,fr4,fr5)
                !>
                !>
                !>
                do j=2,jdim+1
                    !>
                    !> interpolation of the primitive variables q=(r,u,v,w,p)
                    !> j0 face
                    if(iblock%cell(i,j-1,k)%wall_blank .eq. 2.0 .and. j-1 .lt. 2)then
    !                    write(*,*) i,jn,k,wall_blank(jn,j,k),wall_blank(i,jn,k)
                        !>
                        rl =iblock%cell(i,j-1,k)%r
                        ul =iblock%cell(i,j-1,k)%u
                        vl =iblock%cell(i,j-1,k)%v
                        wl =iblock%cell(i,j-1,k)%w
                        pl =iblock%cell(i,j-1,k)%p
                        !>right values
                        !>
                        rr =iblock%cell( i,j-1,k)%r !- fl1(j)
                        ur =iblock%cell( i,j-1,k)%u !- fl2(j)
                        vr =iblock%cell( i,j-1,k)%v !- fl3(j)
                        wr =iblock%cell( i,j-1,k)%w !- fl4(j)
                        pr =iblock%cell( i,j-1,k)%p !- fl5(j)
                    elseif(iblock%cell(i,j,k)%wall_blank .eq. 2.0  .and.  j .gt. 2)then
                        !> jdim face
                        rl =iblock%cell(i,j,k)%r !+ fr1(j-1)
                        ul =iblock%cell(i,j,k)%u !+ fr2(j-1)
                        vl =iblock%cell(i,j,k)%v !+ fr3(j-1)
                        wl =iblock%cell(i,j,k)%w !+ fr4(j-1)
                        pl =iblock%cell(i,j,k)%p !+ fr5(j-1)

                        !>right values
                        !>
                        rr =iblock%cell(i,j,k)%r
                        ur =iblock%cell(i,j,k)%u
                        vr =iblock%cell(i,j,k)%v
                        wr =iblock%cell(i,j,k)%w
                        pr =iblock%cell(i,j,k)%p

                    else
                        !> inner grid
                        rl =iblock%cell(i,j-1,k)%r + fr1(j-1)
                        ul =iblock%cell(i,j-1,k)%u + fr2(j-1)
                        vl =iblock%cell(i,j-1,k)%v + fr3(j-1)
                        wl =iblock%cell(i,j-1,k)%w + fr4(j-1)
                        pl =iblock%cell(i,j-1,k)%p + fr5(j-1)
                        !> inner grid
                        rr =iblock%cell( i,j,k)%r - fl1(j)
                        ur =iblock%cell( i,j,k)%u - fl2(j)
                        vr =iblock%cell( i,j,k)%v - fl3(j)
                        wr =iblock%cell( i,j,k)%w - fl4(j)
                        pr =iblock%cell( i,j,k)%p - fl5(j)
                    end if
!                    write(1422,'(3i4,10e20.12)')i,j,k,fl1(j-1),fl2(j-1),fl3(j-1),fl4(j-1),fl5(j-1),&
!                                                fr1(j),fr2(j),fr3(j),fr4(j),fr5(j)
                    !>
                    !>
                    !>modify the pressure and density
                    !>
                    if(pr<=0.)then
                        write(*,'("nbl=",i5," j-dir=",3i5," pr=",3e20.12)') nbl,i,j,k,pr,iblock%cell( i,j,k)%p,fr5(j)
                        pr = p0
                    end if
                    if(pl<=0.)then
                        write(*,'("nbl=",i5," j-dir=",3i5," pl=",3e20.12)') nbl,i,j-1,k,pl,iblock%cell( i,j-1,k)%p,fl5(j-1)
                        pl = p0
                    end if

                    if(rr<=0.)then
                        write(*,'("nbl=",i5," j-dir=",3i5," rr=",3e20.12)') nbl,i,j,k,rr,iblock%cell( i,j,k)%r,fr1(j)
                        rr = r0
                    end if
                    if(rl<=0.)then
                        write(*,'("nbl=",i5," j-dir=",3i5," rl=",3e20.12)') nbl,i,j-1,k,rl,iblock%cell( i,j-1,k)%r,fl1(j-1)
                        rl = r0
                    end if


                    !>upwind solution
!                    write(1422,'(3i4,5e24.16)')i,j,k,f1(j),f2(j),f3(j),f4(j),f5(j)
!                    write(14222,'(3i4,14e24.16)')i,j,k,iblock%metric(i,j,k)%gi,iblock%metric(i,j,k)%gj,iblock%metric(i,j,k)%gk,iblock%metric(i,j,k)%gg,&
!                    rl,ul,vl,wl,pl,rr,ur,vr,wr,pr
                    call roe(iblock%metric(i,j,k)%gi,iblock%metric(i,j,k)%gj,iblock%metric(i,j,k)%gk,iblock%metric(i,j,k)%gg,&
                                           rl,ul,vl,wl,pl,&
                                           rr,ur,vr,wr,pr,&
                                           f1(j),f2(j),f3(j),f4(j),f5(j))
!                    write(142,'(3i4,5e24.16)')i,j,k,f1(j),f2(j),f3(j),f4(j),f5(j)
                end do
                !>
                !>
                !> wall boundary condition
                !>

                do j=2,jdim
!                    write(14223,'(3i4,5e24.16)') i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,&
!                    iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1  + (f1(j+1) -   f1(j))
                    iblock%variable(i,j,k)%res_2 = iblock%variable(i,j,k)%res_2  + (f2(j+1) -   f2(j))
                    iblock%variable(i,j,k)%res_3 = iblock%variable(i,j,k)%res_3  + (f3(j+1) -   f3(j))
                    iblock%variable(i,j,k)%res_4 = iblock%variable(i,j,k)%res_4  + (f4(j+1) -   f4(j))
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_5  + (f5(j+1) -   f5(j))
!                    if(nbl == 161 .and. icyc ==2) write(14222,'(3i4,10e24.16)')i,j,k,f1(j+1),f1(j),f2(j+1),f2(j),f3(j+1),f3(j),&
!                    f4(j+1),f4(j),f5(j+1),f5(j)
!                    if(nbl == 161 .and.  icyc ==2) write(142223,'(3i4,5e24.16)')i,j,k,f1(j+1)-f1(j),f2(j+1)-f2(j),f3(j+1)-f3(j),&
!                    f4(j+1)-f4(j),f5(j+1)-f5(j)
!                    if(nbl == 108 .and.  icyc ==2) write(1422231,'(3i4,5e24.16)')i,j,k,f1(j+1)-f1(j),f2(j+1)-f2(j),f3(j+1)-f3(j),&
!                    f4(j+1)-f4(j),f5(j+1)-f5(j)
                end do

            end do
        end do
        !--------------------------------------------------------------------
        !> flux in k-directions
        !--------------------------------------------------------------------
        do j=2,jdim
            do i=2,idim
                !>
                !> operations for limiter
                do k=3,kdim
                    dr(k)  =( iblock%cell(i,j,k)%r - iblock%cell(i,j,k-1)%r )
                    du(k)  =( iblock%cell(i,j,k)%u - iblock%cell(i,j,k-1)%u )
                    dv(k)  =( iblock%cell(i,j,k)%v - iblock%cell(i,j,k-1)%v )
                    dw(k)  =( iblock%cell(i,j,k)%w - iblock%cell(i,j,k-1)%w )
                    dp(k)  =( iblock%cell(i,j,k)%p - iblock%cell(i,j,k-1)%p )
                    !> store the pressure and density for gradients by average
                    frr(k)  =iblock%cell(i,j,k)%r
                    fpp(k)  =iblock%cell(i,j,k)%p
                end do
                !> left boundary values
                if(iblock%cell(i,j,1)%wall_blank .eq. 2)then
                    dr(2)  = iblock%cell(i,j,0)%r
                    du(2)  = iblock%cell(i,j,0)%u
                    dv(2)  = iblock%cell(i,j,0)%v
                    dw(2)  = iblock%cell(i,j,0)%w
                    dp(2)  = iblock%cell(i,j,0)%p
                else
                    dr(2)  = iblock%cell(i,j,2)%r - iblock%cell(i,j,1)%r
                    du(2)  = iblock%cell(i,j,2)%u - iblock%cell(i,j,1)%u
                    dv(2)  = iblock%cell(i,j,2)%v - iblock%cell(i,j,1)%v
                    dw(2)  = iblock%cell(i,j,2)%w - iblock%cell(i,j,1)%w
                    dp(2)  = iblock%cell(i,j,2)%p - iblock%cell(i,j,1)%p
                end if
                !> right boundary values
                if(iblock%cell(i,j,kdim+1)%wall_blank .eq. 2)then
                    dr(kdim+1)  = iblock%cell(i,j,kdim+2)%r
                    du(kdim+1)  = iblock%cell(i,j,kdim+2)%u
                    dv(kdim+1)  = iblock%cell(i,j,kdim+2)%v
                    dw(kdim+1)  = iblock%cell(i,j,kdim+2)%w
                    dp(kdim+1)  = iblock%cell(i,j,kdim+2)%p
                else
                    dr(kdim+1)  = iblock%cell(i,j,kdim+1)%r - iblock%cell(i,j,kdim)%r
                    du(kdim+1)  = iblock%cell(i,j,kdim+1)%u - iblock%cell(i,j,kdim)%u
                    dv(kdim+1)  = iblock%cell(i,j,kdim+1)%v - iblock%cell(i,j,kdim)%v
                    dw(kdim+1)  = iblock%cell(i,j,kdim+1)%w - iblock%cell(i,j,kdim)%w
                    dp(kdim+1)  = iblock%cell(i,j,kdim+1)%p - iblock%cell(i,j,kdim)%p
                end if
                frr(1)      =   iblock%cell(i,j,1)%r
                fpp(1)      =   iblock%cell(i,j,1)%p
                frr(2)      =   iblock%cell(i,j,2)%r
                fpp(2)      =   iblock%cell(i,j,2)%p
                frr(kdim)   =   iblock%cell(i,j,kdim)%r
                fpp(kdim)   =   iblock%cell(i,j,kdim)%p
                frr(kdim+1) =   iblock%cell(i,j,kdim+1)%r
                fpp(kdim+1) =   iblock%cell(i,j,kdim+1)%p
                !>
                !>
                !>left boundary
                !>
                dr(0) = iblock%cell(i,j,1)%r - iblock%cell(i,j,0)%r
                du(0) = iblock%cell(i,j,1)%u - iblock%cell(i,j,0)%u
                dv(0) = iblock%cell(i,j,1)%v - iblock%cell(i,j,0)%v
                dw(0) = iblock%cell(i,j,1)%w - iblock%cell(i,j,0)%w
                dp(0) = iblock%cell(i,j,1)%p - iblock%cell(i,j,0)%p
                !>
                !>
                dr(1) = iblock%cell(i,j,2)%r - iblock%cell(i,j,1)%r
                du(1) = iblock%cell(i,j,2)%u - iblock%cell(i,j,1)%u
                dv(1) = iblock%cell(i,j,2)%v - iblock%cell(i,j,1)%v
                dw(1) = iblock%cell(i,j,2)%w - iblock%cell(i,j,1)%w
                dp(1) = iblock%cell(i,j,2)%p - iblock%cell(i,j,1)%p
                !>right boundary
                !>
                dr(kdim+2) = iblock%cell(i,j,kdim+1)%r - iblock%cell(i,j,kdim)%r
                du(kdim+2) = iblock%cell(i,j,kdim+1)%u - iblock%cell(i,j,kdim)%u
                dv(kdim+2) = iblock%cell(i,j,kdim+1)%v - iblock%cell(i,j,kdim)%v
                dw(kdim+2) = iblock%cell(i,j,kdim+1)%w - iblock%cell(i,j,kdim)%w
                dp(kdim+2) = iblock%cell(i,j,kdim+1)%p - iblock%cell(i,j,kdim)%p
                !>
                !>
                dr(kdim+3) = iblock%cell(i,j,kdim+2)%r - iblock%cell(i,j,kdim+1)%r
                du(kdim+3) = iblock%cell(i,j,kdim+2)%u - iblock%cell(i,j,kdim+1)%u
                dv(kdim+3) = iblock%cell(i,j,kdim+2)%v - iblock%cell(i,j,kdim+1)%v
                dw(kdim+3) = iblock%cell(i,j,kdim+2)%w - iblock%cell(i,j,kdim+1)%w
                dp(kdim+3) = iblock%cell(i,j,kdim+2)%p - iblock%cell(i,j,kdim+1)%p
                !>
                !>
                call limit(kdim+1,max_ijk,frr,fpp,du,dv,dw,dp,dr,fl1,fl2,fl3,fl4,fl5,fr1,fr2,fr3,fr4,fr5)
                !>
                !>
                do k  = 2,kdim+1
                    !> interpolation of the primitive variables q=(r,u,v,w,p)
                    !> k0 face
                    if(iblock%cell(i,j,k-1)%wall_blank .eq. 2.0 .and. k-1 .lt. 2)then
                        !>
                        rl =iblock%cell(i,j,k-1)%r
                        ul =iblock%cell(i,j,k-1)%u
                        vl =iblock%cell(i,j,k-1)%v
                        wl =iblock%cell(i,j,k-1)%w
                        pl =iblock%cell(i,j,k-1)%p
                        !>right values
                        rr =iblock%cell( i,j,k-1)%r !- fl1(k)
                        ur =iblock%cell( i,j,k-1)%u !- fl2(k)
                        vr =iblock%cell( i,j,k-1)%v !- fl3(k)
                        wr =iblock%cell( i,j,k-1)%w !- fl4(k)
                        pr =iblock%cell( i,j,k-1)%p !- fl5(k)
                    elseif(iblock%cell(i,j,k)%wall_blank .eq. 2.0  .and.  k .gt. 2)then
                        !> kdim face
                        rl =iblock%cell(i,j,k)%r !+ fr1(k-1)
                        ul =iblock%cell(i,j,k)%u !+ fr2(k-1)
                        vl =iblock%cell(i,j,k)%v !+ fr3(k-1)
                        wl =iblock%cell(i,j,k)%w !+ fr4(k-1)
                        pl =iblock%cell(i,j,k)%p !+ fr5(k-1)
                        !>right values
                        !>
                        rr =iblock%cell(i,j,k)%r
                        ur =iblock%cell(i,j,k)%u
                        vr =iblock%cell(i,j,k)%v
                        wr =iblock%cell(i,j,k)%w
                        pr =iblock%cell(i,j,k)%p
                    else
                        !> inner grid
                        rl =iblock%cell(i,j,k-1)%r+fr1(k-1)
                        ul =iblock%cell(i,j,k-1)%u+fr2(k-1)
                        vl =iblock%cell(i,j,k-1)%v+fr3(k-1)
                        wl =iblock%cell(i,j,k-1)%w+fr4(k-1)
                        pl =iblock%cell(i,j,k-1)%p+fr5(k-1)

                        rr =iblock%cell( i,j,k)%r-fl1(k)
                        ur =iblock%cell( i,j,k)%u-fl2(k)
                        vr =iblock%cell( i,j,k)%v-fl3(k)
                        wr =iblock%cell( i,j,k)%w-fl4(k)
                        pr =iblock%cell( i,j,k)%p-fl5(k)
                    end if
!                    write(1433,'(3i4,10e20.12)')i,j,k,fl1(k-1),fl2(k-1),fl3(k-1),fl4(k-1),fl5(k-1),&
!                                                fr1(k),fr2(k),fr3(k),fr4(k),fr5(k)
                    !>
                    !>
                    !>modify the pressure and density
                    !>
                    if(pr<=0.)then
                        write(*,'("nbl=",i5," k-dir=",3i5," pr=",3e20.12)') nbl,i,j,k,pr,iblock%cell( i,j,k)%p,fr5(k)
                        pr = p0
                    end if
                    if(pl<=0.)then
                        write(*,'("nbl=",i5," k-dir=",3i5," pl=",3e20.12)') nbl,i,j,k-1,pl,iblock%cell(i,j,k-1)%p,fl5(k-1)
                        pl = p0
                    end if

                    if(rr<=0.)then
                        write(*,'("nbl=",i5," k-dir=",3i5," rr=",3e20.12)') nbl,i,j,k,rr,iblock%cell( i,j,k)%r,fr1(k)
                        rr = r0
                    end if
                    if(rl<=0.)then
                        write(*,'("nbl=",i5," k-dir=",3i5," rl=",3e20.12)') nbl,i,j,k-1,rl,iblock%cell(i,j,k-1)%r,fl1(k-1)
                        rl = r0
                    end if
                    !>upwind solution
                    !>
                    !>
!                    if(nbl ==1) write(1433,'(3i4,10e20.12)')i,j,k,rl,ul,vl,wl,pl,rr,ur,vr,wr,pr
!                    write(14333,'(3i4,4e20.12)')i,j,k,iblock%metric(i,j,k)%hi,iblock%metric(i,j,k)%hj,iblock%metric(i,j,k)%hk,iblock%metric(i,j,k)%hh
                    call roe(iblock%metric(i,j,k)%hi,iblock%metric(i,j,k)%hj,iblock%metric(i,j,k)%hk,iblock%metric(i,j,k)%hh,&
                                           rl,ul,vl,wl,pl,&
                                           rr,ur,vr,wr,pr,&
                                           f1(k),f2(k),f3(k),f4(k),f5(k))
!                   if(nbl == 1)  write(143,'(3i4,19e24.16)')i,j,k,iblock%metric(i,j,k)%hi,iblock%metric(i,j,k)%hj,iblock%metric(i,j,k)%hk,iblock%metric(i,j,k)%hh,f1(k),f2(k),f3(k),f4(k),f5(k),&
!                   rl,ul,vl,wl,pl,rr,ur,vr,wr,pr
                end do
                !>
                !>
                !> wall boundary condition
                !>
                !>
                do k=2,kdim
                    iblock%variable (i,j,k)%res_1 = iblock%variable (i,j,k)%res_1  + (f1(k+1) -   f1(k))
                    iblock%variable (i,j,k)%res_2 = iblock%variable (i,j,k)%res_2  + (f2(k+1) -   f2(k))
                    iblock%variable (i,j,k)%res_3 = iblock%variable (i,j,k)%res_3  + (f3(k+1) -   f3(k))
                    iblock%variable (i,j,k)%res_4 = iblock%variable (i,j,k)%res_4  + (f4(k+1) -   f4(k))
                    iblock%variable (i,j,k)%res_5 = iblock%variable (i,j,k)%res_5  + (f5(k+1) -   f5(k))
!                    if(nbl == 1) write(14223,'(3i4,10e24.16)')i,j,k,f1(k+1),f1(k),f2(k+1),f2(k),f3(k+1),f3(k),&
!                    f4(k+1),f4(k),f5(k+1),f5(k)
!                    if(nbl == 3) write(142233,'(3i4,5e24.16)')i,j,k,f1(k+1)-f1(k),f2(k+1)-f2(k),f3(k+1)-f3(k),&
!                    f4(k+1)-f4(k),f5(k+1)-f5(k)

                end do
            end do
        end do
        !>
        !>
        deallocate(fl1,fl2,fl3,fl4,fl5)
        deallocate(fr1,fr2,fr3,fr4,fr5)
        deallocate(frr,fpp)
        deallocate(sl1,sl2,sl3,sl4,sl5)
        deallocate(sr1,sr2,sr3,sr4,sr5)
        deallocate(f1,f2,f3,f4,f5)
        deallocate(du,dv,dw,dp,dr)
        return
    end subroutine muscl


