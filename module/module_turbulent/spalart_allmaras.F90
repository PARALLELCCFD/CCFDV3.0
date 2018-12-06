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
!>  last edit 2018-04-13
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
    !*************************************************************************!
    !                                                                         !
    !  subroutine:   spalart_allmaras.F90                                               !
    !                                                                         !
    !  programmer:   liuxz                                                    !
    !                yuanwu                                                   !
    !                                                                         !
    !  date:         2017/12/24                                               !

    !*************************************************************************!
    !>
    subroutine  spalart_allmaras(nbl)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
		implicit none
        !>
        !>
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !>
        integer ::  nbl,&
                    idim,jdim,kdim,&
                    i,j,k,&
                    ip,jp,kp,&
                    il,jl,kl,&
                    ii,imax
        real(kind = dprec) ::   factor,fact,damp_time,&
                                akarman,cb1,sigma,&
                                cb2,cw1,cw2,cw3,cv1,&
                                ct1,ct2,ct3,ct4,&
                                phi,min_precision,chi,&
                                volume,&
                                xp,yp,zp,&
                                xl,yl,zl,&
                                xa,ya,za,&
                                xc,yc,zc,&
                                ttp,ttm,cnud,cap,cam,anutp,anutm,fnup,fnum,cdp,cdm,&
                                byy,cyy,dyy,&
                                uban,sgnu,app,apm,&
                                cutoff,ss,fv1,fv2,sst,rr,gg,fix,fw,ft2,term1,term2,dist2i,tt,&
                                dfv1,dfv2,drr,dgg,dfw,dft2,&
                                volume_p,volume_l,&
                                ztmp
        !>
        !>
        real(kind = dprec),dimension(:,:,:),allocatable   :: factor_time
        real(kind = dprec),dimension(:,:,:),allocatable   :: tur_value
        real(kind = dprec),dimension(:,:,:),allocatable   :: damp
        real(kind = dprec),dimension(:,:),allocatable   :: temp
        real(kind = dprec),dimension(:,:),allocatable   :: bb
        real(kind = dprec),dimension(:,:),allocatable   :: cc
        real(kind = dprec),dimension(:,:),allocatable   :: dd
        real(kind = dprec),dimension(:,:),allocatable   :: ff
        !>
        !>
        !>
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
        if(dtt .gt. 0.0)then
            phi     = 0.5
        else
            phi     = 0.0
        end if
        !> min_precision is machine zero
        !>
        min_precision    = 1.e-15
        !>
        chi         =   1.341946
        factor      =   10.0
        akarman     =   0.41
        cb1         =   0.1355
        sigma       =   2.0/3.0
        cb2         =   0.622
        cw1         =   cb1/akarman**2+(1.0+cb2)/sigma
        cw2         =   0.30
        cw3         =   2.0
        cv1         =   7.1
        ct1         =   1.0
        ct2         =   2.0
        ct3         =   1.2
        ct4         =   0.5
        !>
        !>
        imax      = 0
        imax      = max(idim,jdim)
        imax      = max(imax,kdim)
        allocate(factor_time(idim,jdim,kdim))
        allocate(tur_value(1:idim+1,1:jdim+1,1:kdim+1))
        allocate(damp(idim,jdim,kdim))
        !>
        allocate(temp(imax,imax))
        allocate(bb(imax,imax))
        allocate(cc(imax,imax))
        allocate(dd(imax,imax))
        allocate(ff(imax,imax))
        !>
        !>
        if(dtt .lt. 0.0)then
            !> time step for turbulent model
            do k=2,kdim
                do j=2,jdim
                    do i=2,idim
                        damp_time          = factor*iblock%cell(i,j,k)%dt
                        factor_time(i,j,k) = min(damp_time,100.0)
                    end do
                end do
            end do
        else
            !>
            !> for the unsteady computing
            !> turbulence model advanced with physical time only
            do k=2,kdim
                do j=2,jdim
                    do i=2,idim
                        factor_time(i,j,k) = dtt
                    end do
                end do
            end do
        end if
        !>
        !> set up the interior points
        do i=2,idim
            do k=2,kdim
                do j=2,jdim
                    tur_value(i,j,k) = iblock%turbulent(i,j,k)%tur_save
!                    if(nbl == 93) write(600,'("inner=",e24.16)') tur_value(i,j,k)
!                    write(600+icyc,'("inner=",e24.16)')tur_value(i,j,k)
!                    write(800+icyc,'("flow=",5e24.16)')iblock%cell(i,j,k)%r,iblock%cell(i,j,k)%u,iblock%cell(i,j,k)%v,&
!                    iblock%cell(i,j,k)%w,iblock%cell(i,j,k)%p
                end do
            end do
        end do
        !>
        !>
        !> set up the boundary points
        !>
        !>set the i=i0 and i=idim boundary
        do k=2,kdim
            do j=2,jdim
                tur_value(1,j,k)      = iblock%turbulent(1,j,k)%tur_save
                tur_value(idim+1,j,k) = iblock%turbulent(idim+1,j,k)%tur_save
!                if(nbl == 93) write(500,'("i=",2e24.16)')tur_value(1,j,k),tur_value(idim+1,j,k)
            end do
        end do
        !>set the j=j0 and j=jdim boundary
        do i=2,idim
            do k=2,kdim
                tur_value(i,1,k)      = iblock%turbulent(i,1,k)%tur_save
                tur_value(i,jdim+1,k) = iblock%turbulent(i,jdim+1,k)%tur_save
!                if(nbl == 93) write(400,'("j=",2e24.16)')tur_value(i,1,k),tur_value(i,jdim+1,k)
            end do
        end do
        !>set the k=k0 and k=kdim boundary
        do i=2,idim
            do j=2,jdim
                tur_value(i,j,1)      = iblock%turbulent(i,j,1)%tur_save
                tur_value(i,j,kdim+1) = iblock%turbulent(i,j,kdim+1)%tur_save
!                if(nbl == 93) write(300,'("k=",2e24.16)')tur_value(i,j,1),tur_value(i,j,kdim+1)
            end do
        end do
        !>
        !>
        !>
        !>
        !> F_eta_eta viscous terms
        do k=3,kdim-1
            do j=2,jdim
                do i=2,idim
                    !> recomputing average metric  for the k-direction
                    !> k and k+1 average metric
                    volume_p=iblock%metric(i,j,k+1)%volume
                    !>
                    xp = iblock%metric(i,j,k+1)%hi * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    yp = iblock%metric(i,j,k+1)%hj * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    zp = iblock%metric(i,j,k+1)%hk * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    !> k-1 and k average metric
                    volume_l=iblock%metric(i,j,k-1)%volume
                    !>
                    xl = iblock%metric(i,j,k)%hi * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    yl = iblock%metric(i,j,k)%hj * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    zl = iblock%metric(i,j,k)%hk * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    !> k and k+1 average metric different the (xp,yp,zp)
                    xa = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>h
                    ya = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>
                    za = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>
                    ttp  = xp*xa + yp*ya + zp*za
                    ttm  = xl*xa + yl*ya + zl*za
                    cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                    cap  = ttp*cnud
                    cam  = ttm*cnud
                    !>
                    anutp= 0.5*(tur_value(i,j,k+1)+tur_value(i,j,k))
                    anutm= 0.5*(tur_value(i,j,k-1)+tur_value(i,j,k))
                    !>
                    fnup = 0.5*(iblock%cell(i,j,k+1)%viscous/iblock%cell(i,j,k+1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                    fnum = 0.5*(iblock%cell(i,j,k-1)%viscous/iblock%cell(i,j,k-1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                    !>
                    cdp  = (fnup+(1.+cb2)*anutp)*ttp*rmre/sigma
                    cdm  = (fnum+(1.+cb2)*anutm)*ttm*rmre/sigma
                    byy  = -max(cdm+cam,0.0)
                    cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                    dyy  = -max(cdp+cap,0.0)
                    !>
                    !>
                    iblock%turbulent(i,j,k)%viscous = -byy*tur_value(i,j,k-1) - cyy*tur_value(i,j,k) - dyy*tur_value(i,j,k+1)
!                    if(nbl ==161 .and. icyc ==2) write(790,'("axis=",3i5,7e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,&!cdp,cdm,&
!                     tur_value(i,j,k-1),tur_value(i,j,k),tur_value(i,j,k+1),&
!                    byy,cyy,dyy
                enddo
            enddo
        enddo
        !>
        !> k0 boundary points
        !>
        k=2
        kp = min(3,kdim)
        do j=2,jdim
            do i=2,idim
                !> recomputing average metric  for the k-direction
                !> k and k+1 average metric
                volume_p=iblock%metric(i,j,kp)%volume
                !>
                xp = iblock%metric(i,j,k+1)%hi * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                yp = iblock%metric(i,j,k+1)%hj * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                zp = iblock%metric(i,j,k+1)%hk * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                !> k-1 and k average metric
!                volume_l=iblock%metric(i,j,2)%volume
                volume_l=iblock%metric(i,j,1)%volume
                !>
                xl = iblock%metric(i,j,k)%hi * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                yl = iblock%metric(i,j,k)%hj * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                zl = iblock%metric(i,j,k)%hk * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                !> k and k+1 average metric different the (xp,yp,zp)
                xa = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                ya = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                za = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                ttp  = xp*xa + yp*ya + zp*za
                ttm  = xl*xa + yl*ya + zl*za
                cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                cap  = ttp*cnud
                cam  = ttm*cnud
                !>
                anutp= 0.5*(tur_value(i,j,k+1)+tur_value(i,j,k))
                anutm= 0.5*(tur_value(i,j,k-1)+tur_value(i,j,k))
                !>
                fnup = 0.5*(iblock%cell(i,j,k+1)%viscous/iblock%cell(i,j,k+1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                fnum = 0.5*(iblock%cell(i,j,k-1)%viscous/iblock%cell(i,j,k-1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                !>
                cdp  = (fnup+(1.0 + cb2)*anutp)*ttp*rmre/sigma
                cdm  = (fnum+(1.0 + cb2)*anutm)*ttm*rmre/sigma
                byy  = -max(cdm+cam,0.0)
                cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                dyy  = -max(cdp+cap,0.0)
                !>
                !>
                iblock%turbulent(i,j,k)%viscous = -byy*tur_value(i,j,k-1) - cyy*tur_value(i,j,k) - dyy*tur_value(i,j,k+1)
!                if(nbl ==161 .and. icyc ==2) write(80001,'("axis=",3i5,7e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,byy,cyy,dyy,tur_value(i,j,k-1),tur_value(i,j,k),tur_value(i,j,k+1)
            end do
        end do
        !>
        !> kdim boundary points
        !>
        k  = kdim
        kl = kdim-1
        do j=2,jdim
            do i=2,idim
                !> recomputing average metric  for the k-direction
                !> j and j+1 average metric
                volume_p=iblock%metric(i,j,kdim+1)%volume
                !>
                xp = iblock%metric(i,j,kdim+1)%hi * iblock%metric(i,j,kdim+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                yp = iblock%metric(i,j,kdim+1)%hj * iblock%metric(i,j,kdim+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                zp = iblock%metric(i,j,kdim+1)%hk * iblock%metric(i,j,kdim+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                !> k-1 and k average metric
                volume_l=iblock%metric(i,j,kdim-1)%volume
                !>
                xl = iblock%metric(i,j,k)%hi * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                yl = iblock%metric(i,j,k)%hj * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                zl = iblock%metric(i,j,k)%hk * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                !> j and j+1 average metric different the (xp,yp,zp)
                xa = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                ya = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                za = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                ttp  = xp*xa + yp*ya + zp*za
                ttm  = xl*xa + yl*ya + zl*za
                cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                cap  = ttp*cnud
                cam  = ttm*cnud
                !>
                anutp= 0.5*(tur_value(i,j,k+1)+tur_value(i,j,k))
                anutm= 0.5*(tur_value(i,j,k-1)+tur_value(i,j,k))
                !>
                fnup = 0.5*(iblock%cell(i,j,k+1)%viscous/iblock%cell(i,j,k+1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                fnum = 0.5*(iblock%cell(i,j,k-1)%viscous/iblock%cell(i,j,k+1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                !>
                cdp  = (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                cdm  = (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                byy  = -max(cdm+cam,0.0)
                cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                dyy  = -max(cdp+cap,0.0)
                !>
                !>
                iblock%turbulent(i,j,k)%viscous = -byy*tur_value(i,j,k-1) - cyy*tur_value(i,j,k) - dyy*tur_value(i,j,k+1)
!                if(nbl ==161 .and. icyc ==2) write(80002,'("axis=",3i5,7e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,byy,cyy,dyy,tur_value(i,j,k-1),tur_value(i,j,k),tur_value(i,j,k+1)
            end do
        end do
!        if(nbl == 7)write(*,'("run at sa viscous calculate #2")')
        !>
        !> advective terms in eta
        !>
!        open(unit=1807023,file='turbulent_viscous_eta.dat',status='unknown')
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    !>
                    !>
                    xc = 0.5*( iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh )/iblock%metric(i,j,k)%volume
                    !>
                    yc = 0.5*( iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh )/iblock%metric(i,j,k)%volume
                    !>
                    zc = 0.5*( iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh )/iblock%metric(i,j,k)%volume
                    !> the time term
                    !>
                    !>tc = 0.5*( iblock%metric(i+1,j,k)%ft*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff )/iblock%metric(i,j,k)%volume
                    !>
                    uban = xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v +zc*iblock%cell(i,j,k)%w !+ tc
                    sgnu = sign(1.0 , uban)
                    app  = 0.5*(1.0 + sgnu)
                    apm  = 0.5*(1.0 - sgnu)
!                   cdp  = iblock%turbulent(i,j,k)%viscous
                    iblock%turbulent(i,j,k)%viscous = iblock%turbulent(i,j,k)%viscous - &
                                                      uban * (app * (tur_value(i,j,k  )-tur_value(i,j,k-1)) + &
                                                              apm * (tur_value(i,j,k+1)-tur_value(i,j,k  )))
!                    if(nbl ==161 .and. icyc ==2) write(800,'("axis=",3i5,4e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,uban,cdp,sgnu
                end do
            end do
        end do
        !>
        !>
        !>
        !>
        !>
        !> F_xi_xi viscous terms
        do j=3,jdim-1
            do k=2,kdim
                do i=2,idim
                    !> recomputing average metric  for the j-direction
                    !> j and j+1 average metric
                    volume_p=iblock%metric(i,j+1,k)%volume
                    !>
                    xp = iblock%metric(i,j+1,k)%gi * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    yp = iblock%metric(i,j+1,k)%gj * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    zp = iblock%metric(i,j+1,k)%gk * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    !> j-1 and j average metric
                    volume_l=iblock%metric(i,j-1,k)%volume
                    !>
                    xl = iblock%metric(i,j,k)%gi * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    yl = iblock%metric(i,j,k)%gj * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    zl = iblock%metric(i,j,k)%gk * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    !> j and j+1 average metric different the (xp,yp,zp)
                    xa = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    ya = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    za = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    ttp  = xp*xa + yp*ya + zp*za
                    ttm  = xl*xa + yl*ya + zl*za
                    cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                    cap  = ttp*cnud
                    cam  = ttm*cnud
                    !>
                    anutp= 0.5*(tur_value(i,j+1,k)+tur_value(i,j,k))
                    anutm= 0.5*(tur_value(i,j-1,k)+tur_value(i,j,k))
                    !>
                    fnup = 0.5*(iblock%cell(i,j+1,k)%viscous/iblock%cell(i,j+1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                    fnum = 0.5*(iblock%cell(i,j-1,k)%viscous/iblock%cell(i,j-1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                    !>
                    cdp  = (fnup+(1.+cb2)*anutp)*ttp*rmre/sigma
                    cdm  = (fnum+(1.+cb2)*anutm)*ttm*rmre/sigma
                    byy  = -max(cdm+cam,0.0)
                    cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                    dyy  = -max(cdp+cap,0.0)
                    !>
                    !>
                    iblock%turbulent(i,j,k)%viscous = iblock%turbulent(i,j,k)%viscous -byy*tur_value(i,j-1,k) - cyy*tur_value(i,j,k) - dyy*tur_value(i,j+1,k)
                enddo
            enddo
        enddo
!        if(nbl == 7)write(*,'("run at sa viscous calculate #3")')
        !>
        !> j0 boundary points
        !>
        j=2
        jp = min(3,jdim)
        do k=2,kdim
            do i=2,idim
                !> recomputing average metric  for the j-direction
                !> j and j+1 average metric
                volume_p=iblock%metric(i,jp,k)%volume
                !>
                xp = iblock%metric(i,j+1,k)%gi * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                yp = iblock%metric(i,j+1,k)%gj * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                zp = iblock%metric(i,j+1,k)%gk * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                !> j-1 and j average metric
!                volume_l=iblock%metric(i,2,k)%volume
                volume_l=iblock%metric(i,1,k)%volume
                !>
                xl = iblock%metric(i,j,k)%gi * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                yl = iblock%metric(i,j,k)%gj * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                zl = iblock%metric(i,j,k)%gk * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                !> j and j+1 average metric different the (xp,yp,zp)
                xa = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                ya = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                za = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                ttp  = xp*xa + yp*ya + zp*za
                ttm  = xl*xa + yl*ya + zl*za
                cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                cap  = ttp*cnud
                cam  = ttm*cnud
                !>
                anutp= 0.5*(tur_value(i,j+1,k)+tur_value(i,j,k))
                anutm= 0.5*(tur_value(i,j-1,k)+tur_value(i,j,k))
                !>
                fnup = 0.5*(iblock%cell(i,j+1,k)%viscous/iblock%cell(i,j+1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                fnum = 0.5*(iblock%cell(i,j-1,k)%viscous/iblock%cell(i,j-1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                !>
                cdp  = (fnup+(1.0 + cb2)*anutp)*ttp*rmre/sigma
                cdm  = (fnum+(1.0 + cb2)*anutm)*ttm*rmre/sigma
                byy  = -max(cdm+cam,0.0)
                cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                dyy  = -max(cdp+cap,0.0)
                !>
                !>
!                write(18070233,'("axis=",3i5,7e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,byy,cyy,dyy,tur_value(i,j-1,k),tur_value(i,j,k),tur_value(i,j+1,k)
                iblock%turbulent(i,j,k)%viscous = iblock%turbulent(i,j,k)%viscous -byy*tur_value(i,j-1,k) - cyy*tur_value(i,j,k) - dyy*tur_value(i,j+1,k)

            end do
        end do
        !>
        !> jdim boundary points
        !>
        j  = jdim
        jl = jdim-1
        do k=2,kdim
            do i=2,idim
                !> recomputing average metric  for the j-direction
                !> j and j+1 average metric
                volume_p=iblock%metric(i,jdim+1,k)%volume
                !>
                xp = iblock%metric(i,jdim+1,k)%gi * iblock%metric(i,jdim+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                yp = iblock%metric(i,jdim+1,k)%gj * iblock%metric(i,jdim+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                zp = iblock%metric(i,jdim+1,k)%gk * iblock%metric(i,jdim+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                !> j-1 and j average metric
                volume_l=iblock%metric(i,jdim-1,k)%volume
                !>
                xl = iblock%metric(i,j,k)%gi * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                yl = iblock%metric(i,j,k)%gj * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                zl = iblock%metric(i,j,k)%gk * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                !> j and j+1 average metric different the (xp,yp,zp)
                xa = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                ya = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                za = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                ttp  = xp*xa + yp*ya + zp*za
                ttm  = xl*xa + yl*ya + zl*za
                cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                cap  = ttp*cnud
                cam  = ttm*cnud
                !>
                anutp= 0.5*(tur_value(i,j+1,k)+tur_value(i,j,k))
                anutm= 0.5*(tur_value(i,j-1,k)+tur_value(i,j,k))
                !>
                fnup = 0.5*(iblock%cell(i,j+1,k)%viscous/iblock%cell(i,j+1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                fnum = 0.5*(iblock%cell(i,j-1,k)%viscous/iblock%cell(i,j-1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                !>
                cdp  = (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                cdm  = (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                byy  = -max(cdm+cam,0.0)
                cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                dyy  = -max(cdp+cap,0.0)
                !>
                !>
!                write(18070234,'("axis=",3i5,7e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,byy,cyy,dyy,tur_value(i,j-1,k),tur_value(i,j,k),tur_value(i,j+1,k)
                iblock%turbulent(i,j,k)%viscous = iblock%turbulent(i,j,k)%viscous - byy*tur_value(i,j-1,k) - cyy*tur_value(i,j,k) - dyy*tur_value(i,j+1,k)

            end do
        end do
!        if(nbl == 7)write(*,'("run at sa viscous calculate #4")')
        !>
        !> advective terms in xi
        !>
!        open(unit=1807024,file='turbulent_viscous_xi.dat',status='unknown')
        do i=2,idim
            do k=2,kdim
                do j=2,jdim
                    !>
                    !>
                    xc = 0.5*( iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg )/iblock%metric(i,j,k)%volume
                    !>
                    yc = 0.5*( iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg )/iblock%metric(i,j,k)%volume
                    !>
                    zc = 0.5*( iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg )/iblock%metric(i,j,k)%volume
                    !> the time term
                    !>
                    !>tc = 0.5*( iblock%metric(i+1,j,k)%ft*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff )/iblock%metric(i,j,k)%volume
                    !>
                    uban = xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v +zc*iblock%cell(i,j,k)%w !+ tc
                    sgnu = sign(1.0 , uban)
                    app  = 0.5*(1.0 + sgnu)
                    apm  = 0.5*(1.0 - sgnu)
!                    if(nbl==1) write(18070241,'("sgnu=",3i4,6e24.16)') i,j,k,xc,yc,zc,iblock%cell(i,j,k)%u,iblock%cell(i,j,k)%v,iblock%cell(i,j,k)%w
                    cdm  = iblock%turbulent(i,j,k)%viscous
                    iblock%turbulent(i,j,k)%viscous = iblock%turbulent(i,j,k)%viscous - &
                                                      uban * (app * (tur_value(i,j  ,k)-tur_value(i,j-1,k)) + &
                                                              apm * (tur_value(i,j+1,k)-tur_value(i,j  ,k)))
!                    if(nbl ==161 .and. icyc ==2) write(810,'("axis=",3i5,4e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,uban,cdm,sgnu
                end do
            end do
        end do
        !>
        !>
        !>
!        if(nbl == 7)write(*,'("run at sa viscous calculate #5")')
        !> F_zeta_zeta viscous terms
        if(idim .gt. 3)then
            do i=3,idim-1
                do k=2,kdim
                    do j=2,jdim
                        !> recomputing average metric  for the i-direction
                        !> i and i+1 average metric
                        volume_p=iblock%metric(i+1,j,k)%volume
                        !>
                        xp = iblock%metric(i+1,j,k)%fi * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                        !>
                        yp = iblock%metric(i+1,j,k)%fj * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                        !>
                        zp = iblock%metric(i+1,j,k)%fk * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                        !>
                        !> i-1 and i average metric
                        volume_l=iblock%metric(i-1,j,k)%volume
                        !>
                        xl = iblock%metric(i,j,k)%fi * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                        !>
                        yl = iblock%metric(i,j,k)%fj * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                        !>
                        zl = iblock%metric(i,j,k)%fk * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                        !>
                        !> i and i+1 average metric different the (xp,yp,zp)
                        xa = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        ya = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        za = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        ttp  = xp*xa + yp*ya + zp*za
                        ttm  = xl*xa + yl*ya + zl*za
                        cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                        cap  = ttp*cnud
                        cam  = ttm*cnud
                        !>
                        anutp= 0.5*(tur_value(i+1,j,k)+tur_value(i,j,k))
                        anutm= 0.5*(tur_value(i-1,j,k)+tur_value(i,j,k))
                        !>
                        fnup = 0.5*(iblock%cell(i+1,j,k)%viscous/iblock%cell(i+1,j,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                        fnum = 0.5*(iblock%cell(i-1,j,k)%viscous/iblock%cell(i-1,j,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                        !>
                        cdp  = (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                        cdm  = (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                        byy  = -max(cdm+cam,0.0)
                        cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                        dyy  = -max(cdp+cap,0.0)
                        !>
                        !>
                        iblock%turbulent(i,j,k)%viscous = iblock%turbulent(i,j,k)%viscous -byy*tur_value(i-1,j,k) - cyy*tur_value(i,j,k) - dyy*tur_value(i+1,j,k)
                    enddo
                enddo
            enddo
            !>
            !> i0 boundary points
            !>
            i=2
            ip = min(3,idim)
            do k=2,kdim
                do j=2,jdim
                    !> recomputing average metric  for the i-direction
                    !> i and i+1 average metric
                    volume_p=iblock%metric(ip,j,k)%volume
                    !>
                    xp = iblock%metric(i+1,j,k)%fi * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    yp = iblock%metric(i+1,j,k)%fj * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    zp = iblock%metric(i+1,j,k)%fk * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    !> i-1 and i average metric
!                    volume_l=iblock%metric(2,j,k)%volume
                    volume_l=iblock%metric(1,j,k)%volume
                    !>
                    xl = iblock%metric(i,j,k)%fi * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    yl = iblock%metric(i,j,k)%fj * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    zl = iblock%metric(i,j,k)%fk * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    !> i and i+1 average metric different the (xp,yp,zp)
                    xa = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    ya = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    za = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    ttp  = xp*xa + yp*ya + zp*za
                    ttm  = xl*xa + yl*ya + zl*za
                    cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                    cap  = ttp*cnud
                    cam  = ttm*cnud
                    !>
                    anutp= 0.5*(tur_value(i+1,j,k)+tur_value(i,j,k))
                    anutm= 0.5*(tur_value(i-1,j,k)+tur_value(i,j,k))
                    !>
                    fnup = 0.5*(iblock%cell(i+1,j,k)%viscous/iblock%cell(i+1,j,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                    fnum = 0.5*(iblock%cell(i-1,j,k)%viscous/iblock%cell(i-1,j,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                    !>
                    cdp  = (fnup+(1.+cb2)*anutp)*ttp*rmre/sigma
                    cdm  = (fnum+(1.+cb2)*anutm)*ttm*rmre/sigma
                    byy  = -max(cdm+cam,0.0)
                    cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                    dyy  = -max(cdp+cap,0.0)
                    !>
                    !>
!                    write(18070235,'("axis=",3i5,7e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,byy,cyy,dyy,tur_value(i-1,j,k),tur_value(i,j,k),tur_value(i+1,j,k)
                    iblock%turbulent(i,j,k)%viscous = iblock%turbulent(i,j,k)%viscous -byy*tur_value(i-1,j,k) - cyy*tur_value(i,j,k) - dyy*tur_value(i+1,j,k)
                end do
            end do
            !>
            !> idim boundary points
            !>
            i  = idim
            il = idim-1
            do k=2,kdim
                do j=2,jdim
                    !> recomputing average metric  for the i-direction
                    !> i and i+1 average metric
                    volume_p=iblock%metric(idim+1,j,k)%volume
                    !>
                    xp = iblock%metric(idim+1,j,k)%fi * iblock%metric(idim+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    yp = iblock%metric(idim+1,j,k)%fj * iblock%metric(idim+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    zp = iblock%metric(idim+1,j,k)%fk * iblock%metric(idim+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    !> i-1 and i average metric
                    volume_l=iblock%metric(idim-1,j,k)%volume
                    !>
                    xl = iblock%metric(i,j,k)%fi * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    yl = iblock%metric(i,j,k)%fj * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    zl = iblock%metric(i,j,k)%fk * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    !> i and i+1 average metric different the (xp,yp,zp)
                    xa = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    ya = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    za = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    ttp  = xp*xa + yp*ya + zp*za
                    ttm  = xl*xa + yl*ya + zl*za
                    cnud = -cb2*tur_value(i,j,k)*rmre/sigma
                    cap  = ttp*cnud
                    cam  = ttm*cnud
                    !>
                    anutp= 0.5*(tur_value(i+1,j,k)+tur_value(i,j,k))
                    anutm= 0.5*(tur_value(i-1,j,k)+tur_value(i,j,k))
                    !>
                    fnup = 0.5*(iblock%cell(i+1,j,k)%viscous/iblock%cell(i+1,j,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                    fnum = 0.5*(iblock%cell(i-1,j,k)%viscous/iblock%cell(i-1,j,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r )
                    !>
                    cdp  = (fnup+(1.+cb2)*anutp)*ttp*rmre/sigma
                    cdm  = (fnum+(1.+cb2)*anutm)*ttm*rmre/sigma
                    byy  = -max(cdm+cam,0.0)
                    cyy  =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                    dyy  = -max(cdp+cap,0.0)
                    !>
                    !>
!                    write(18070236,'("axis=",3i5,7e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,byy,cyy,dyy,tur_value(i-1,j,k),tur_value(i,j,k),tur_value(i+1,j,k)
                    iblock%turbulent(i,j,k)%viscous =iblock%turbulent(i,j,k)%viscous -byy*tur_value(i-1,j,k) - cyy*tur_value(i,j,k) - dyy*tur_value(i+1,j,k)
                end do
            end do
            !>
            !> advective terms in zeta
            !>
!            open(unit=1807025,file='turbulent_viscous_zeta.dat',status='unknown')
            do i=2,idim
                do k=2,kdim
                    do j=2,jdim
                        !>
                        !>
                        xc = 0.5*( iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff )/iblock%metric(i,j,k)%volume
                        !>
                        yc = 0.5*( iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff )/iblock%metric(i,j,k)%volume
                        !>
                        zc = 0.5*( iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff )/iblock%metric(i,j,k)%volume
                        !> the time term
                        !>
                        !>tc = 0.5*( iblock%metric(i+1,j,k)%ft*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff )/iblock%metric(i,j,k)%volume
                        !>
                        uban = xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v +zc*iblock%cell(i,j,k)%w !+ tc
                        sgnu = sign(1.0 , uban)
                        app  = 0.5*(1.0 + sgnu)
                        apm  = 0.5*(1.0 - sgnu)
                        iblock%turbulent(i,j,k)%viscous = iblock%turbulent(i,j,k)%viscous - &
                                                          uban * (app * (tur_value(i  ,j,k)-tur_value(i-1,j,k)) + &
                                                                  apm * (tur_value(i+1,j,k)-tur_value(i  ,j,k)))
!                        if(nbl ==161 .and. icyc ==2 ) write(820,'("axis=",3i5,5e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,xc,yc,zc,uban
                    end do
                end do
            end do
        end if
        !>
        !>
        !>
        !>
        !>
        !>
        !> use the damp to temporarily the distance variable
        !> later the damp used to store LHS implicit contribution
        !>
!        if(nbl == 7)write(*,'("run at sa viscous calculate #6")')
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    !damp(i,j,k)  =  0.1186612621696583E+02!abs(iblock%turbulent(i,j,k)%distance)
                    damp(i,j,k)  =  abs(iblock%turbulent(i,j,k)%distance)
!                    if(icyc==1 .and. nbl==161) write(9000+nbl,'("axis=",3i5,e24.16)') i,j,k,damp(i,j,k)
!                    if(nbl == 93) write(700,'("viscous=",2e24.16)') iblock%turbulent(i,j,k)%viscous,damp(i,j,k)
!                    write(1807232,'("damp=",3i4,e24.16)') i-1,j-1,k-1,damp(i,j,k)
                end do
            end do
        end do
        !>
        !>
!        open(unit=1807031,file='turbulent_viscous_damp.dat',status='unknown')
        do i=2,idim
            do j=2,jdim
                do k=2,kdim
                    if(iblock%turbulent(i,j,k)%distance .lt. 0.0)then
                        cutoff = 0.0
                    else
                        cutoff = 1.0
                    end if
                    !>
                    cutoff = 1.0
                    !>
                    ss  = iblock%turbulent(i,j,k)%vorticity
                    chi = tur_value(i,j,k)/(iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    !>
                    fv1 = chi**3/(chi**3 +cv1**3)
                    fv2 = 1.0-(chi/(1.0 +chi*fv1))
                    !>
                    sst = ss +tur_value(i,j,k)*fv2/((1.0/rmre)*(akarman*damp(i,j,k))**2)
                    sst = max(sst,min_precision)
                    rr  = tur_value(i,j,k)/((1.0/rmre)*sst*(akarman*damp(i,j,k))**2)
                    rr  = min(rr,10.0)
                    gg  = rr+cw2*(rr**6-rr)
                    gg  = max(gg,min_precision)
                    !>
                    !> fix for single precision.
                    !>
                    fw    = gg*((1.0+cw3**6)/(gg**6+cw3**6))**(1./6.)
                    ft2   = ct3*exp(-ct4*chi**2)
                    term1 = cb1*(1.0-ft2)*ss
                    term2 = cb1*((1.0-ft2)*fv2+ft2) / akarman**2-cw1*fw
                    dist2i= 1.0/((1.0/rmre)*damp(i,j,k)**2   + 1.e-20)
                    tt    = cutoff*term1*tur_value(i,j,k) + term2*tur_value(i,j,k)**2*dist2i
                    !>  Store quantity to be added to certain implicit LHS terms:
                    damp(i,j,k) = 2.0 * term2*tur_value(i,j,k) * dist2i
                    dfv1        = (fv1-fv1*fv1)*3.0/tur_value(i,j,k)
                    dfv2        = (fv2-1.0)/tur_value(i,j,k)+((1.0-fv2)**2)*(fv1/tur_value(i,j,k)+dfv1)
                    dft2        = -(2.0*ct4*tur_value(i,j,k)/((iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)**2))*ft2
                    drr         = rr/tur_value(i,j,k)-rr*rr*(fv2/tur_value(i,j,k)+dfv2)
                    dgg         = (1.0-cw2+6.*cw2*(rr**5))*drr
                    gg          = max(gg,min_precision*10.0)
                    !> fix for single precision.
                    dfw         = ((1.0 + cw3**6)/( gg**6 + cw3**6 ))**(1.0/6.0) - ((1.0 + cw3**6 )**(1.0/6.0)/((gg**6 + cw3**6 )**(7.0/6.0)))*gg**6
                    dfw         = dfw*dgg
                    damp(i,j,k) = damp(i,j,k)+dist2i*(tur_value(i,j,k)**2)*(cb1/(akarman**2)*(dfv2-ft2*dfv2-fv2*dft2+dft2)-cw1*dfw)
                    !> add to RHS:
                    !write(1807031,'(3i4,7e24.16)') i-1,j-1,k-1,tt,ss,damp(i,j,k),dist2i,term1,term2,chi
                    iblock%turbulent(i,j,k)%viscous  =  iblock%turbulent(i,j,k)%viscous  +  tt
!                    if(nbl ==161 .and. icyc==2)write(830,'(3i4,4e24.16)') i-1,j-1,k-1,iblock%turbulent(i,j,k)%viscous,damp(i,j,k),tt,ss
                end do
            end do
        end do
        !>
        !>
        !>
        !>
        !> implicit F_eta_eta viscous terms
        !>
        !>
!        if(nbl == 7)write(*,'("run at sa viscous calculate #7")')
        do i=2,idim
            !>
            !> interior points
            do k=3,kdim-1
                do j=2,jdim
                    !>
                    volume_p=iblock%metric(i,j,k+1)%volume
                    !> k and k+1 average metric
                    xp = iblock%metric(i,j,k+1)%hi * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    yp = iblock%metric(i,j,k+1)%hj * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    zp = iblock%metric(i,j,k+1)%hk * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    !> k-1 and k average metric
                    volume_l=iblock%metric(i,j,k-1)%volume
                    !>
                    xl = iblock%metric(i,j,k)%hi * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    yl = iblock%metric(i,j,k)%hj * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    zl = iblock%metric(i,j,k)%hk * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    !> k and k+1 average metric different the (xp,yp,zp)
                    xa = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>
                    ya = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>
                    za = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>
                    ttp =   xp*xa+yp*ya+zp*za
                    ttm =   xl*xa+yl*ya+zl*za
                    !>
                    cnud    = -cb2*tur_value(i,j,k)*rmre/sigma
                    cap     = ttp*cnud
                    cam     = ttm*cnud
                    anutp   = 0.5*(tur_value(i,j,k+1)+tur_value(i,j,k))
                    anutm   = 0.5*(tur_value(i,j,k-1)+tur_value(i,j,k))
                    fnup    = 0.5*(iblock%cell(i,j,k+1)%viscous/iblock%cell(i,j,k+1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    fnum    = 0.5*(iblock%cell(i,j,k-1)%viscous/iblock%cell(i,j,k-1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    cdp     = (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                    cdm     = (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                    !>
                    !>
                    bb(j,k) = -max(cdm+cam,0.0)
                    cc(j,k) =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                    dd(j,k) = -max(cdp+cap,0.0)
                enddo
                !>
                !>
                do j=2,jdim
                    !>
                    xc = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>
                    yc = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>
                    zc = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                    !>
                    !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    uban = xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                    sgnu = sign(1.0,uban)
                    app  = 0.5*(1.0+sgnu)
                    apm  = 0.5*(1.0-sgnu)
                    bb(j,k)=bb(j,k) - uban*app
                    !> Add part of source terms to diagonal LHS in this sweep
                    if(real(damp(i,j,k)) .lt. 0.) then
                        cc(j,k)=cc(j,k) + uban*(app-apm) - damp(i,j,k)
                    else
                        cc(j,k)=cc(j,k) + uban*(app-apm)
                    end if
                    dd(j,k)=dd(j,k) + uban*apm
                enddo
                !>
                !>
                do j=2,jdim
                    fact    = factor_time(i,j,k)
                    bb(j,k) = bb(j,k)*fact
                    cc(j,k) = cc(j,k)*fact+1.0*(1.0+phi)
                    dd(j,k) = dd(j,k)*fact
                    ff(j,k) = iblock%turbulent(i,j,k)%viscous*fact
!                    write(18070321,'(3i4,5e24.16)') i-1,j-1,k-1,ff(j,k),fact,bb(j,k),cc(j,k),dd(j,k)
                enddo
            enddo
            !>
            !> k0 boundary points
            k  = 2
            kp = min(3,kdim)
            do j=2,jdim
                volume_p = iblock%metric(i,j,kp)%volume
                !> k and k+1 average metric
                xp = iblock%metric(i,j,k+1)%hi * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                yp = iblock%metric(i,j,k+1)%hj * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                zp = iblock%metric(i,j,k+1)%hk * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                !> k-1 and k average metric
!                volume_l=iblock%metric(i,j,2)%volume
                volume_l=iblock%metric(i,j,1)%volume
                !>
                xl = iblock%metric(i,j,k)%hi * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                yl = iblock%metric(i,j,k)%hj * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                zl = iblock%metric(i,j,k)%hk * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                !> k and k+1 average metric different the (xp,yp,zp)
                xa = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                ya = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                za = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                !>
                ttp  = xp*xa+yp*ya+zp*za
                ttm  = xl*xa+yl*ya+zl*za
                !>
                cnud    =   -cb2*tur_value(i,j,k)*rmre/sigma
                cap     =   ttp*cnud
                cam     =   ttm*cnud
                !>
                anutp   =   0.5*(tur_value(i,j,k+1)+tur_value(i,j,k))
                anutm   =   0.5*(tur_value(i,j,k-1)+tur_value(i,j,k))
                !>
                fnup    =   0.5*(iblock%cell(i,j,k+1)%viscous/iblock%cell(i,j,k+1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                fnum    =   0.5*(iblock%cell(i,j,k-1)%viscous/iblock%cell(i,j,k-1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                !>
                cdp     =   (fnup+(1.0 + cb2)*anutp)*ttp*rmre/sigma
                cdm     =   (fnum+(1.0 + cb2)*anutm)*ttm*rmre/sigma
                bb(j,k) =   -max(cdm+cam,0.0)
                cc(j,k) =    max(cdp+cap,0.0) + max(cdm+cam,0.0)
                dd(j,k) =   -max(cdp+cap,0.0)
            enddo
            !>
            !>
            do j=2,jdim
                !>
                xc = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                yc = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                zc = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                !>
                uban    =   xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                sgnu    =   sign(1.0,uban)
                app     =   0.5*(1.+sgnu)
                apm     =   0.5*(1.-sgnu)
                bb(j,k) =   bb(j,k) - uban*app
                !> add part of source terms to diagonal LHS in this sweep
                if(real(damp(i,j,k)) .lt. 0.) then
                    cc(j,k)=cc(j,k) + uban*(app-apm) - damp(i,j,k)
                else
                    cc(j,k)=cc(j,k) + uban*(app-apm)
                end if
                dd(j,k) =   dd(j,k) + uban*apm
            enddo
            !>
            !>
            do j=2,jdim
                fact    =   factor_time(i,j,k)
                bb(j,k) =   bb(j,k)*fact
                cc(j,k) =   cc(j,k)*fact+1.0*(1.0+phi)
                dd(j,k) =   dd(j,k)*fact
                ff(j,k) =   iblock%turbulent(i,j,k)%viscous*fact
!                write(18070322,'(3i4,5e24.16)') i-1,j-1,k-1,ff(j,k),fact,bb(j,k),cc(j,k),dd(j,k)
            enddo
            !>
            !> kdim boundary points
            k = kdim
            !>
            do j=2,jdim
                volume_p=iblock%metric(i,j,kdim+1)%volume
                !>
                xp = iblock%metric(i,j,k+1)%hi * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                yp = iblock%metric(i,j,k+1)%hj * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                zp = iblock%metric(i,j,k+1)%hk * iblock%metric(i,j,k+1)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                !> k-1 and k average metric
                volume_l=iblock%metric(i,j,kdim-1)%volume
                !>
                xl = iblock%metric(i,j,k)%hi * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                yl = iblock%metric(i,j,k)%hj * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                zl = iblock%metric(i,j,k)%hk * iblock%metric(i,j,k)%hh / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                !> k and k+1 average metric different the (xp,yp,zp)
                xa = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                ya = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                za = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh)/iblock%metric(i,j,k)%volume
                !>
                ttp=xp*xa+yp*ya+zp*za
                ttm=xl*xa+yl*ya+zl*za
                !>
                cnud    =   -cb2*tur_value(i,j,k)*rmre/sigma
                cap     =   ttp*cnud
                cam     =   ttm*cnud
                anutp   =   0.5*(tur_value(i,j,k+1)+ tur_value(i,j,k))
                anutm   =   0.5*(tur_value(i,j,k-1) + tur_value(i,j,k))
                !>
                fnup    =   0.5*(iblock%cell(i,j,k+1)%viscous/iblock%cell(i,j,k+1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                fnum    =   0.5*(iblock%cell(i,j,k-1)%viscous/iblock%cell(i,j,k-1)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                cdp     =   (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                cdm     =   (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                bb(j,k) =   -max(cdm+cam,0.0)
                cc(j,k) =    max(cdp+cap,0.0) + max(cdm+cam,0.0)
                dd(j,k) =   -max(cdp+cap,0.0)
!                write(1807032301,'(3i4,3e24.16)') i-1,j-1,k-1,bb(j,k),cc(j,k),dd(j,k)
            enddo
            !>
            !>
            do j=2,jdim
                xc = 0.5*(iblock%metric(i,j,k+1)%hi*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh )/iblock%metric(i,j,k)%volume
                !>
                yc = 0.5*(iblock%metric(i,j,k+1)%hj*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh )/iblock%metric(i,j,k)%volume
                !>
                zc = 0.5*(iblock%metric(i,j,k+1)%hk*iblock%metric(i,j,k+1)%hh + iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh )/iblock%metric(i,j,k)%volume
                !>
                !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                !>
                uban    =   xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                sgnu    =   sign(1.0,uban)
                app     =   0.5*(1.+sgnu)
                apm     =   0.5*(1.-sgnu)
                bb(j,k) =   bb(j,k) - uban*app
                !> add part of source terms to diagonal LHS in this sweep
                if(real(damp(i,j,k)) .lt. 0.) then
                    cc(j,k) = cc(j,k) + uban*(app-apm) - damp(i,j,k)
                else
                    cc(j,k) = cc(j,k) + uban*(app-apm)
                end if
                dd(j,k) =   dd(j,k) + uban*apm
!                write(1807032302,'(3i4,7e24.16)') i-1,j-1,k-1,bb(j,k),cc(j,k),dd(j,k),uban,damp(i,j,k),app,apm
            enddo
            !>
            !>
            do j=2,jdim
                fact    =   factor_time(i,j,k)
                bb(j,k) =   bb(j,k)*fact
                cc(j,k) =   cc(j,k)*fact+1.0*(1.0+phi)
                dd(j,k) =   dd(j,k)*fact
                ff(j,k) =   iblock%turbulent(i,j,k)%viscous*fact
!                write(18070323,'(3i4,5e24.16)') i-1,j-1,k-1,ff(j,k),fact,bb(j,k),cc(j,k),dd(j,k)
            enddo
            !>
            !>
            if (overlap .eq. 1) then
                do k=2,kdim
                    do j=2,jdim
                        ff(j,k) =   ff(j,k)*iblock%cell(i,j,k)%blank
                        bb(j,k) =   bb(j,k)*iblock%cell(i,j,k)%blank
                        dd(j,k) =   dd(j,k)*iblock%cell(i,j,k)%blank
                        cc(j,k) =   cc(j,k)*iblock%cell(i,j,k)%blank+( 1.0 - iblock%cell(i,j,k)%blank )
                    enddo
                enddo
            end if
            !>
            !> solve a scalar tridiagonal system of equations
            do j=2,jdim
                temp(j,2) = dd(j,2)/cc(j,2)
                ff(j,2)   = ff(j,2)/cc(j,2)
            end do
            !>
            do k=3,kdim
                do j=2,jdim
                    ztmp       = 1.0/(cc(j,k)-bb(j,k)*temp(j,k-1))
                    temp(j,k) = dd(j,k)*ztmp
                    ff(j,k) = (ff(j,k)-bb(j,k)*ff(j,k-1))*ztmp
                   end do
            end do
            !>
            do ii=3,kdim
                k=kdim + 2 - ii
                do j=2,jdim
                    ff(j,k)=ff(j,k)-temp(j,k)*ff(j,k+1)
                end do
            end do
            !>
            !>
            do k=2,kdim
                do j=2,jdim
                    iblock%turbulent(i,j,k)%viscous =   ff(j,k)
!                    if(nbl ==161 .and. icyc==2) write(840,'("axis=",3i5,e24.16)') i,j,k,ff(j,k)
                enddo
            enddo
        enddo
!        if(nbl == 7)write(*,'("run at sa viscous calculate #8")')
        !>
        !> implicit F_xi_xi viscous terms
        !>
        !>
        !>
!        open(unit=1807033,file='turbulent_viscous_X_eta.dat',status='unknown')
        do i=2,idim
            !>
            !> interior points
            do j=3,jdim-1
                do k=2,kdim
                    !>
                    volume_p=iblock%metric(i,j+1,k)%volume
                    !> j and j+1 average metric
                    xp = iblock%metric(i,j+1,k)%gi * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    yp = iblock%metric(i,j+1,k)%gj * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    zp = iblock%metric(i,j+1,k)%gk * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    !> j-1 and j average metric
                    volume_l=iblock%metric(i,j-1,k)%volume
                    !>
                    xl = iblock%metric(i,j,k)%gi * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    yl = iblock%metric(i,j,k)%gj * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    zl = iblock%metric(i,j,k)%gk * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    !> j and j+1 average metric different the (xp,yp,zp)
                    xa = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    ya = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    za = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    ttp =   xp*xa+yp*ya+zp*za
                    ttm =   xl*xa+yl*ya+zl*za
                    !>
                    cnud    = -cb2*tur_value(i,j,k)*rmre/sigma
                    cap     = ttp*cnud
                    cam     = ttm*cnud
                    anutp   = 0.5*(tur_value(i,j+1,k)+tur_value(i,j,k))
                    anutm   = 0.5*(tur_value(i,j-1,k)+tur_value(i,j,k))
                    fnup    = 0.5*(iblock%cell(i,j+1,k)%viscous/iblock%cell(i,j+1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    fnum    = 0.5*(iblock%cell(i,j-1,k)%viscous/iblock%cell(i,j-1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    cdp     = (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                    cdm     = (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                    !>
                    !>
                    bb(k,j) = -max(cdm+cam,0.0)
                    cc(k,j) =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                    dd(k,j) = -max(cdp+cap,0.0)
                enddo
                !>
                !>
                do k=2,kdim
                    !>
                    xc = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    yc = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    zc = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                    !>
                    !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    uban = xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                    sgnu = sign(1.0,uban)
                    app  = 0.5*(1.0+sgnu)
                    apm  = 0.5*(1.0-sgnu)
                    bb(k,j)=bb(k,j) - uban*app
                    cc(k,j)=cc(k,j) + uban*(app-apm)
                    dd(k,j)=dd(k,j) + uban*apm
                enddo
                !>
                !>
                do k=2,kdim
                    fact    = factor_time(i,j,k)
                    bb(k,j) = bb(k,j)*fact
                    cc(k,j) = cc(k,j)*fact+1.0*(1.0+phi)
                    dd(k,j) = dd(k,j)*fact
                    ff(k,j) = iblock%turbulent(i,j,k)%viscous*(1.0+phi)
!                    write(18070324,'(3i4,5e24.16)') i-1,j-1,k-1,ff(k,j),fact,bb(k,j),cc(k,j),dd(k,j)
                enddo
            enddo
            !>
            !> j0 boundary points
            j  = 2
            jp = min(3,jdim)
            do k=2,kdim
                volume_p = iblock%metric(i,jp,k)%volume
                !> j and j+1 average metric
                xp = iblock%metric(i,j+1,k)%gi * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                yp = iblock%metric(i,j+1,k)%gj * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                zp = iblock%metric(i,j+1,k)%gk * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                !> j-1 and i average metric
!                volume_l=iblock%metric(i,2,k)%volume
                volume_l=iblock%metric(i,1,k)%volume
                !>
                xl = iblock%metric(i,j,k)%gi * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                yl = iblock%metric(i,j,k)%gj * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                zl = iblock%metric(i,j,k)%gk * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                !> j and j+1 average metric different the (xp,yp,zp)
                xa = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                ya = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                za = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                !>
                ttp  = xp*xa+yp*ya+zp*za
                ttm  = xl*xa+yl*ya+zl*za
                !>
                cnud    =   -cb2*tur_value(i,j,k)*rmre/sigma
                cap     =   ttp*cnud
                cam     =   ttm*cnud
                !>
                anutp   =   0.5*(tur_value(i,j+1,k)+tur_value(i,j,k))
                anutm   =   0.5*(tur_value(i,j-1,k)+tur_value(i,j,k))
                !>
                fnup    =   0.5*(iblock%cell(i,j+1,k)%viscous/iblock%cell(i,j+1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                fnum    =   0.5*(iblock%cell(i,j-1,k)%viscous/iblock%cell(i,j-1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                !>
                cdp     =   (fnup+(1.0 + cb2)*anutp)*ttp*rmre/sigma
                cdm     =   (fnum+(1.0 + cb2)*anutm)*ttm*rmre/sigma
                bb(k,j) =   -max(cdm+cam,0.0)
                cc(k,j) =    max(cdp+cap,0.0) + max(cdm+cam,0.0)
                dd(k,j) =   -max(cdp+cap,0.0)
            enddo
            !>
            !>
            do k=2,kdim
                !>
                xc = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                yc = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                zc = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                !>
                uban    =   xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                sgnu    =   sign(1.0,uban)
                app     =   0.5*(1.+sgnu)
                apm     =   0.5*(1.-sgnu)
                bb(k,j) =   bb(k,j) - uban*app
                cc(k,j) =   cc(k,j) + uban*(app-apm)
                dd(k,j) =   dd(k,j) + uban*apm
            enddo
            !>
            !>
            do k=2,kdim
                fact    =   factor_time(i,j,k)
                bb(k,j) =   bb(k,j)*fact
                cc(k,j) =   cc(k,j)*fact+1.0*(1.0+phi)
                dd(k,j) =   dd(k,j)*fact
                ff(k,j) =   iblock%turbulent(i,j,k)%viscous*(1.0+phi)
!                write(18070325,'(3i4,5e24.16)') i-1,j-1,k-1,ff(k,j),fact,bb(k,j),cc(k,j),dd(k,j)
            enddo
            !>
            !> jdim boundary points
            j = jdim
            !>
            do k=2,kdim
                volume_p=iblock%metric(i,jdim+1,k)%volume
                !>
                xp = iblock%metric(i,j+1,k)%gi * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                yp = iblock%metric(i,j+1,k)%gj * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                zp = iblock%metric(i,j+1,k)%gk * iblock%metric(i,j+1,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                !>
                !> j-1 and j average metric
                volume_l=iblock%metric(i,jdim-1,k)%volume
                !>
                xl = iblock%metric(i,j,k)%gi * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                yl = iblock%metric(i,j,k)%gj * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                zl = iblock%metric(i,j,k)%gk * iblock%metric(i,j,k)%gg / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                !>
                !> j and j+1 average metric different the (xp,yp,zp)
                xa = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                ya = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                za = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg)/iblock%metric(i,j,k)%volume
                !>
                ttp=xp*xa+yp*ya+zp*za
                ttm=xl*xa+yl*ya+zl*za
                !>
                cnud    =   -cb2*tur_value(i,j,k)*rmre/sigma
                cap     =   ttp*cnud
                cam     =   ttm*cnud
                anutp   =   0.5*(tur_value(i,j+1,k)+ tur_value(i,j,k))
                anutm   =   0.5*(tur_value(i,j-1,k) + tur_value(i,j,k))
                !>
                fnup    =   0.5*(iblock%cell(i,j+1,k)%viscous/iblock%cell(i,j+1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                fnum    =   0.5*(iblock%cell(i,j-1,k)%viscous/iblock%cell(i,j-1,k)%r + iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                cdp     =   (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                cdm     =   (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                bb(k,j) =   -max(cdm+cam,0.0)
                cc(k,j) =    max(cdp+cap,0.0) + max(cdm+cam,0.0)
                dd(k,j) =   -max(cdp+cap,0.0)
            enddo
            !>
            !>
            do k=2,kdim
                xc = 0.5*(iblock%metric(i,j+1,k)%gi*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg )/iblock%metric(i,j,k)%volume
                !>
                yc = 0.5*(iblock%metric(i,j+1,k)%gj*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg )/iblock%metric(i,j,k)%volume
                !>
                zc = 0.5*(iblock%metric(i,j+1,k)%gk*iblock%metric(i,j+1,k)%gg + iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg )/iblock%metric(i,j,k)%volume
                !>
                !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                !>
                uban    =   xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                sgnu    =   sign(1.0,uban)
                app     =   0.5*(1.+sgnu)
                apm     =   0.5*(1.-sgnu)
                bb(k,j) =   bb(k,j) - uban*app
                cc(k,j) =   cc(k,j) + uban*(app-apm)
                dd(k,j) =   dd(k,j) + uban*apm
            enddo
            !>
            !>
            do k=2,kdim
                fact    =   factor_time(i,j,k)
                bb(k,j) =   bb(k,j)*fact
                cc(k,j) =   cc(k,j)*fact+1.0*(1.0+phi)
                dd(k,j) =   dd(k,j)*fact
                ff(k,j) =   iblock%turbulent(i,j,k)%viscous*(1.0+phi)
!                write(18070326,'(3i4,5e24.16)') i-1,j-1,k-1,ff(k,j),fact,bb(k,j),cc(k,j),dd(k,j)
            enddo
            !>
            !>
            if (overlap .eq. 1) then
                do k=2,kdim
                    do j=2,jdim
                        ff(k,j) =   ff(k,j)*iblock%cell(i,j,k)%blank
                        bb(k,j) =   bb(k,j)*iblock%cell(i,j,k)%blank
                        dd(k,j) =   dd(k,j)*iblock%cell(i,j,k)%blank
                        cc(k,j) =   cc(k,j)*iblock%cell(i,j,k)%blank+( 1.0 - iblock%cell(i,j,k)%blank )
                    enddo
                enddo
            end if
            !>
            !> solve a scalar tridiagonal system of equations
            !>
            do k=2,kdim
                temp(k,2) = dd(k,2)/cc(k,2)
                ff(k,2)   = ff(k,2)/cc(k,2)
            end do
            !>
            do j=3,jdim
                do k=2,kdim
                    ztmp       = 1.0/(cc(k,j)-bb(k,j)*temp(k,j-1))
                    temp(k,j) = dd(k,j)*ztmp
                    ff(k,j) = (ff(k,j)-bb(k,j)*ff(k,j-1))*ztmp
                   end do
            end do
            !>
            do ii=3,jdim
                j=jdim + 2 - ii
                do k=2,kdim
                    ff(k,j)=ff(k,j)-temp(k,j)*ff(k,j+1)
                end do
            end do
            !>
            !>
            do j=2,jdim
                do k=2,kdim
                    iblock%turbulent(i,j,k)%viscous =   ff(k,j)
!                    if(nbl ==161 .and. icyc==2) write(850,'("axis=",3i5,e24.16)') i,j,k,ff(k,j)
!                    write(1807033,'(3i4,e24.16)') i-1,j-1,k-1,ff(k,j)
                enddo
            enddo
        enddo
        !>
        !>
        !>
        !> implicit F_zeta_zeta viscous terms
        !>
        !>
!        if(nbl == 7)write(*,'("run at sa viscous calculate #9")')
!        open(unit=1807034,file='turbulent_viscous_Z_eta.dat',status='unknown')
        if(idim  .gt. 3) then
            do j=2,jdim
                !>
                !> interior points
                do i=3,idim-1
                    do k=2,kdim
                        !>
                        volume_p=iblock%metric(i+1,j,k)%volume
                        !> i and i+1 average metric
                        xp = iblock%metric(i+1,j,k)%fi * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                        !>
                        yp = iblock%metric(i+1,j,k)%fj * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                        !>
                        zp = iblock%metric(i+1,j,k)%fk * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                        !>
                        !> i-1 and i average metric
                        volume_l=iblock%metric(i-1,j,k)%volume
                        !>
                        xl = iblock%metric(i,j,k)%fi * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                        !>
                        yl = iblock%metric(i,j,k)%fj * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                        !>
                        zl = iblock%metric(i,j,k)%fk * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                        !>
                        !> i and i+1 average metric different the (xp,yp,zp)
                        xa = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        ya = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        za = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        ttp =   xp*xa+yp*ya+zp*za
                        ttm =   xl*xa+yl*ya+zl*za
                        !>
                        cnud    = -cb2*tur_value(i,j,k)*rmre/sigma
                        cap     = ttp*cnud
                        cam     = ttm*cnud
                        anutp   = 0.5*(tur_value(i+1,j,k)+tur_value(i,j,k))
                        anutm   = 0.5*(tur_value(i-1,j,k)+tur_value(i,j,k))
                        fnup    = 0.5*(iblock%cell(i+1,j,k)%viscous/iblock%cell(i+1,j,k)%r +iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                        fnum    = 0.5*(iblock%cell(i-1,j,k)%viscous/iblock%cell(i-1,j,k)%r +iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                        cdp     = (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                        cdm     = (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                        !>
                        !>
                        bb(k,i) = -max(cdm+cam,0.0)
                        cc(k,i) =  max(cdp+cap,0.0) + max(cdm+cam,0.0)
                        dd(k,i) = -max(cdp+cap,0.0)
                    enddo
                    !>
                    !>
                    do k=2,kdim
                        !>
                        xc = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        yc = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        zc = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        !>
                        !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                        uban = xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                        sgnu = sign(1.0,uban)
                        app  = 0.5*(1.0+sgnu)
                        apm  = 0.5*(1.0-sgnu)
                        bb(k,i)=bb(k,i) - uban*app
                        cc(k,i)=cc(k,i) + uban*(app-apm)
                        dd(k,i)=dd(k,i) + uban*apm
                    enddo
                    !>
                    !>
                    do k=2,kdim
                        fact    = factor_time(i,j,k)
                        bb(k,i) = bb(k,i)*fact
                        cc(k,i) = cc(k,i)*fact+1.0*(1.0+phi)
                        dd(k,i) = dd(k,i)*fact
                        ff(k,i) = iblock%turbulent(i,j,k)%viscous*(1.0+phi)
!                        write(18070327,'(3i4,6e24.16)') i-1,j-1,k-1,iblock%turbulent(i,j,k)%viscous,ff(k,i),fact,bb(k,i),cc(k,i),dd(k,i)
                    enddo
                enddo
                !>
                !> i0 boundary points
                i  = 2
                ip = min(3,idim)
                do k=2,kdim
                    volume_p = iblock%metric(ip,j,k)%volume
                    !> i and i+1 average metric
                    xp = iblock%metric(i+1,j,k)%fi * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    yp = iblock%metric(i+1,j,k)%fj * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    zp = iblock%metric(i+1,j,k)%fk * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    !> i-1 and i average metric
!                    volume_l=iblock%metric(2,j,k)%volume
                    volume_l=iblock%metric(1,j,k)%volume
                    !>
                    xl = iblock%metric(i,j,k)%fi * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    yl = iblock%metric(i,j,k)%fj * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    zl = iblock%metric(i,j,k)%fk * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    !> i and i+1 average metric different the (xp,yp,zp)
                    xa = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    ya = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    za = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    !>
                    ttp  = xp*xa+yp*ya+zp*za
                    ttm  = xl*xa+yl*ya+zl*za
                    !>
                    cnud    =   -cb2*tur_value(i,j,k)*rmre/sigma
                    cap     =   ttp*cnud
                    cam     =   ttm*cnud
                    !>
                    anutp   =   0.5*(tur_value(i+1,j,k)+tur_value(i,j,k))
                    anutm   =   0.5*(tur_value(i-1,j,k)+tur_value(i,j,k))
                    !>
                    fnup    = 0.5*(iblock%cell(i+1,j,k)%viscous/iblock%cell(i+1,j,k)%r +iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    fnum    = 0.5*(iblock%cell(i-1,j,k)%viscous/iblock%cell(i-1,j,k)%r +iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    !>
                    cdp     =   (fnup+(1.0 + cb2)*anutp)*ttp*rmre/sigma
                    cdm     =   (fnum+(1.0 + cb2)*anutm)*ttm*rmre/sigma
                    bb(k,i) =   -max(cdm+cam,0.)
                    cc(k,i) =    max(cdp+cap,0.) + max(cdm+cam,0.)
                    dd(k,i) =   -max(cdp+cap,0.)
                enddo
                !>
                !>
                do k=2,kdim
                    !>
                    xc = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    yc = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    zc = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    uban    =   xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                    sgnu    =   sign(1.0,uban)
                    app     =   0.5*(1.+sgnu)
                    apm     =   0.5*(1.-sgnu)
                    bb(k,i) =   bb(k,i) - uban*app
                    cc(k,i) =   cc(k,i) + uban*(app-apm)
                    dd(k,i) =   dd(k,i) + uban*apm
                enddo
                !>
                !>
                do k=2,kdim
                    fact    =   factor_time(i,j,k)
                    bb(k,i) =   bb(k,i)*fact
                    cc(k,i) =   cc(k,i)*fact+1.0*(1.0+phi)
                    dd(k,i) =   dd(k,i)*fact
                    ff(k,i) =   iblock%turbulent(i,j,k)%viscous*(1.0+phi)
!                    write(18070328,'(3i4,5e24.16)') i-1,j-1,k-1,ff(k,i),fact,bb(k,i),cc(k,i),dd(k,i)
                enddo
                !>
                !> idim boundary points
                i = idim
                !>
                do k=2,kdim
                    volume_p=iblock%metric(idim+1,j,k)%volume
                    !>
                    xp = iblock%metric(i+1,j,k)%fi * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    yp = iblock%metric(i+1,j,k)%fj * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    zp = iblock%metric(i+1,j,k)%fk * iblock%metric(i+1,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_p ))
                    !>
                    !> i-1 and i average metric
                    volume_l=iblock%metric(idim-1,j,k)%volume
                    !>
                    xl = iblock%metric(i,j,k)%fi * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    yl = iblock%metric(i,j,k)%fj * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    zl = iblock%metric(i,j,k)%fk * iblock%metric(i,j,k)%ff / ( 0.5*( iblock%metric(i,j,k)%volume + volume_l ))
                    !>
                    !> i and i+1 average metric different the (xp,yp,zp)
                    xa = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    ya = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    za = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    ttp=xp*xa+yp*ya+zp*za
                    ttm=xl*xa+yl*ya+zl*za
                    !>
                    cnud    =   -cb2*tur_value(i,j,k)*rmre/sigma
                    cap     =   ttp*cnud
                    cam     =   ttm*cnud
                    anutp   =   0.5*(tur_value(i+1,j,k)+tur_value(i,j,k))
                    anutm   =   0.5*(tur_value(i-1,j,k)+tur_value(i,j,k))
                    !>
                    fnup    = 0.5*(iblock%cell(i+1,j,k)%viscous/iblock%cell(i+1,j,k)%r +iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    fnum    = 0.5*(iblock%cell(i-1,j,k)%viscous/iblock%cell(i-1,j,k)%r +iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    cdp     =   (fnup+(1.0+cb2)*anutp)*ttp*rmre/sigma
                    cdm     =   (fnum+(1.0+cb2)*anutm)*ttm*rmre/sigma
                    bb(k,i) =   -max(cdm+cam,0.0)
                    cc(k,i) =    max(cdp+cap,0.0) + max(cdm+cam,0.0)
                    dd(k,i) =   -max(cdp+cap,0.0)
                enddo
                !>
                !>
                do k=2,kdim
                    xc = 0.5*(iblock%metric(i+1,j,k)%fi*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    yc = 0.5*(iblock%metric(i+1,j,k)%fj*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    zc = 0.5*(iblock%metric(i+1,j,k)%fk*iblock%metric(i+1,j,k)%ff + iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    !>tc = 0.5*(iblock%metric(i+1,j,k)%ft*iblock%metric(i,j,k)%ff + iblock%metric(i,j,k)%ft*iblock%metric(i,j,k)%ff)/iblock%metric(i,j,k)%volume
                    !>
                    uban    =   xc*iblock%cell(i,j,k)%u + yc*iblock%cell(i,j,k)%v + zc*iblock%cell(i,j,k)%w!+tc
                    sgnu    =   sign(1.0,uban)
                    app     =   0.5*(1.+sgnu)
                    apm     =   0.5*(1.-sgnu)
                    bb(k,i) =   bb(k,i) - uban*app
                    cc(k,i) =   cc(k,i) + uban*(app-apm)
                    dd(k,i) =   dd(k,i) + uban*apm
                enddo
                !>
                !>
                do k=2,kdim
                    fact    =   factor_time(i,j,k)
                    bb(k,i) =   bb(k,i)*fact
                    cc(k,i) =   cc(k,i)*fact+1.0*(1.0+phi)
                    dd(k,i) =   dd(k,i)*fact
                    ff(k,i) =   iblock%turbulent(i,j,k)%viscous*(1.0+phi)
!                    write(18070329,'(3i4,5e24.16)') i-1,j-1,k-1,ff(k,i),fact,bb(k,i),cc(k,i),dd(k,i)
                enddo
                !>
                !>
                if (overlap .eq. 1) then
                    do i=2,idim
                        do k=2,kdim
                            ff(k,i) =   ff(k,i)*iblock%cell(i,j,k)%blank
                            bb(k,i) =   bb(k,i)*iblock%cell(i,j,k)%blank
                            dd(k,i) =   dd(k,i)*iblock%cell(i,j,k)%blank
                            cc(k,i) =   cc(k,i)*iblock%cell(i,j,k)%blank+(1.-iblock%cell(i,j,k)%blank)
                        enddo
                    enddo
                end if
                !>
                !> solve a scalar tridiagonal system of equations
                do k=2,kdim
                    temp(k,2) = dd(k,2)/cc(k,2)
                    ff(k,2)   = ff(k,2)/cc(k,2)
                end do
                !>
                do i=3,idim
                    do k=2,kdim
                        ztmp       = 1.0/(cc(k,i)-bb(k,i)*temp(k,i-1))
                        temp(k,i) = dd(k,i)*ztmp
                        ff(k,i) = (ff(k,i)-bb(k,i)*ff(k,i-1))*ztmp
                    end do
                end do
                !>
                do ii=3,idim
                    i=idim + 2 - ii
                    do k=2,kdim
                        ff(k,i)=ff(k,i)-temp(k,i)*ff(k,i+1)
                    end do
                end do
                !>
                !>
                do i=2,idim
                    do k=2,kdim
                        iblock%turbulent(i,j,k)%viscous = ff(k,i)
!                        if(nbl ==1)write(860,'("axis=",3i6,e24.16)') i,j,k,ff(k,i)
!                        write(1807034,'(3i4,e24.16)') i-1,j-1,k-1,ff(k,i)
                    enddo
                enddo
            enddo
        end if
        !>
        !>
        !>
        !>update the turbulent viscous
        !>
        !>
!        if(nbl == 7)write(*,'("run at sa viscous calculate #10")')
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    if((real(tur_value(i,j,k)+iblock%turbulent(i,j,k)%viscous)) .lt. 1.0e-12) then
                        tur_value(i,j,k)                =   1.e-12
                        iblock%turbulent(i,j,k)%viscous =   0.0
                    else
                        tur_value(i,j,k)                =   tur_value(i,j,k) + iblock%turbulent(i,j,k)%viscous
                    end if
                enddo
            enddo
        enddo
        !>
        !>
        !> update viscous  and save tur_value in tur_save
        !>
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    chi                                 =   tur_value(i,j,k)/(iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r)
                    fv1                                 =   chi**3/(chi**3+cv1**3)
                    iblock%turbulent(i,j,k)%viscous     =   fv1*tur_value(i,j,k)*iblock%cell(i,j,k)%r
                    iblock%turbulent(i,j,k)%tur_save    =   tur_value(i,j,k)
!                  if(icyc==1 .and. nbl==161) write(8000+nbl,'("axis=",3i5,3e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,tur_value(i,j,k),iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r
!                  if(icyc==2 .and. nbl==161) write(7000+nbl,'("axis=",3i5,3e24.16)') i,j,k,iblock%turbulent(i,j,k)%viscous,tur_value(i,j,k),iblock%cell(i,j,k)%viscous/iblock%cell(i,j,k)%r
                enddo
            enddo
        enddo
        !>
        !>deallocate the array
        deallocate(factor_time)
        deallocate(tur_value)
        deallocate(damp)
        deallocate(temp)
        deallocate(bb)
        deallocate(cc)
        deallocate(dd)
        deallocate(ff)
        !>
        !>
        return
    end subroutine  spalart_allmaras
    !>
    !>
