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
    !*************************************************************************!
    !                                                                         !
    !  module:       baldwin_lomax.f90                                               !
    !                                                                         !
    !  programmer:   liuxz                                                    !
    !                yuanwu                                                   !
    !                                                                         !
    !  date:         2017/12/24                                               !

    !*************************************************************************!
    !>
    subroutine  baldwin_lomax(nbl)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
		implicit none
        !>
        !>
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !>
        !>
        integer :: nbl
        integer :: ibci0,&
                   ibcidim,&
                   ibcj0,&
                   ibcjdim,&
                   ibck0,&
                   ibckdim,&
                   nnn,&
                   i,j,k,&
                   ip,jp,kp,&
                   ipp,jpp,kpp,&
                   idim,jdim,kdim,&
                   iloop,iloop1,&
                   inmax,inmax1,in1,in2,infmax,&
                   inswtch
        !>
        !>
        real(kind = dprec) :: umax,umin,uvw,fblmax,ymax,fwake,fwake2,udif
        real(kind = dprec) :: ent,entmin
        real(kind = dprec) :: aplus,clauser,ccp,aplusi,xmi,ckout
        real(kind = dprec) :: c2,c,c2b,c2bp,t,yplus,term,amixl,utmax,utmin
        real(kind = dprec) :: s1,s2,suth
        real(kind = dprec),dimension(:),allocatable     :: utot
        real(kind = dprec),dimension(:),allocatable     :: damp
        real(kind = dprec),dimension(:),allocatable     :: eomu_inner
        real(kind = dprec),dimension(:),allocatable     :: eomu_outer
        real(kind = dprec),dimension(:),allocatable     :: fbl
        real(kind = dprec),dimension(:,:,:),allocatable :: emos
        !>
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        !>
        !>
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        !>
        allocate(utot(ijkmax))
        allocate(damp(ijkmax))
        allocate(eomu_inner(ijkmax))
        allocate(eomu_outer(ijkmax))
        allocate(fbl(ijkmax))
        allocate(emos(idim,jdim,kdim))
        !>
        !>
        emos   = 0.0
        aplus  = 26.e0
        ccp    = 1.6e0
        clauser= 0.0168e0
        !>
        aplusi = 1.e0/aplus
        xmi    = 1.e0/mach
        ckout  = reue*clauser*ccp/mach
        !>rmre
        suth   = 198.6/tinf
        !> determine which walls to use, based on bcj,bck,bci
        !> determine which walls to use the turbulent model
        !> the wall_blank array equit 2.0 when the bc is wall bc
        !> and no slip bc
        ibci0   = 0
        ibcidim = 0
        do k=2,kdim
            do j=2,jdim
                if (real(iblock%cell(1     ,j,k)%wall_blank ) .gt. 1.0)    ibci0=1
                if (real(iblock%cell(idim+1,j,k)%wall_blank)  .gt. 1.0)    ibcidim=1
            end do
        end do
        !>
        !>
        ibcj0   = 0
        ibcjdim = 0
        do k=2,kdim
            do i=2,idim
                if (real(iblock%cell(i,1     ,k)%wall_blank) .gt. 1.0)    ibcj0=1
                if (real(iblock%cell(i,jdim+1,k)%wall_blank) .gt. 1.0)    ibcjdim=1
            end do
        end do
        !>
        !>
        ibck0   = 0
        ibckdim = 0
        do j=2,jdim
            do i=2,idim
                if (real(iblock%cell(i,j,1     )%wall_blank) .gt. 1.0)    ibck0=1
                if (real(iblock%cell(i,j,kdim+1)%wall_blank) .gt. 1.0)    ibckdim=1
            end do
        end do
        !    write(*,'(6i5)') ibc0,ibcidim,ibcj0,ibcjdim,ibck0,ibckdim
        !>
        !>
        !>
        if (ivisc_k .gt. 1) then
            !>***************************************************************
            !>      evaluate turbulent viscosity along normal to k=0 wall
            !>***************************************************************
            !>
            !>   defaults to using min face for distance if neither min nor max contains wall
            if (ibck0 .eq. 1 .or. ibckdim .eq. 0) then

                !>
                !>loop through loop through k stations
                !>
                do  j=2,jdim
                    !>
                    !>loop through j stations
                    !>
                    do  i=2,idim
                        !>
                        !>create working arrays in i-direction
                        !>bound turbulent region
                        !>
                        iloop   = 0.80*kdim
                        inmax   = iloop
                        inmax1  = inmax-1
                        !>
                        !>bound search region for fmax
                        !>
                        in1 = inmax1*0.20+1
                        !>
                        !>
                        in2 = inmax1*0.80+1
                        !>
                        !>
                        !>chimera scheme modification....limit upper bound of search region to
                        !>the lower edge of a hole.  if hole extends to lower search bound, set
                        !>eddy viscosity to zero at all points this station.
                        if(overlap .eq. 1)then
                            !>
                            !> the overlap processing
                            !> you can rewrite teh code at here to computing the
                            !> overlap grid
                            !>
                        end if
                        !>
                        !>
                        do k=1,iloop
                            utot (k) =  sqrt( iblock%cell(i,j,k)%u**2 + iblock%cell(i,j,k)%v**2 + iblock%cell(i,j,k)%w**2)
                        end do
!!                    dist(1)       =  0.
!                    utot(1)       =  utot(2)
                        !>
                        !>
                        if(real(iblock%cell(i,j,1)%wall_blank) .eq. 2.0 )then
                            !
                            !> modify wall values in no-slip region
                            !>
                            !>
    !                        utot(1)  = sqrt( u(1,j,k)**2 + v(1,j,k)**2 + w(1,j,k)**2)
                            iblock%turbulent(i,j,1)%vorticity = (utot(2)-utot(1))/abs(iblock%turbulent(i,j,2)%distance_tur_5)
                        end if
                        !>
                        !>convenient groupings and initial values at k=1
                        !>
                        !>
                        yplus = abs(iblock%cell(i,j,1)%r*iblock%turbulent(i,j,1)%vorticity*reue*xmi/iblock%cell(i,j,1)%viscous)
                        yplus = sqrt(yplus)/26.e0
    !                    yplus   = yplusoy*disn(2)*aplus*2.0
                        if (real(iblock%cell(i,j,1)%wall_blank) .eq. 1.0 ) then
                            !>wake damping
                            term     = 1.-exp(-50.e0)
                            do kp= 2,inmax
                                damp(kp) = term
                            end do
                        else
                            !>wall damping
                            do  kp=2,inmax
                                term = yplus*iblock%turbulent(i,j,kp)%distance_tur_5*sqrt(iblock%cell(i,j,kp)%r/iblock%cell(i,j,1)%r)*(iblock%cell(i,j,1)%viscous/iblock%cell(i,j,kp)%viscous)
                                damp(kp) = 1.0e0 -exp(-term)
                            end do
                        end if
                        !>
                        !>
                        fbl(1)    = 0.0e0
                        eomu_inner(1)  = 0.0e0
                        do kp=2,inmax
                            fbl(kp)   = iblock%turbulent(i,j,kp)%vorticity*iblock%turbulent(i,j,kp)%distance_tur_5*damp(kp)
                            amixl     = 0.4*iblock%turbulent(i,j,kp)%distance_tur_5*damp(kp)
                            eomu_inner(kp) = reue/mach*iblock%cell(i,j,kp)%r*iblock%turbulent(i,j,kp)%vorticity*amixl*amixl
                        end do
                        !>
                        !>
                        !>
                        utmax = max(real(inmax),utot(1))
                        utmin = min(real(inmax),utot(1))
                        !>
                        !>
                        !>
                        if (real(iblock%cell(i,j,1)%wall_blank) .eq. 2.0 ) utmin = utot(1)
                        !>locate max value and location of fbl
                        !>
                        !>
                        fblmax = 1.e-10
                        infmax = in1
                        do kp =in1,in2
                            if (real(fbl(kp)) .gt. real(fblmax)) then
                                fblmax = fbl(kp)
                                infmax = kp
                            end if
                            !>
                            !>first maximum (degani-schiff)
                            !>
                            !>if (nvisc .eq. 3 .and. real(fbl(ip)) .lt. 0.9e0*real(fblmax))  exit
                        end do
                        !>
                        !>
                        !> curve fit of f near infmax to improve fmax location
                        !>
                        !>
                        ymax   = iblock%turbulent(i,j,infmax)%distance_tur_5
                        !>
                        !>
                        !>
                        ymax   = max(ymax,1.e-20)
                        udif   = utmax - utmin
                        fwake  = ymax*fblmax
                        fwake2 = udif*udif*ymax/fblmax
                        fwake  = min(fwake2,fwake)
                        !>
                        !>
                        !>
                        do  kp=2,inmax
                            eomu_outer(kp) = 1.e0/(1.e0+5.5e0*(0.30*iblock%turbulent(i,j,kp)%distance_tur_5/ymax)**6)
                            eomu_outer(kp) = ckout*eomu_outer(kp)*iblock%cell(i,j,kp)%r*fwake
                        end do
                        !>
                        !>
                        eomu_outer(1)  = eomu_outer(2)
                        !>
                        !> baldwin-lomax model
                        !>
                        do kp=2,in2
                            if (real(eomu_inner(kp)) .gt. real(eomu_outer(kp))) exit
                        end do
                        inswtch  = kp -1
                        do  kpp=inswtch,inmax
                            eomu_inner(kpp) = eomu_outer(kpp)
                        end do
                        !>
                        !>
                        !>
                        !>fill eomu array corresponding to k=k0 wall
                        !>
                        do k=2,iloop
                            iblock%turbulent(i,j,k)%viscous = eomu_inner(k)
                        end do
                        !>
                        !>
                        if(iloop .lt. kdim)then
                            do k=iloop+1,kdim
                            iblock%turbulent(i,j,k)%viscous = 0.0
                            end do
                        end if
                    !> j=2,jdim
                    end do
                !> k=2,kdim
                end do
            !>
            !> at i0 wall boundary
            end if
            !>
            !>
            !>
            !>
            !>***************************************************************
            !> evaluate turbulent viscosity along normal to k=kdim wall
            !>***************************************************************
            !>
            !>   defaults to using min face for distance if neither min nor max contains wall
            if (ibckdim .eq. 1) then

                !>
                !>loop through loop through k stations
                !>
                do  j=2,jdim
                    !>
                    !>loop through j stations
                    !>
                    do  i=2,idim
                        !>
                        !>create working arrays in i-direction
                        !>bound turbulent region
                        !>
                        iloop   = 0.80*kdim
                        iloop1  = iloop -1
                        inmax   = iloop
                        !>
                        iloop   = kdim - iloop1
                        iloop   = max(1,iloop)
                        inmax   = iloop+1
                        !>
                        !>bound search region for fmax
                        !>
                        in1 = iloop1*0.20+1
                        !>
                        !>
                        in2 = iloop1*0.80+1
                        !>
                        in1 = kdim - in1
                        !>
                        in2 = kdim - in2
                        !>
                        !>chimera scheme modification....limit upper bound of search region to
                        !>the lower edge of a hole.  if hole extends to lower search bound, set
                        !>eddy viscosity to zero at all points this station.
                        if(overlap .eq. 1)then
                            !>
                            !> the overlap processing
                            !> you can rewrite teh code at here to computing the
                            !> overlap grid
                            !>
                        end if
                        !>
                        !>
                        do k=kdim+1,iloop+1,-1
    !                        dist (i) =  abs( cell_dist_tur(i,j,k,2) )
                            utot (k) =  sqrt( iblock%cell(i,j,k)%u**2 + iblock%cell(i,j,k)%v**2 + iblock%cell(i,j,k)%w**2)
                        end do

                        !>
                        !>
                        if(real(iblock%cell(i,j,kdim+1)%wall_blank) .eq. 2.0 )then
                            !
                            !> modify wall values in no-slip region
                            !>
                            !>
    !                        utot(1)  = sqrt( u(1,j,k)**2 + v(1,j,k)**2 + w(1,j,k)**2)
                            iblock%turbulent(i,j,kdim+1)%vorticity = (utot(kdim)-utot(kdim+1))/abs(iblock%turbulent(i,j,kdim)%distance_tur_6)
                        end if
                        !>
                        !>convenient groupings and initial values at k=1
                        !>
                        !>
                        yplus = abs(iblock%cell(i,j,kdim+1)%r*iblock%turbulent(i,j,kdim+1)%vorticity*reue*xmi/iblock%cell(i,j,kdim+1)%viscous)
                        yplus = sqrt(yplus)/26.e0
    !                    yplus   = yplusoy*disn(2)*aplus*2.0
                        if (real(iblock%cell(i,j,kdim+1)%wall_blank) .eq. 1.0 ) then
                            !>wake damping
                            term     = 1.-exp(-50.e0)
                            do kp= kdim,inmax,-1
                                damp(kp) = term
                            end do
                        else
                            !>wall damping
                            do  kp=kdim,inmax,-1
                                term = yplus*iblock%turbulent(i,j,kp)%distance_tur_6*sqrt(iblock%cell(i,j,kp)%r/iblock%cell(i,j,kdim+1)%r)*(iblock%cell(i,j,kdim+1)%viscous/iblock%cell(i,j,kp)%viscous)
                                damp(kp) = 1.0e0 -exp(-term)
                            end do
                        end if
                        !>
                        !>
                        fbl(kdim+1)    = 0.0e0
                        eomu_inner(kdim+1)  = 0.0e0
                        do kp=kdim,inmax,-1
                            fbl(kp)   = iblock%turbulent(i,j,kp)%vorticity*iblock%turbulent(i,j,kp)%distance_tur_6*damp(kp)
                            amixl     = 0.4*iblock%turbulent(i,j,kp)%distance_tur_6*damp(kp)
                            eomu_inner(kp) = reue/mach*iblock%cell(i,j,kp)%r*iblock%turbulent(i,j,kp)%vorticity*amixl*amixl
                        end do
                        !>
                        !>
                        !>
                        utmax = max(real(iloop1),utot(inmax))
                        utmin = min(real(iloop1),utot(inmax))
                        !>
                        !>
                        !>
                        if (real(iblock%cell(i,j,kdim+1)%wall_blank) .eq. 2.0) utmin = utot(kdim+1)
                        !>locate max value and location of fbl
                        !>
                        !>
                        fblmax = 1.e-10
                        infmax = in1+2
                        do kp =in1+2,in2+2,-1
                            if (real(fbl(kp)) .gt. real(fblmax)) then
                                fblmax = fbl(kp)
                                infmax = kp
                            end if
                            !>
                            !>first maximum (degani-schiff)
                            !>
                            !>if (nvisc .eq. 3 .and. real(fbl(ip)) .lt. 0.9e0*real(fblmax))  exit
                        end do
                        !>
                        !>
                        !> curve fit of f near infmax to improve fmax location
                        !>
                        !>
                        ymax   = iblock%turbulent(i,j,infmax)%distance_tur_6
                        !>
                        !>
                        !>
                        ymax   = max(ymax,1.e-20)
                        udif   = utmax - utmin
                        fwake  = ymax*fblmax
                        fwake2 = udif*udif*ymax/fblmax
                        fwake  = min(fwake2,fwake)
                        !>
                        !>
                        !>
                        do  kp=kdim,inmax,-1
                            eomu_outer(kp) = 1.e0/(1.e0+5.5e0*(0.30*iblock%turbulent(i,j,kp)%distance_tur_6/ymax)**6)
                            eomu_outer(kp) = ckout*eomu_outer(kp)*iblock%cell(i,j,kp)%r*fwake
                        end do
                        !>
                        !>
                        eomu_outer(kdim+1)  = eomu_outer(kdim)
                        !>
                        !> baldwin-lomax model
                        !>
                        do  kp=kdim,in2+1,-1
                            if (real(eomu_inner(kp)) .gt. real(eomu_outer(kp))) exit
                        end do
                        inswtch  = kp + 1
                        do  kpp=inswtch,inmax,-1
                            eomu_inner(kpp) = eomu_outer(kpp)
                        end do
                        !>
                        !>
                        !>
                        !>fill eomu array corresponding to k=kdim wall
                        !>
                        if(ibck0 .eq. 1)then
                            do k=kdim,iloop+1,-1
                                s1  = iblock%turbulent(i,j,k)%distance_tur_5**2
                                s2  = iblock%turbulent(i,j,k)%distance_tur_6**2
                                iblock%turbulent(i,j,k)%viscous = (s1*eomu_inner(k) + s2*iblock%turbulent(i,j,k)%viscous)/(s1+s2)
                            end do
                        else
                            !>
                            !>
                            do k=kdim,iloop+1,-1
                                iblock%turbulent(i,j,k)%viscous = eomu_inner(k)
                            end do
                            if(iloop .gt. 1)then
                                do k=iloop,2,-1
                                    iblock%turbulent(i,j,k)%viscous = 0.0
                                end do
                            end if
                        end if
                    !> j=2,jdim
                    end do
                !> k=2,kdim
                end do
            !>
            !> at idim wall boundary
            end if
        !> nvisc .eq. 2
        end if
!    write(*,*) 'in the k-directions'
        !>
        !>
        if (ivisc_j .gt. 1) then
            !>***************************************************************
            !>      evaluate turbulent viscosity along normal to j=0 wall
            !>***************************************************************
            !>
            !>   defaults to using min face for distance if neither min nor max contains wall
            if (ibcj0 .eq. 1 .or. ibcjdim .eq. 0) then

                !>
                !>loop through loop through k stations
                !>
                do  i=2,idim
                    !>
                    !>loop through i stations
                    !>
                    do  k=2,kdim
                        !>
                        !>create working arrays in i-direction
                        !>bound turbulent region
                        !>
                        iloop   = 0.80*jdim
                        inmax   = iloop
                        inmax1  = inmax-1
                        !>
                        !>bound search region for fmax
                        !>
                        in1 = inmax1*0.20+1
                        !>
                        !>
                        in2 = inmax1*0.80+1
                        !>
                        !>
                        !>chimera scheme modification....limit upper bound of search region to
                        !>the lower edge of a hole.  if hole extends to lower search bound, set
                        !>eddy viscosity to zero at all points this station.
                        if(overlap .eq. 1)then
                            !>
                            !> the overlap processing
                            !> you can rewrite teh code at here to computing the
                            !> overlap grid
                            !>
                        end if
                        !>
                        !>
                        do j=1,iloop
    !                        dist (i) =  abs( cell_dist_tur(i,j,k,3) )
                            utot (j) =  sqrt( iblock%cell(i,j,k)%u**2 + iblock%cell(i,j,k)%v**2 + iblock%cell(i,j,k)%w**2)
                        end do
    !!                    dist(1)       =  0.
    !                    utot(1)       =  utot(2)
                        !>
                        !>
                        if(real(iblock%cell(i,1,k)%wall_blank) .eq. 2.0 )then
                            !
                            !> modify wall values in no-slip region
                            !>
                            !>
    !                        utot(1)  = sqrt( u(1,j,k)**2 + v(1,j,k)**2 + w(1,j,k)**2)
                            iblock%turbulent(i,1,k)%vorticity = (utot(2)-utot(1))/abs(iblock%turbulent(i,2,k)%distance_tur_3)
    !                        write(*,'("ijk=",3i4," distance=",e24.16)') j,k,i,abs(cell_dist_tur(i,2,k,3))
                        end if
                        !>
                        !>convenient groupings and initial values at k=1
                        !>
                        !>
                        yplus = abs(iblock%cell(i,1,k)%r*iblock%turbulent(i,1,k)%vorticity*reue*xmi/iblock%cell(i,1,k)%viscous)
                        yplus = sqrt(yplus)/26.e0
                        !>
    !                    write(*,'("ijk=",3i4," r=",e20.12," xmi=",e20.12," reue=",e20.12)') j,k,i,r(i,1,k),reue,xmi
                        if (real(iblock%cell(i,1,k)%wall_blank) .eq. 1.0 ) then
                            !>wake damping
                            term     = 1.-exp(-50.e0)
                            do jp= 2,inmax
                                damp(jp) = term
                            end do
                        else
                            !>wall damping
                            do  jp=2,inmax
                                term = yplus*iblock%turbulent(i,jp,k)%distance_tur_3*sqrt(iblock%cell(i,jp,k)%r/iblock%cell(i,1,k)%r)*(iblock%cell(i,1,k)%viscous/iblock%cell(i,jp,k)%viscous)
                                damp(jp) = 1.0e0 -exp(-term)
                            end do
                        end if
                        !>
                        !>
                        fbl(1)    = 0.0e0
                        eomu_inner(1)  = 0.0e0
                        do jp=2,inmax
                            fbl(jp)         = iblock%turbulent(i,jp,k)%vorticity*iblock%turbulent(i,jp,k)%distance_tur_3*damp(jp)
                            amixl           = 0.4*iblock%turbulent(i,jp,k)%distance_tur_3*damp(jp)
                            eomu_inner(jp)  = reue*xmi*iblock%cell(i,jp,k)%r*iblock%turbulent(i,jp,k)%vorticity*amixl*amixl
                        end do
                        !>
                        !>
                        !>
                        utmax = max(real(inmax),utot(1))
                        utmin = min(real(inmax),utot(1))
                        !>
                        !>
                        !>
                        if (real(iblock%cell(i,1,k)%wall_blank) .eq. 2.0 ) utmin = utot(1)
                        !>locate max value and location of fbl
                        !>
                        !>
                        fblmax = 1.e-10
                        infmax = in1
                        do jp =in1,in2
                            if (real(fbl(jp)) .gt. real(fblmax)) then
                                fblmax = fbl(jp)
                                infmax = jp
                            end if
                            !>
                            !>first maximum (degani-schiff)
                            !>
                            !>if (nvisc .eq. 3 .and. real(fbl(ip)) .lt. 0.9e0*real(fblmax))  exit
                        end do
                        !>
                        !>
                        !> curve fit of f near infmax to improve fmax location
                        !>
                        !>
                        ymax   = iblock%turbulent(i,infmax,k)%distance_tur_3
                        !>
                        !>
                        !>
                        ymax   = max(ymax,1.e-20)
                        udif   = utmax - utmin
                        fwake  = ymax*fblmax
                        fwake2 = udif*udif*ymax/fblmax
                        fwake  = min(fwake2,fwake)
                        !>
                        !>
                        !>
                        do  jp=2,inmax
                            eomu_outer(jp) = 1.e0/(1.e0+5.5e0*(0.30*iblock%turbulent(i,jp,k)%distance_tur_3/ymax)**6)
                            eomu_outer(jp) = ckout*eomu_outer(jp)*iblock%cell(i,jp,k)%r*fwake
                        end do
                        !>
                        !>
                        eomu_outer(1)  = eomu_outer(2)
                        !>
                        !> baldwin-lomax model
                        !>
                        do  jp=2,in2
                            if (real(eomu_inner(jp)) .gt. real(eomu_outer(jp))) exit
                        end do
                        inswtch  = jp -1
                        do  jpp=inswtch,inmax
                            eomu_inner(jpp) = eomu_outer(jpp)
                        end do
                        !>
                        !>
                        !>
                        !>fill eomu array corresponding to i=i0 wall
                        !>
                        do j=2,iloop
    !                        cell_viscous(i,j,k) = eomu_inner(j)
                            emos(i,j,k) = eomu_inner(j)
                        end do
                        !>
                        !>
                        if(iloop .lt. jdim)then
                            do j=iloop+1,jdim
    !                        cell_viscous(i,j,k) = 0.0
                            emos(i,j,k) = 0.0
                            end do
                        end if
                    !> j=2,jdim
                    end do
                !> k=2,kdim
                end do
    !            open(unit=1122,file='viscous-j0.dat',status='unknown')
    !            do k =2,kdim
    !                do j=2,jdim
    !                    do i=2,idim
    !                        write(1122,'("i=",i4," j=",i4," k=",i4," viscous=",e24.16)') i-1,j-1,k-1,emos(i,j,k)
    !                    end do
    !                end do
    !            end do
            !>
            !> at i0 wall boundary
            end if
            !>
            !>
            !>
            !>
            !>***************************************************************
            !> evaluate turbulent viscosity along normal to j=jdim wall
            !>***************************************************************
            !>
            !>   defaults to using min face for distance if neither min nor max contains wall
            if (ibcjdim .eq. 1) then

                !>
                !>loop through loop through k stations
                !>
                do  k=2,kdim
                    !>
                    !>loop through j stations
                    !>
                    do  i=2,idim
                        !>
                        !>create working arrays in i-direction
                        !>bound turbulent region
                        !>
                        iloop   = 0.80*jdim
                        iloop1  = iloop -1
                        inmax   = iloop
                        !>
                        iloop   = jdim - iloop1
                        iloop   = max(1,iloop)
                        inmax   = iloop+1
                        !>
                        !>bound search region for fmax
                        !>
                        in1 = iloop1*0.20+1
                        !>
                        !>
                        in2 = iloop1*0.80+1
                        !>
                        in1 = jdim - in1
                        !>
                        in2 = jdim - in2
                        !>
                        !>chimera scheme modification....limit upper bound of search region to
                        !>the lower edge of a hole.  if hole extends to lower search bound, set
                        !>eddy viscosity to zero at all points this station.
                        if(overlap .eq. 1)then
                            !>
                            !> the overlap processing
                            !> you can rewrite teh code at here to computing the
                            !> overlap grid
                            !>
                        end if
                        !>
                        !>
                        do j=jdim+1,iloop+1,-1
    !                        dist (i) =  abs( cell_dist_tur(i,j,k,2) )
                            utot (j) =  sqrt( iblock%cell(i,j,k)%u**2 + iblock%cell(i,j,k)%v**2 + iblock%cell(i,j,k)%w**2)
                        end do

                        !>
                        !>
                        if(real(iblock%cell(i,jdim+1,k)%wall_blank) .eq. 2.0 )then
                            !
                            !> modify wall values in no-slip region
                            !>
                            !>
    !                        utot(1)  = sqrt( u(1,j,k)**2 + v(1,j,k)**2 + w(1,j,k)**2)
                            iblock%turbulent(i,jdim+1,k)%vorticity = (utot(jdim)-utot(jdim+1))/abs(iblock%turbulent(i,jdim,k)%distance_tur_4)
                        end if
                        !>
                        !>convenient groupings and initial values at k=1
                        !>
                        !>
                        yplus = abs(iblock%cell(i,jdim+1,k)%r*iblock%turbulent(i,jdim+1,k)%vorticity*reue*xmi/iblock%cell(i,jdim+1,k)%viscous)
                        yplus = sqrt(yplus)/26.e0
    !                    yplus   = yplusoy*disn(2)*aplus*2.0
                        if (real(iblock%cell(i,jdim+1,k)%wall_blank) .eq. 1.0 ) then
                            !>wake damping
                            term     = 1.-exp(-50.e0)
                            do jp=  jdim,inmax,-1
                                damp(jp) = term
                            end do
                        else
                            !>wall damping
                            do  jp=jdim,inmax,-1
                                term = yplus*iblock%turbulent(i,jp,k)%distance_tur_4*sqrt(iblock%cell(i,jp,k)%r/iblock%cell(i,jdim+1,k)%r)*(iblock%cell(i,jdim+1,k)%viscous/iblock%cell(i,jp,k)%viscous)
                                damp(jp) = 1.0e0 -exp(-term)
                            end do
                        end if
                        !>
                        !>
                        fbl(jdim+1)    = 0.0e0
                        eomu_inner(jdim+1)  = 0.0e0
                        do jp=jdim,inmax,-1
                            fbl(jp)   = iblock%turbulent(i,jp,k)%vorticity*iblock%turbulent(i,jp,k)%distance_tur_4*damp(jp)
                            amixl     = 0.4*iblock%turbulent(i,jp,k)%distance_tur_4*damp(jp)
                            eomu_inner(jp) = reue*xmi*iblock%cell(i,jp,k)%r*iblock%turbulent(i,jp,k)%vorticity*amixl*amixl
                        end do
                        !>
                        !>
                        !>
                        utmax = max(real(iloop1),utot(inmax))
                        utmin = min(real(iloop1),utot(inmax))
                        !>
                        !>
                        !>
                        if (real(iblock%cell(i,jdim+1,k)%wall_blank) .eq. 2.0) utmin = utot(jdim+1)
                        !>locate max value and location of fbl
                        !>
                        !>
                        fblmax = 1.e-10
                        infmax = in1+2
                        do jp =in1+2,in2+2,-1
                            if (real(fbl(jp)) .gt. real(fblmax)) then
                                fblmax = fbl(jp)
                                infmax = jp
                            end if
                            !>
                            !>first maximum (degani-schiff)
                            !>
                            !>if (nvisc .eq. 3 .and. real(fbl(ip)) .lt. 0.9e0*real(fblmax))  exit
                        end do
                        !>
                        !>
                        !> curve fit of f near infmax to improve fmax location
                        !>
                        !>
                        ymax   = iblock%turbulent(i,infmax,k)%distance_tur_4
                        !>
                        !>
                        !>
                        ymax   = max(ymax,1.e-20)
                        udif   = utmax - utmin
                        fwake  = ymax*fblmax
                        fwake2 = udif*udif*ymax/fblmax
                        fwake  = min(fwake2,fwake)
                        !>
                        !>
                        !>
                        do  jp=jdim,inmax,-1
                            eomu_outer(jp) = 1.e0/(1.e0+5.5e0*(0.30*iblock%turbulent(i,jp,k)%distance_tur_4/ymax)**6)
                            eomu_outer(jp) = ckout*eomu_outer(jp)*iblock%cell(i,jp,k)%r*fwake
                        end do
                        !>
                        !>
                        eomu_outer(jdim+1)  = eomu_outer(jdim)
                        !>
                        !> baldwin-lomax model
                        !>
                        do  jp=jdim,in2+1,-1
                            if (real(eomu_inner(jp)) .gt. real(eomu_outer(jp))) exit
                        end do
                        inswtch  = jp + 1
                        do  jpp=inswtch,inmax,-1
                            eomu_inner(jpp) = eomu_outer(jpp)
                        end do
                        !>
                        !>
                        !>
                        !>fill eomu array corresponding to i=idim wall
                        !>
                        if(ibcj0 .eq. 1)then
                            do j=jdim,iloop+1,-1
                                s1  = iblock%turbulent(i,j,k)%distance_tur_3**2
                                s2  = iblock%turbulent(i,j,k)%distance_tur_4**2
    !                            cell_viscous(i,j,k) = (s1*eomu_inner(j) + s2*cell_viscous(i,j,k))/(s1+s2)
                                emos(i,j,k) = (s1*eomu_inner(j) + s2*emos(i,j,k))/(s1+s2)
                            end do
                        else
                            !>
                            !>
                            do j=jdim,iloop+1,-1
    !                            cell_viscous(i,j,k) = eomu_inner(j)
                                emos(i,j,k) = eomu_inner(j)
                            end do
                            if(iloop .gt. 1)then
                                do j=iloop,2,-1
    !                                cell_viscous(i,j,k) = 0.0
                                    emos(i,j,k) = 0.0
                                end do
                            end if
                        end if
                    !> j=2,jdim
                    end do
                !> k=2,kdim
                end do
            !>
            !> at idim wall boundary
            end if

            if(ivisc_k .gt. 1)then
                !>
                !> form the composite eddy-viscosity from k&j walls
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            !>
                            !>
                            if(ibck0 .eq. 0 .and. ibckdim .eq. 1)then
                                s1 = iblock%turbulent(i,j,k)%distance_tur_6**2
                            else if(ibck0 .eq.1 .and.  ibckdim .eq. 1)then
                                s1 = min(iblock%turbulent(i,j,k)%distance_tur_5,iblock%turbulent(i,j,k)%distance_tur_6)
                                s1 = s1**2
                            else
                                s1 = iblock%turbulent(i,j,k)%distance_tur_5**2
                            end if
                            !>
                            !>
                            if(ibcj0 .eq. 0 .and. ibcjdim .eq. 1)then
                                s2 = iblock%turbulent(i,j,k)%distance_tur_4**2
                            else if(ibcj0 .eq. 1 .and. ibcjdim .eq. 1)then
                                s2 = min(iblock%turbulent(i,j,k)%distance_tur_3,iblock%turbulent(i,j,k)%distance_tur_4)
                                s2 = s2**2
                            else
                                s2 = iblock%turbulent(i,j,k)%distance_tur_3**2
                            end if
                            !>
                            !>
                            iblock%turbulent(i,j,k)%viscous  = (s2*iblock%turbulent(i,j,k)%viscous + s1*emos(i,j,k))/(s1+s2)
                        end do
                    end do
                end do
                !>
                !>
                !> smooth interior values
                !>
                !>
                do nnn =1,3
                    !>
                    !>
                    do k=2,kdim
                        do j=2,jdim
                            do i=2,idim
                                emos(i,j,k) = iblock%turbulent(i,j,k)%viscous
                            end do
                        end do
                    end do
                    !>
                    !>
                    !>
                    do k=3,kdim-1
                        do j=3,jdim-1
                            do i=3,idim-1
                                iblock%turbulent(i,j,k)%viscous= 0.125e0*(emos(i-1,j-1,k-1) + emos(i-1,j+1,k-1)+&
                                                                              emos(i+1,j-1,k-1) + emos(i+1,j+1,k-1)+&
                                                                              emos(i-1,j-1,k+1) + emos(i+1,j-1,k+1)+&
                                                                              emos(i-1,j+1,k+1) + emos(i+1,j+1,k+1))
                            end do
                        end do
                    end do
                !>
                end do
    !    !> nvisc .eq. 2

            else
                !>single turbulent wall
                !> install cell_viscous array as emos array
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            iblock%turbulent(i,j,k)%viscous  = emos(i,j,k)
                        end do
                    end do
                end do
    !            open(unit=112,file='viscous.dat',status='unknown')
    !            do k =2,kdim
    !                do j=2,jdim
    !                    do i=2,idim
    !                        write(112,'("i=",i4," j=",i4," k=",i4," viscous=",e24.16)') i-1,j-1,k-1,cell_viscous(i,j,k)
    !                    end do
    !                end do
    !            end do
            !>end if the i&k face composite
            end if
            !>
            deallocate(utot)
            deallocate(damp)
            deallocate(eomu_inner)
            deallocate(eomu_outer)
            deallocate(fbl)
            deallocate(emos)
            !>
            !>
            return
        !> end if the nvisc_j .gt. 1
        end if
    !    write(*,*) 'in the j-directions'
        !>
        !>
        !>
        if (ivisc_i .gt. 1) then
            !>***************************************************************
            !>      evaluate turbulent viscosity along normal to i=0 wall
            !>***************************************************************
            !>
            !>   defaults to using min face for distance if neither min nor max contains wall
            if (ibci0 .eq. 1 .or. ibcidim .eq. 0) then

                !>
                !>loop through loop through k stations
                !>
                do  k=2,kdim
                    !>
                    !>loop through j stations
                    !>
                    do  j=2,jdim
                        !>
                        !>create working arrays in i-direction
                        !>bound turbulent region
                        !>
                        iloop   = 0.80*idim
                        inmax   = iloop
                        inmax1  = inmax-1
                        !>
                        !>bound search region for fmax
                        !>
                        in1 = inmax1*0.20+1
                        !>
                        !>
                        in2 = inmax1*0.80+1
                        !>
                        !>
                        !>chimera scheme modification....limit upper bound of search region to
                        !>the lower edge of a hole.  if hole extends to lower search bound, set
                        !>eddy viscosity to zero at all points this station.
                        if(overlap .eq. 1)then
                            !>
                            !> the overlap processing
                            !> you can rewrite teh code at here to computing the
                            !> overlap grid
                            !>
                        end if
                        !>
                        !>
                        do i=1,iloop
    !                        dist (i) =  abs( cell_dist_tur(i,j,k,1) )
                            utot (i) =  sqrt( iblock%cell(i,j,k)%u**2 + iblock%cell(i,j,k)%v**2 + iblock%cell(i,j,k)%w**2)
                        end do
    !!                    dist(1)       =  0.
    !                    utot(1)       =  utot(2)
                        !>
                        !>
                        if(real(iblock%cell(1,j,k)%wall_blank) .eq. 2.0 )then
                            !
                            !> modify wall values in no-slip region
                            !>
                            !>
    !                        utot(1)  = sqrt( u(1,j,k)**2 + v(1,j,k)**2 + w(1,j,k)**2)
                            iblock%turbulent(1,j,k)%vorticity = (utot(2)-utot(1))/abs(iblock%turbulent(2,j,k)%distance_tur_1)
                        end if
                        !>
                        !>convenient groupings and initial values at k=1
                        !>
                        !>
                        yplus = abs(iblock%cell(1,j,k)%r*iblock%turbulent(1,j,k)%vorticity*reue*xmi/iblock%cell(1,j,k)%viscous)
                        yplus = sqrt(yplus)/26.e0
    !                    yplus   = yplusoy*disn(2)*aplus*2.0
                        if (real(iblock%cell(1,j,k)%wall_blank) .eq. 1.0 ) then
                            !>wake damping
                            term     = 1.-exp(-50.e0)
                            do ip= 2,inmax
                                damp(ip) = term
                            end do
                        else
                            !>wall damping
                            do  ip=2,inmax
                                term = yplus*iblock%turbulent(ip,j,k)%distance_tur_1*sqrt(iblock%cell(ip,j,k)%r/iblock%cell(1,j,k)%r)*(iblock%cell(1,j,k)%viscous/iblock%cell(ip,j,k)%viscous)
                                damp(ip) = 1.0e0 -exp(-term)
                            end do
                        end if
                        !>
                        !>
                        fbl(1)    = 0.0e0
                        eomu_inner(1)  = 0.0e0
                        do ip=2,inmax
                            fbl(ip)   = iblock%turbulent(ip,j,k)%vorticity*iblock%turbulent(ip,j,k)%distance_tur_1*damp(ip)
                            amixl     = 0.4*iblock%turbulent(ip,j,k)%distance_tur_1*damp(ip)
                            eomu_inner(ip) = reue/mach*iblock%cell(ip,j,k)%r*iblock%turbulent(ip,j,k)%vorticity*amixl*amixl
                        end do
                        !>
                        !>
                        !>
                        utmax = max(real(inmax),utot(1))
                        utmin = min(real(inmax),utot(1))
                        !>
                        !>
                        !>
                        if (real(iblock%cell(1,j,k)%wall_blank) .eq. 2.0 ) utmin = utot(1)
                        !>locate max value and location of fbl
                        !>
                        !>
                        fblmax = 1.e-10
                        infmax = in1
                        do ip =in1,in2
                            if (real(fbl(ip)) .gt. real(fblmax)) then
                                fblmax = fbl(ip)
                                infmax = ip
                            end if
                            !>
                            !>first maximum (degani-schiff)
                            !>
                            !>if (nvisc .eq. 3 .and. real(fbl(ip)) .lt. 0.9e0*real(fblmax))  exit
                        end do
                        !>
                        !>
                        !> curve fit of f near infmax to improve fmax location
                        !>
                        !>
                        ymax   = iblock%turbulent(infmax,j,k)%distance_tur_1
                        !>
                        !>
                        !>
                        ymax   = max(ymax,1.e-20)
                        udif   = utmax - utmin
                        fwake  = ymax*fblmax
                        fwake2 = udif*udif*ymax/fblmax
                        fwake  = min(fwake2,fwake)
                        !>
                        !>
                        !>
                        do  ip=2,inmax
                            eomu_outer(ip) = 1.e0/(1.e0+5.5e0*(0.30*iblock%turbulent(ip,j,k)%distance_tur_1/ymax)**6)
                            eomu_outer(ip) = ckout*eomu_outer(ip)*iblock%cell(ip,j,k)%r*fwake
                        end do
                        !>
                        !>
                        eomu_outer(1)  = eomu_outer(2)
                        !>
                        !> baldwin-lomax model
                        !>
                        do ip=2,in2
                            if (real(eomu_inner(ip)) .gt. real(eomu_outer(ip))) exit
                        end do
                        inswtch  = ip -1
                        do  ipp=inswtch,inmax
                            eomu_inner(ipp) = eomu_outer(ipp)
                        end do
                        !>
                        !>
                        !>
                        !>fill eomu array corresponding to i=i0 wall
                        !>
                        do i=2,iloop
                            emos(i,j,k) = eomu_inner(i)
                        end do
                        !>
                        !>
                        if(iloop .lt. idim)then
                            do i=iloop+1,idim
                            emos(i,j,k) = 0.0
                            end do
                        end if
                    !> j=2,jdim
                    end do
                !> k=2,kdim
                end do
            !>
            !> at i0 wall boundary
            end if
            !>
            !>
            !>
            !>
            !>***************************************************************
            !> evaluate turbulent viscosity along normal to i=idim wall
            !>***************************************************************
            !>
            !>   defaults to using min face for distance if neither min nor max contains wall
    !        write(*,*) 'in the i-directions'
            if (ibcidim .eq. 1) then

                !>
                !>loop through loop through k stations
                !>
                do  k=2,kdim
                    !>
                    !>loop through j stations
                    !>
                    do  j=2,jdim
                        !>
                        !>create working arrays in i-direction
                        !>bound turbulent region
                        !>
                        iloop   = 0.80*idim
                        iloop1  = iloop -1
                        inmax   = iloop
                        !>
                        iloop   = idim - iloop1
                        iloop   = max(1,iloop)
                        inmax   = iloop+1
                        !>
                        !>bound search region for fmax
                        !>
                        in1 = iloop1*0.20+1
                        !>
                        !>
                        in2 = iloop1*0.80+1
                        !>
                        in1 = idim - in1
                        !>
                        in2 = idim - in2
                        !>
                        !>chimera scheme modification....limit upper bound of search region to
                        !>the lower edge of a hole.  if hole extends to lower search bound, set
                        !>eddy viscosity to zero at all points this station.
                        if(overlap .eq. 1)then
                            !>
                            !> the overlap processing
                            !> you can rewrite teh code at here to computing the
                            !> overlap grid
                            !>
                        end if
                        !>
                        !>
                        do i=idim+1,iloop+1,-1
    !                        dist (i) =  abs( cell_dist_tur(i,j,k,2) )
                            utot (i) =  sqrt( iblock%cell(i,j,k)%u**2 + iblock%cell(i,j,k)%v**2 + iblock%cell(i,j,k)%w**2)
                        end do

                        !>
                        !>
                        if(real(iblock%cell(idim+1,j,k)%wall_blank) .eq. 2.0 )then
                            !
                            !> modify wall values in no-slip region
                            !>
                            !>
    !                        utot(1)  = sqrt( u(1,j,k)**2 + v(1,j,k)**2 + w(1,j,k)**2)
                            iblock%turbulent(idim+1,j,k)%vorticity = (utot(idim)-utot(idim+1))/abs(iblock%turbulent(idim,j,k)%distance_tur_2)
                        end if
                        !>
                        !>convenient groupings and initial values at k=1
                        !>
                        !>
                        yplus = abs(iblock%cell(idim+1,j,k)%r*iblock%turbulent(idim+1,j,k)%vorticity*reue*xmi/iblock%cell(idim+1,j,k)%viscous)
                        yplus = sqrt(yplus)/26.e0
    !                    yplus   = yplusoy*disn(2)*aplus*2.0
                        if (real(iblock%cell(idim+1,j,k)%wall_blank) .eq. 1.0 ) then
                            !>wake damping
                            term     = 1.-exp(-50.e0)
                            do ip= idim,inmax,-1
                                damp(ip) = term
                            end do
                        else
                            !>wall damping
                            do  ip=idim,inmax,-1
                                term = yplus*iblock%turbulent(ip,j,k)%distance_tur_2*sqrt(iblock%cell(ip,j,k)%r /iblock%cell(idim+1,j,k)%r)*(iblock%cell(idim+1,j,k)%viscous/iblock%cell(ip,j,k)%viscous)
                                damp(ip) = 1.0e0 -exp(-term)
                            end do
                        end if
                        !>
                        !>
                        fbl(idim+1)    = 0.0e0
                        eomu_inner(idim+1)  = 0.0e0
                        do ip=idim,inmax,-1
                            fbl(ip)   = iblock%turbulent(ip,j,k)%vorticity*iblock%turbulent(ip,j,k)%distance_tur_2*damp(ip)
                            amixl     = 0.4*iblock%turbulent(ip,j,k)%distance_tur_2*damp(ip)
                            eomu_inner(ip) = reue/mach*iblock%cell(ip,j,k)%r*iblock%turbulent(ip,j,k)%vorticity*amixl*amixl
                        end do
                        !>
                        !>
                        !>
                        utmax = max(real(iloop1),utot(inmax))
                        utmin = min(real(iloop1),utot(inmax))
                        !>
                        !>
                        !>
                        if (real(iblock%cell(idim+1,j,k)%wall_blank) .eq. 2.0) utmin = utot(idim+1)
                        !>locate max value and location of fbl
                        !>
                        !>
                        fblmax = 1.e-10
                        infmax = in1+2
                        do ip =in1+2,in2+2,-1
                            if (real(fbl(ip)) .gt. real(fblmax)) then
                                fblmax = fbl(ip)
                                infmax = ip
                            end if
                            !>
                            !>first maximum (degani-schiff)
                            !>
                            !>if (nvisc .eq. 3 .and. real(fbl(ip)) .lt. 0.9e0*real(fblmax))  exit
                        end do
                        !>
                        !>
                        !> curve fit of f near infmax to improve fmax location
                        !>
                        !>
                        ymax   = iblock%turbulent(infmax,j,k)%distance_tur_2
                        !>
                        !>
                        !>
                        ymax   = max(ymax,1.e-20)
                        udif   = utmax - utmin
                        fwake  = ymax*fblmax
                        fwake2 = udif*udif*ymax/fblmax
                        fwake  = min(fwake2,fwake)
                        !>
                        !>
                        !>
                        do  ip=idim,inmax,-1
                            eomu_outer(ip) = 1.e0/(1.e0+5.5e0*(0.30*iblock%turbulent(ip,j,k)%distance_tur_2/ymax)**6)
                            eomu_outer(ip) = ckout*eomu_outer(ip)*iblock%cell(ip,j,k)%r*fwake
                        end do
                        !>
                        !>
                        eomu_outer(idim+1)  = eomu_outer(idim)
                        !>
                        !> baldwin-lomax model
                        !>
                        do  ip=idim,in2+1,-1
                            if (real(eomu_inner(ip)) .gt. real(eomu_outer(ip))) exit
                        end do
                        inswtch  = ip + 1
                        do  ipp=inswtch,inmax,-1
                            eomu_inner(ipp) = eomu_outer(ipp)
                        end do
                        !>
                        !>
                        !>
                        !>fill eomu array corresponding to i=idim wall
                        !>
                        if(ibci0 .eq. 1)then
                            do i=idim,iloop+1,-1
                                s1  = iblock%turbulent(i,j,k)%distance_tur_1**2
                                s2  = iblock%turbulent(i,j,k)%distance_tur_2**2
                                emos(i,j,k) = (s1*eomu_inner(i) + s2*emos(i,j,k))/(s1+s2)
                            end do
                        else
                            !>
                            !>
                            do i=idim,iloop+1,-1
                                emos(i,j,k) = eomu_inner(i)
                            end do
                            if(iloop .gt. 1)then
                                do i=iloop,2,-1
                                    emos(i,j,k) = 0.0
                                end do
                            end if
                        end if
                    !> j=2,jdim
                    end do
                !> k=2,kdim
                end do
            !>
            !> at idim wall boundary
            end if
        !> nvisc_i .gt. 2
        end if
        !>
        !>
        !>
        if(ivisc_k .gt. 1)then
            !>
            !>
            !>form the composite eddy-viscoucity from k&i walls
            do k=2,kdim
                do j=2,jdim
                    do i=2,idim
                        !>
                        !>
                        if(ibck0 .eq. 0 .and. ibckdim .eq. 1)then
                            s1   = iblock%turbulent(i,j,k)%distance_tur_6**2
                        else if(ibck0 .eq. 1 .and. ibckdim .eq. 1)then
                            s1   = min(iblock%turbulent(i,j,k)%distance_tur_5,iblock%turbulent(i,j,k)%distance_tur_6)
                            s1   = s1**2
                        else
                            s1   = iblock%turbulent(i,j,k)%distance_tur_5**2
                        end if
                        !>
                        !>
                        if(ibci0 .eq. 0 .and. ibcidim .eq. 1)then
                            s2   = iblock%turbulent(i,j,k)%distance_tur_2**2
                        else if(ibci0 .eq. 1 .and. ibcidim .eq. 1)then
                            s2   = min(iblock%turbulent(i,j,k)%distance_tur_1,iblock%turbulent(i,j,k)%distance_tur_2)
                            s2   = s2**2
                        else
                            s2   = iblock%turbulent(i,j,k)%distance_tur_1**2
                        end if
                        iblock%turbulent(i,j,k)%viscous  = (iblock%turbulent(i,j,k)%viscous*s2 + emos(i,j,k)*s1)/(s1 + s2)
                    end do
                end do
            end do
            !>
            !>
            !> smooth interior values
            !>
            !>
            do nnn =1,3
                !>
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            emos(i,j,k) = iblock%turbulent(i,j,k)%viscous
                        end do
                    end do
                end do
                !>
                !>
                !>
                do k=3,kdim-1
                    do j=3,jdim-1
                        do i=3,idim-1
                            iblock%turbulent(i,j,k)%viscous= 0.125e0*(emos(i-1,j-1,k-1) + emos(i-1,j+1,k-1)+&
                                                                          emos(i+1,j-1,k-1) + emos(i+1,j+1,k-1)+&
                                                                          emos(i-1,j-1,k+1) + emos(i+1,j-1,k+1)+&
                                                                          emos(i-1,j+1,k+1) + emos(i+1,j+1,k+1))
                        end do
                    end do
                end do
                !>
            end do
            !>
    !        open(unit=112,file='viscous.dat',status='unknown')
    !        do k =2,kdim
    !            do j=2,jdim
    !                do i=2,idim
    !                    write(112,'("i=",i4," j=",i4," k=",i4," viscous=",e24.16)') i-1,j-1,k-1,iblock%turbulent(i,j,k)%viscous
    !                end do
    !            end do
    !        end do

        else
            !>single turbulent wall
            !> install cell_viscous array as emos array
            do k=2,kdim
                do j=2,jdim
                    do i=2,idim
                        iblock%turbulent(i,j,k)%viscous  = emos(i,j,k)
                    end do
                end do
            end do
        !>end if the i&k face composite
        end if
        !>
        !>
        deallocate(utot)
        deallocate(damp)
        deallocate(eomu_inner)
        deallocate(eomu_outer)
        deallocate(fbl)
        deallocate(emos)
!    write(*,*) 'end the subroutine !'
        !>
        !>
        !>
        return
    end subroutine  baldwin_lomax
    !>
    !>
