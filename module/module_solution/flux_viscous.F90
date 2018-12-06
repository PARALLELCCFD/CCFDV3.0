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
!             __________________________________________________________________________________
!----------------------------------------------------------------------------------------------
!>  subroutine input_files_cgns
!>  last edit 2016-05-24
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
    subroutine flux_viscous(nbl)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
        use turbulent_module
        use metric_module
        use variable_module
        use bc_module
        use nodes_paras
        !>
        implicit none
        !>
        !>
        type(blocks_type),pointer  :: iblock
        type(bc_types),pointer     :: ibc
        type(overlap_type),pointer :: mesh
        !>
        !>
        real(kind =dprec),dimension(:,:,:),pointer :: c2,q2,viscous_ave_l,viscous_ave_t
        real(kind =dprec),dimension(:),pointer :: f2,f3,f4,f5
        real(kind =dprec) :: suth,cons,&
                             t,&
                             grp,gm1pr,&
                             gtx,gty,gtz,&
                             gtt,&
                             du,dv,dw,&
                             cc,&
                             dc,dq,wq,&
                             zta1,zta2,&
                             pran
        !>
        real*8 :: errorrn,errorrun,errorrvn,errorrwn,errorren
        integer,dimension(:),pointer :: iwall,jwall,kwall
        !>
        integer::nbl,&
                 iii,jjj,kkk,&
                 ii,jj,kk,&
                 ip,jp,kp,&
                 i,j,k,&
                 ipw,jpw,kpw,&
                 idim,jdim,kdim,&
                 ab,bb,max_ijk,&
                 minr1,maxr1,minr2,maxr2
        !>
        !> sutherland's constant
        !>
        !>
        suth  = 198.6/tinf
        grp   = gamma/prandtl
        gm1pr = (gamma-1.0)*prandtl
        cons  = rmre*1.0
        !>
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        !>
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        max_ijk = 0
        max_ijk = max(idim,jdim)
        max_ijk = max(max_ijk,kdim)+1
        !>
        !>
        allocate(f2(max_ijk))
        allocate(f3(max_ijk))
        allocate(f4(max_ijk))
        allocate(f5(max_ijk))
        !>
        allocate(c2(0:idim+2,0:jdim+2,0:kdim+2))
        allocate(q2(0:idim+2,0:jdim+2,0:kdim+2))
        allocate(viscous_ave_l(0:idim+2,0:jdim+2,0:kdim+2))
        allocate(viscous_ave_t(0:idim+2,0:jdim+2,0:kdim+2))
        !>
        !>
        allocate(iwall(num_bc))
        allocate(jwall(num_bc))
        allocate(kwall(num_bc))
        !>
        ipw = 0
        jpw = 0
        kpw = 0
        do ii =1,num_bc
            ibc => mesh%bcs(ii)
            if( ibc%block_index .ne. nbl) cycle
            if( ibc%bc_type .ne. 'wallinviscid' .and. ibc%bc_type .ne. 'wallviscid') cycle
            !>
            !>the wall bc at the i-face
            if(ibc%istart .eq. ibc%iend)then
                ipw = ipw + 1
                iwall(ipw) = ii
            end if
            !>the wall bc at the i-face
            if(ibc%jstart .eq. ibc%jend)then
                jpw = jpw + 1
                jwall(jpw) = ii
            end if
            !>the wall bc at the i-face
            if(ibc%kstart .eq. ibc%kend)then
                kpw = kpw + 1
                kwall(kpw) = ii
            end if
        end do
        !>
        !>
        do k = 0,kdim+2
            do j = 0,jdim+2
                do i = 0,idim+2
                    c2(i,j,k) = gamma*iblock%cell(i,j,k)%p / iblock%cell(i,j,k)%r
                    q2(i,j,k) = 0.5* (iblock%cell(i,j,k)%u**2 + &
                                      iblock%cell(i,j,k)%v**2 + &
                                      iblock%cell(i,j,k)%w**2 )
                end do
            end do
        end do
        !>
        !>
        !>
        !>*****************************************************
        !> viscous flux in i-direction
        !> calculate the right-hand residual contributions in
        !> the i-direction due to the viscous term
        !>*****************************************************
        !>
        !>
        !>
        do k=2,kdim
            do j=2,jdim
                do i=2,idim-1
                    !> average laminar viscosity
                    viscous_ave_l(i+1,j,k) = 0.5*(iblock%cell(i,j,k)%viscous + iblock%cell(i+1,j,k)%viscous)
                    !> average turbulent viscosity
                    if(ivisc_i .ge. 2)then
                        viscous_ave_t(i+1,j,k) =  0.5*(iblock%turbulent(i,j,k)%viscous+ iblock%turbulent(i+1,j,k)%viscous)
                    else
                        viscous_ave_t(i+1,j,k) =  0.0
                    end if
                end do
            end do
        end do
        !> average viscous coefficient
        !> i-face
        do k=2,kdim
            do j=2,jdim
                !>
                !> the cell at i0 -face
                if(iblock%cell(1,j,k)%wall_blank .gt. 1.0)then
                    t = gamma*iblock%cell(1,j,k)%p / iblock%cell(1,j,k)%r
                    viscous_ave_l(2,j,k) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_i .ge. 2)then
                        viscous_ave_t(2,j,k) = iblock%turbulent(1,j,k)%viscous
                    else
                        viscous_ave_t(2,j,k) = 0.0
                    endif
                else
                    t = 0.5*(gamma*iblock%cell(1,j,k)%p / iblock%cell(1,j,k)%r  + gamma*iblock%cell(2,j,k)%p / iblock%cell(2,j,k)%r )
                    viscous_ave_l(2,j,k) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_i .ge. 2)then
                        viscous_ave_t(2,j,k) = 0.5*(iblock%turbulent(1,j,k)%viscous + iblock%turbulent(2,j,k)%viscous)
!                        write(*,'(" i0 viscous="3i4,3e24.16)') ivisc_i,j,k,viscous_ave(2,j,k),t*sqrt(t)*(1.0+suth)/(t+suth),0.5*(iblock%turbulent(1,j,k)%viscous + iblock%turbulent(2,j,k)%viscous)
                    else
                        viscous_ave_t(2,j,k) = 0.0
!                        write(*,'(" i0 viscous="3i4,2e24.16)') ivisc_i,j,k,viscous_ave(2,j,k),t*sqrt(t)*(1.0+suth)/(t+suth)
                    endif
                endif
                !>
                !> the cell at idim - face
                if(iblock%cell(idim+1,j,k)%wall_blank .gt. 1.0)then
                    t = gamma*iblock%cell(idim+1,j,k)%p/iblock%cell(idim+1,j,k)%r
                    viscous_ave_l(idim+1,j,k) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_i .ge. 2)then
                        viscous_ave_t(idim+1,j,k) = iblock%turbulent(idim+1,j,k)%viscous
                    else
                        viscous_ave_t(idim+1,j,k) = 0.0
                    end if
                else
                    t = 0.5*(gamma*iblock%cell(idim+1,j,k)%p / iblock%cell(idim+1,j,k)%r  + gamma*iblock%cell(idim,j,k)%p / iblock%cell(idim,j,k)%r )
                    viscous_ave_l(idim+1,j,k) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_i .ge. 2)then
                        viscous_ave_t(idim+1,j,k) = 0.5*(iblock%turbulent(idim,j,k)%viscous + iblock%turbulent(idim+1,j,k)%viscous)
!                        write(*,'(" idim viscous="3i4,3e24.16)') ivisc_i,j,k,viscous_ave(idim+1,j,k),t*sqrt(t)*(1.0+suth)/(t+suth),0.5*(iblock%turbulent(idim,j,k)%viscous + iblock%turbulent(idim+1,j,k)%viscous)
                    else
                        viscous_ave_t(idim+1,j,k) = 0.0
!                        write(*,'(" idim viscous="3i4,3e24.16)') ivisc_i,j,k,viscous_ave(idim+1,j,k),t*sqrt(t)*(1.0+suth)/(t+suth)
                    end if
                endif
            end do
        end do
        !>
        !>
        !> calculate the flux contributions in i-direction
        do k=2,kdim
            do j=2,jdim

                do  i = 2,idim+1
                    !>
                    gtt  = 2.0/(iblock%metric(i-1,j,k)%volume + iblock%metric(i,j,k)%volume)
                    gtx  = iblock%metric(i,j,k)%fi*iblock%metric(i,j,k)%ff*gtt
                    gty  = iblock%metric(i,j,k)%fj*iblock%metric(i,j,k)%ff*gtt
                    gtz  = iblock%metric(i,j,k)%fk*iblock%metric(i,j,k)%ff*gtt
                    !>
                    !>
                    du  = iblock%cell(i,j,k)%u - iblock%cell(i-1,j,k)%u
                    dv  = iblock%cell(i,j,k)%v - iblock%cell(i-1,j,k)%v
                    dw  = iblock%cell(i,j,k)%w - iblock%cell(i-1,j,k)%w
                    dc  = c2(i,j,k)-c2(i-1,j,k)
                    dq  = q2(i,j,k)-q2(i-1,j,k)
                    !>
                    !>
                    wq = 0.5*((iblock%cell(i,j,k)%u + iblock%cell(i-1,j,k)%u)*gtx + &
                              (iblock%cell(i,j,k)%v + iblock%cell(i-1,j,k)%v)*gty + &
                              (iblock%cell(i,j,k)%w + iblock%cell(i-1,j,k)%w)*gtz   )

                    cc = cons*(viscous_ave_l(i,j,k)+viscous_ave_t(i,j,k))/gtt
                    !>
                    !>
                    zta1 = gtx*gtx+gty*gty+gtz*gtz
                    zta2 = (gtx*du+gty*dv+gtz*dw)/3.0
                    !>
                    !>
                    pran   = (1.0+(prandtl/prandtl_tur)*(viscous_ave_t(i,j,k)/viscous_ave_l(i,j,k)))/((viscous_ave_l(i,j,k)+viscous_ave_t(i,j,k))/viscous_ave_l(i,j,k))
                    !>

                    f2(i) = cc * (zta1*du                +gtx*zta2)
                    f3(i) = cc * (zta1*dv                +gty*zta2)
                    f4(i) = cc * (zta1*dw                +gtz*zta2)
                    f5(i) = cc * (zta1*(dq+pran*dc/(gm1*prandtl))+ (wq-0.0)*zta2)
!                    if(nbl==161 .and. icyc ==2) write(1521,'("i=",3i4,4e24.16)') i,j,k,f2(i),f3(i),f4(i),f5(i)
                    !>
                end do
                !>
                !>wall boundary
                do iii =1,ipw
                    !>
                    !>
                    ibc => mesh%bcs(iwall(iii))
                    !>
                    minr1 = min(ibc%jstart,ibc%jend)
                    maxr1 = max(ibc%jstart,ibc%jend)
                    minr2 = min(ibc%kstart,ibc%kend)
                    maxr2 = max(ibc%kstart,ibc%kend)
                    if(ibc%istart .le. 2)then
                        ii = ibc%istart
                    else
                        ii = ibc%iend + 1
                    end if
                    if( (j .ge. minr1  .and.  j .le. maxr1 ) .and. ( k .ge. minr2 .and. k .le. maxr2 ) )then
                        !>
                        !>
                        gtt  = 2.0/(iblock%metric(ii-1,j,k)%volume + iblock%metric(ii,j,k)%volume)
!                        if(ii .eq. 2)write(0925,'("i0 face volume=",i5,2e24.16)') nbl,iblock%metric(ii-1,j,k)%volume,iblock%metric(ii,j,k)%volume
!                        if(ii .gt. 2)write(0925,'("idim face volume=",i5,2e24.16)') nbl,iblock%metric(ii-1,j,k)%volume,iblock%metric(ii,j,k)%volume
                        gtx  =iblock%metric(ii,j,k)%fi*iblock%metric(ii,j,k)%ff*gtt
                        gty  =iblock%metric(ii,j,k)%fj*iblock%metric(ii,j,k)%ff*gtt
                        gtz  =iblock%metric(ii,j,k)%fk*iblock%metric(ii,j,k)%ff*gtt
                        !>
                        !>
                        du  = 2.0*( iblock%cell(ii,j,k)%u - iblock%cell(ii-1,j,k)%u)
                        dv  = 2.0*( iblock%cell(ii,j,k)%v - iblock%cell(ii-1,j,k)%v)
                        dw  = 2.0*( iblock%cell(ii,j,k)%w - iblock%cell(ii-1,j,k)%w)
                        dq  = 2.0*( q2(ii,j,k)-q2(ii-1,j,k))
                        dc  = 2.0*( c2(ii,j,k)-c2(ii-1,j,k))
                        !>
                        !>
                        if(ibc%istart .le. 2)then
                            wq = iblock%cell(ii-1,j,k)%u*gtx + iblock%cell(ii-1,j,k)%v*gty + iblock%cell(ii-1,j,k)%w*gtz
                        else
                            wq = iblock%cell(ii  ,j,k)%u*gtx + iblock%cell(ii  ,j,k)%v*gty + iblock%cell(ii  ,j,k)%w*gtz
                        end if

                        cc   =cons*(viscous_ave_l(ii,j,k)+viscous_ave_t(ii,j,k))/gtt
                        !>
                        !>
                        zta1 = gtx*gtx+gty*gty+gtz*gtz
                        zta2 = (gtx*du+gty*dv+gtz*dw)/3.0
                        !>
                        !>
                        !>
                        pran   = (1.0+(prandtl/prandtl_tur)*(viscous_ave_t(ii,j,k)/viscous_ave_l(ii,j,k)))/((viscous_ave_l(ii,j,k)+viscous_ave_t(ii,j,k))/viscous_ave_l(ii,j,k))
                        !>
                        f2(ii) = cc * (zta1*du                +gtx*zta2)
                        f3(ii) = cc * (zta1*dv                +gty*zta2)
                        f4(ii) = cc * (zta1*dw                +gtz*zta2)
                        f5(ii) = cc * (zta1*(dq+pran*dc/(gm1*prandtl))+ (wq-0.0)*zta2)
                    end if

                end do
                !>
                !> viscous flux = fv(i+1/2) - fv(i-1/2)
                do i = 2,idim
!                    if(nbl==161 .and. icyc ==2) write(180,'(3i4,4e24.16)') i,j,k,f2(i)-f2(i+1),f3(i)-f3(i+1),f4(i)-f4(i+1),f5(i)-f5(i+1)
                    iblock%variable(i,j,k)%res_2 = iblock%variable(i,j,k)%res_2 - (f2(i+1)-f2(i))
                    iblock%variable(i,j,k)%res_3 = iblock%variable(i,j,k)%res_3 - (f3(i+1)-f3(i))
                    iblock%variable(i,j,k)%res_4 = iblock%variable(i,j,k)%res_4 - (f4(i+1)-f4(i))
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_5 - (f5(i+1)-f5(i))
                end do
                !>
                !>
            end do
        end do
        !>
        !>
        !>
        !>*****************************************************
        !> viscous flux in j-direction
        !> calculate the right-hand residual contributions in
        !> the j-direction due to the viscous term
        !>*****************************************************
        !> average viscous coefficient
        do k=2,kdim
            do i=2,idim
                do j=2,jdim-1
                     viscous_ave_l(i,j+1,k) =  0.5*(iblock%cell(i,j,k)%viscous + iblock%cell(i,j+1,k)%viscous)
                    if(ivisc_j .ge. 2)then
                        viscous_ave_t(i,j+1,k) =  0.5*(iblock%turbulent(i,j,k)%viscous+ iblock%turbulent(i,j+1,k)%viscous)
                    else
                        viscous_ave_t(i,j+1,k) =  0.0
                    end if
                end do
            end do
        end do
        !>
        !> j -face
        do k=2,kdim
            do i=2,idim
                !>
                !> the cell at j0 -face
                if(iblock%cell(i,1,k)%wall_blank .gt. 1.0)then
                    t = gamma*iblock%cell(i,1,k)%p/iblock%cell(i,1,k)%r
                    viscous_ave_l(i,2,k) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_j .ge. 2)then
                        viscous_ave_t(i,2,k) = iblock%turbulent(i,1,k)%viscous
                    else
                        viscous_ave_t(i,2,k) = 0.0
                    endif
                else
                    t = 0.5*(gamma*iblock%cell(i,1,k)%p/iblock%cell(i,1,k)%r  + gamma*iblock%cell(i,2,k)%p/iblock%cell(i,2,k)%r)
                    viscous_ave_l(i,2,k) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_j .ge. 2)then
                        viscous_ave_t(i,2,k) = 0.5*(iblock%turbulent(i,1,k)%viscous + iblock%turbulent(i,2,k)%viscous)
                    else
                        viscous_ave_t(i,2,k) = 0.0
                    endif
                endif
                !>
                !> the cell at jdim - face
                if(iblock%cell(i,jdim+1,k)%wall_blank .gt. 1.0)then
                    t = gamma*iblock%cell(i,jdim+1,k)%p/iblock%cell(i,jdim+1,k)%r
                    viscous_ave_l(i,jdim+1,k) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_j .ge. 2)then
                        viscous_ave_t(i,jdim+1,k) = iblock%turbulent(i,jdim+1,k)%viscous
                    else
                        viscous_ave_t(i,jdim+1,k) = 0.0
                    end if
                else
                    t = 0.5*(gamma*iblock%cell(i,jdim+1,k)%p/iblock%cell(i,jdim+1,k)%r  + gamma*iblock%cell(i,jdim,k)%p/iblock%cell(i,jdim,k)%r)
                    viscous_ave_l(i,jdim+1,k) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_j .ge. 2)then
                        viscous_ave_t(i,jdim+1,k) = 0.5*(iblock%turbulent(i,jdim,k)%viscous + iblock%turbulent(i,jdim+1,k)%viscous)
                    else
                        viscous_ave_t(i,jdim+1,k) = 0.0
                    end if
                endif
                !>
            end do
        end do
        !>
        !>
        !> calculate the flux contributions in j-direction
        !>
        do i=2,idim
            do k=2,kdim
                do j = 2,jdim+1
                    !>
                    !>
                    gtt  = 2.0/(iblock%metric(i,j-1,k)%volume + iblock%metric(i,j,k)%volume)
                    gtx  = iblock%metric(i,j,k)%gi*iblock%metric(i,j,k)%gg*gtt
                    gty  = iblock%metric(i,j,k)%gj*iblock%metric(i,j,k)%gg*gtt
                    gtz  = iblock%metric(i,j,k)%gk*iblock%metric(i,j,k)%gg*gtt
                    !>
                    !>
                    du  = iblock%cell(i,j,k)%u - iblock%cell(i,j-1,k)%u
                    dv  = iblock%cell(i,j,k)%v - iblock%cell(i,j-1,k)%v
                    dw  = iblock%cell(i,j,k)%w - iblock%cell(i,j-1,k)%w
                    dc  = c2(i,j,k)-c2(i,j-1,k)
                    dq  = q2(i,j,k)-q2(i,j-1,k)
                    !>
                    !>
                    wq  =0.5*((iblock%cell(i,j,k)%u + iblock%cell(i,j-1,k)%u)*gtx  + &
                              (iblock%cell(i,j,k)%v + iblock%cell(i,j-1,k)%v)*gty  + &
                              (iblock%cell(i,j,k)%w + iblock%cell(i,j-1,k)%w)*gtz    )

                    cc  =cons*(viscous_ave_l(i,j,k)+viscous_ave_t(i,j,k))/gtt
!                    if(nbl ==161 ) write(6000+nbl,'("cc=",3i4,3e24.16)') i,j,k,cc,gtt,viscous_ave_l(i,j,k)+viscous_ave_t(i,j,k)
                    !>
                    !>
                    !>
                    !>
                    zta1 = gtx*gtx+gty*gty+gtz*gtz
                    zta2 = (gtx*du+gty*dv+gtz*dw)/3.0
                    !>
                    !>
                    pran   = (1.0+(prandtl/prandtl_tur)*(viscous_ave_t(i,j,k)/viscous_ave_l(i,j,k)))/((viscous_ave_l(i,j,k)+viscous_ave_t(i,j,k))/viscous_ave_l(i,j,k))
                    !>
                    !>
                    f2(j) = cc * (zta1*du                +gtx*zta2)
                    f3(j) = cc * (zta1*dv                +gty*zta2)
                    f4(j) = cc * (zta1*dw                +gtz*zta2)
                    f5(j) = cc * (zta1*(dq+pran*dc/(gm1*prandtl))+ (wq-0.0)*zta2)
!                    write(1521,'("f=",3i4,5e24.16)') i,j,k,cc,f2(j),f3(j),f4(j),f5(j)
!                    write(1521,'("p=",3i4,8e24.16)') i,j,k,gtx,gty,gtz,wq,zta2,dc,dq,zta1
!                    write(1521,'("c=",3i4,4e24.16)') i,j,k,cc,cons,viscous_ave(i,j,k),gtt
                end do
                !>
                !>wall boundary
                do jjj =1,jpw
                    !>
                    !>
                    ibc => mesh%bcs(jwall(jjj))
                    !>
                    minr1 = min(ibc%istart,ibc%iend)
                    maxr1 = max(ibc%istart,ibc%iend)
                    minr2 = min(ibc%kstart,ibc%kend)
                    maxr2 = max(ibc%kstart,ibc%kend)
                    if(ibc%jstart .le. 2)then
                        jj = ibc%jstart
                    else
                        jj = ibc%jend + 1
                    end if
                    if( (i .ge. minr1  .and.  i .le. maxr1 ) .and. ( k .ge. minr2 .and. k .le. maxr2 ) )then
                        !>
                        !>
                        gtt  = 2.0/(iblock%metric(i,jj-1,k)%volume + iblock%metric(i,jj,k)%volume)
                        gtx  =iblock%metric(i,jj,k)%gi*iblock%metric(i,jj,k)%gg*gtt
                        gty  =iblock%metric(i,jj,k)%gj*iblock%metric(i,jj,k)%gg*gtt
                        gtz  =iblock%metric(i,jj,k)%gk*iblock%metric(i,jj,k)%gg*gtt
                        !>
                        !>
                        du  = 2.0*( iblock%cell(i,jj,k)%u - iblock%cell(i,jj-1,k)%u)
                        dv  = 2.0*( iblock%cell(i,jj,k)%v - iblock%cell(i,jj-1,k)%v)
                        dw  = 2.0*( iblock%cell(i,jj,k)%w - iblock%cell(i,jj-1,k)%w)
                        dq  = 2.0*( q2(i,jj,k) - q2(i,jj-1,k) )
                        dc  = 2.0*( c2(i,jj,k) - c2(i,jj-1,k) )
                        !>
                        !>
                        if(ibc%jstart .le. 2)then
                            wq = iblock%cell(i,jj-1,k)%u*gtx + iblock%cell(i,jj-1,k)%v*gty + iblock%cell(i,jj-1,k)%w*gtz
                        else
                            wq = iblock%cell(i,jj  ,k)%u*gtx + iblock%cell(i,jj  ,k)%v*gty + iblock%cell(i,jj  ,k)%w*gtz
                        end if
                        !>
                        !>
                        cc   =cons*(viscous_ave_l(i,jj,k)+viscous_ave_t(i,jj,k))/gtt
                        !>
                        zta1 = gtx*gtx+gty*gty+gtz*gtz
                        zta2 = (gtx*du+gty*dv+gtz*dw)/3.0
!                        if(nbl ==161 .and. icyc ==2 ) write(6000+nbl,'("cc=",3i4,5e24.16)') i,j,k,cc,gtt,viscous_ave_l(i,jj,k)+viscous_ave_t(i,jj,k),zta1,zta2
                        !>
                        !>
                        !>
                        pran   = (1.0+(prandtl/prandtl_tur)*(viscous_ave_t(i,jj,k)/viscous_ave_l(i,jj,k)))/((viscous_ave_l(i,jj,k)+viscous_ave_t(i,jj,k))/viscous_ave_l(i,jj,k))
                        !>
                        !>
                        f2(jj) = cc * (zta1*du                +gtx*zta2)
                        f3(jj) = cc * (zta1*dv                +gty*zta2)
                        f4(jj) = cc * (zta1*dw                +gtz*zta2)
                        f5(jj) = cc * (zta1*(dq+pran*dc/(gm1*prandtl))+ (wq-0.0)*zta2)
!                        if(nbl==161 .and. icyc ==2) write(1522,'("j=",3i4,4e24.16)') i,jj,k,f2(jj),f3(jj),f4(jj),f5(jj)
!                        write(1522,'("c=",3i4,4e24.16)') i,jj,k,cc,cons,viscous_ave(i,jj,k),gtt
!                        write(1522,'("p=",3i4,9e24.16)') i,j,k,gtx,gty,gtz,wq,zta2,dc,dq,zta1,cc
                    end if

                end do
                !>
                !> viscous flux = gv(j+1/2) - gv(j-1/2)
                do j = 2,jdim
                    !>
                    iblock%variable(i,j,k)%res_2 =iblock%variable(i,j,k)%res_2 - (f2(j+1) - f2(j))
                    iblock%variable(i,j,k)%res_3 =iblock%variable(i,j,k)%res_3 - (f3(j+1) - f3(j))
                    iblock%variable(i,j,k)%res_4 =iblock%variable(i,j,k)%res_4 - (f4(j+1) - f4(j))
                    iblock%variable(i,j,k)%res_5 =iblock%variable(i,j,k)%res_5 - (f5(j+1) - f5(j))
!                    if(nbl==161 .and. icyc ==2) write(190,'(3i4,4e24.16)') i,j,k,f2(j)-f2(j+1),f3(j)-f3(j+1),f4(j)-f4(j+1),f5(j)-f5(j+1)
                end do
            end do
        end do
        !>
        !>
        !>*****************************************************
        !> viscous flux in k-direction
        !> calculate the right-hand residual contributions in
        !> the k-direction due to the viscous term
        !>*****************************************************
        !>
        !> average viscous coefficient
        do j=2,jdim
            do i=2,idim
                do k=2,kdim-1
                    viscous_ave_l(i,j,k+1) =  0.5*(iblock%cell(i,j,k)%viscous + iblock%cell(i,j,k+1)%viscous)
                    if(ivisc_k .ge. 2)then
                        viscous_ave_t(i,j,k+1) =  0.5*(iblock%turbulent(i,j,k)%viscous+ iblock%turbulent(i,j,k+1)%viscous)
                    else
                        viscous_ave_t(i,j,k+1) =  0.0
                    end if
                end do
            end do
        end do
        !> k -face
        do j=2,jdim
            do i=2,idim
                !>
                !> the cell at k0 -face
                if(iblock%cell(i,j,1)%wall_blank .gt. 1.0)then
                    t = gamma*iblock%cell(i,j,1)%p/iblock%cell(i,j,1)%r
                    viscous_ave_l(i,j,2) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_k .ge. 2)then
                        viscous_ave_t(i,j,2) = iblock%turbulent(i,j,1)%viscous
                    else
                        viscous_ave_t(i,j,2) = 0.0
                    endif
                else
                    t = 0.5*(gamma*iblock%cell(i,j,1)%p/iblock%cell(i,j,1)%r  + gamma*iblock%cell(i,j,2)%p/iblock%cell(i,j,2)%r)
                    viscous_ave_l(i,j,2) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_k .ge. 2)then
                        viscous_ave_t(i,j,2) = 0.5*(iblock%turbulent(i,j,1)%viscous + iblock%turbulent(i,j,2)%viscous)
                    else
                        viscous_ave_t(i,j,2) = 0.0
                    endif
                endif
                !>
                !> the cell at kdim - face
                if(iblock%cell(i,j,kdim+1)%wall_blank .gt. 1.0)then
                    t = gamma*iblock%cell(i,j,kdim+1)%p/iblock%cell(i,j,kdim+1)%r
                    viscous_ave_l(i,j,kdim+1) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_k .ge. 2)then
                        viscous_ave_t(i,j,kdim+1) = iblock%turbulent(i,j,kdim+1)%viscous
                    else
                        viscous_ave_t(i,j,kdim+1) = 0.0
                    end if
                else
                    t = 0.5*(gamma*iblock%cell(i,j,kdim+1)%p/iblock%cell(i,j,kdim+1)%r  + gamma*iblock%cell(i,j,kdim)%p/iblock%cell(i,j,kdim)%r)
                    viscous_ave_l(i,j,kdim+1) = t*sqrt(t)*(1.0+suth)/(t+suth)
                    if(ivisc_k .ge. 2)then
                        viscous_ave_t(i,j,kdim+1) = 0.5*(iblock%turbulent(i,j,kdim)%viscous + iblock%turbulent(i,j,kdim+1)%viscous)
                    else
                        viscous_ave_t(i,j,kdim+1) = 0.0
                    end if
                endif
            end do
        end do
        !>
        !>
        !> calculate the flux contributions in k-direction
        do j=2,jdim
            do i=2,idim

                do k = 2,kdim+1
                    !>
                    !>
                    gtt  = 2.0/(iblock%metric(i,j,k-1)%volume + iblock%metric(i,j,k)%volume)
                    gtx  =iblock%metric(i,j,k)%hi*iblock%metric(i,j,k)%hh*gtt
                    gty  =iblock%metric(i,j,k)%hj*iblock%metric(i,j,k)%hh*gtt
                    gtz  =iblock%metric(i,j,k)%hk*iblock%metric(i,j,k)%hh*gtt
                    !>
                    !>
                    !>
                    du  =  iblock%cell(i,j,k)%u- iblock%cell(i,j,k-1)%u
                    dv  =  iblock%cell(i,j,k)%v- iblock%cell(i,j,k-1)%v
                    dw  =  iblock%cell(i,j,k)%w- iblock%cell(i,j,k-1)%w
                    !>
                    !>
                    dc  = 1.0*(c2(i,j,k)-c2(i,j,k-1))
                    dq  = 1.0*(q2(i,j,k)-q2(i,j,k-1))
                    wq  = 0.5*1.0*((iblock%cell(i,j,k)%u + iblock%cell(i,j,k-1)%u)*gtx + &
                                   (iblock%cell(i,j,k)%v + iblock%cell(i,j,k-1)%v)*gty + &
                                   (iblock%cell(i,j,k)%w + iblock%cell(i,j,k-1)%w)*gtz   )

                    cc  =cons*(viscous_ave_l(i,j,k)+viscous_ave_t(i,j,k))/gtt
                    !>
                    !>
                    zta1 = gtx*gtx+gty*gty+gtz*gtz
                    zta2 = (gtx*du+gty*dv+gtz*dw)/3.0
                    !>
                    !>
                    pran   = (1.0+(prandtl/prandtl_tur)*(viscous_ave_t(i,j,k)/viscous_ave_l(i,j,k)))/((viscous_ave_l(i,j,k)+viscous_ave_t(i,j,k))/viscous_ave_l(i,j,k))
                    !>
                    !>
                    f2(k) = cc * (zta1*du                +gtx*zta2)
                    f3(k) = cc * (zta1*dv                +gty*zta2)
                    f4(k) = cc * (zta1*dw                +gtz*zta2)
                    f5(k) = cc * (zta1*(dq+pran*dc/(gm1*prandtl))+ (wq-0.0)*zta2)

                end do
                !>
                !>wall boundary
                do kkk =1,kpw
                    !>
                    !>
                    ibc => mesh%bcs(kwall(kkk))
                    !>
                    minr1 = min(ibc%istart,ibc%iend)
                    maxr1 = max(ibc%istart,ibc%iend)
                    minr2 = min(ibc%jstart,ibc%jend)
                    maxr2 = max(ibc%jstart,ibc%jend)
                    if(ibc%kstart .le. 2)then
                        kk = ibc%kstart
                    else
                        kk = ibc%kend + 1
                    end if
                    if( (i .ge. minr1  .and.  i .le. maxr1 ) .and. ( j .ge. minr2 .and. j .le. maxr2 ) )then
                        !>
                        !>
                        gtt  = 2.0/(iblock%metric(i,j,kk-1)%volume + iblock%metric(i,j,kk)%volume)
                        gtx  =iblock%metric(i,j,kk)%hi*iblock%metric(i,j,kk)%hh*gtt
                        gty  =iblock%metric(i,j,kk)%hj*iblock%metric(i,j,kk)%hh*gtt
                        gtz  =iblock%metric(i,j,kk)%hk*iblock%metric(i,j,kk)%hh*gtt
                        !>
                        du  = 2.0*( iblock%cell(i,j,kk)%u - iblock%cell(i,j,kk-1)%u)
                        dv  = 2.0*( iblock%cell(i,j,kk)%v - iblock%cell(i,j,kk-1)%v)
                        dw  = 2.0*( iblock%cell(i,j,kk)%w - iblock%cell(i,j,kk-1)%w)
                        !>
                        !>
                        dq  = 2.0*(q2(i,j,kk)-q2(i,j,kk-1))
                        dc  = 2.0*(c2(i,j,kk)-c2(i,j,kk-1))

                        if(ibc%kstart .le. 2)then
                            wq = iblock%cell(i,j,kk-1)%u*gtx + iblock%cell(i,j,kk-1)%v*gty + iblock%cell(i,j,kk-1)%w*gtz
                        else
                            wq = iblock%cell(i,j,kk  )%u*gtx + iblock%cell(i,j,kk  )%v*gty + iblock%cell(i,j,kk  )%w*gtz
                        end if


                        cc   =cons*(viscous_ave_l(i,j,kk)+viscous_ave_t(i,j,kk))/gtt
                        !>
                        !>
                        zta1 = gtx*gtx+gty*gty+gtz*gtz
                        zta2 = (gtx*du+gty*dv+gtz*dw)/3.0
                        !>
                        !>
                        pran   = (1.0+(prandtl/prandtl_tur)*(viscous_ave_t(i,j,kk)/viscous_ave_l(i,j,kk)))/((viscous_ave_l(i,j,kk)+viscous_ave_t(i,j,kk))/viscous_ave_l(i,j,kk))
                        !>
                        !>
                        f2(kk) = cc * (zta1*du                +gtx*zta2)
                        f3(kk) = cc * (zta1*dv                +gty*zta2)
                        f4(kk) = cc * (zta1*dw                +gtz*zta2)
                        f5(kk) = cc * (zta1*(dq+pran*dc/(gm1*prandtl))+ (wq-0.0)*zta2)
                    end if

                end do
                !>
                !> viscous flux = hv(k+1/2) - hv(k-1/2)
                do k = 2,kdim
!                    if(nbl==161 )write(17011610+icyc,'(3i4,5e24.16)') i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
                    iblock%variable(i,j,k)%res_2 =iblock%variable(i,j,k)%res_2 - (f2(k+1) - f2(k))
                    iblock%variable(i,j,k)%res_3 =iblock%variable(i,j,k)%res_3 - (f3(k+1) - f3(k))
                    iblock%variable(i,j,k)%res_4 =iblock%variable(i,j,k)%res_4 - (f4(k+1) - f4(k))
                    iblock%variable(i,j,k)%res_5 =iblock%variable(i,j,k)%res_5 - (f5(k+1) - f5(k))
!                    if(nbl==161 .and. icyc ==2)write(170,'(3i4,4e24.16)') i,j,k,f2(k)-f2(k+1),f3(k)-f3(k+1),f4(k)-f4(k+1),f5(k)-f5(k+1)
!                    if(nbl==95)write(17053,'(3i4,4e24.16)') i,j,k,-f2(k)+f2(k+1),-f3(k)+f3(k+1),-f4(k)+f4(k+1),-f5(k)+f5(k+1)
!                    if(nbl==95)write(170531,'(3i4,4e24.16)') i,j,k,iblock%variable(i,j,k)%res_2,iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
!                    if(nbl==151 )write(1701510+icyc,'(3i4,5e24.16)') i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
!                    if(nbl==161 )write(1701610+icyc,'(3i4,5e24.16)') i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
!                    if(nbl==147 )write(1701470+icyc,'(3i4,5e24.16)') i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
!                    if(nbl==108 )write(1701080+icyc,'(3i4,5e24.16)') i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
!                    if(nbl==154 )write(1701540+icyc,'(3i4,5e24.16)') i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
!                    if(nbl==157 )write(1701570+icyc,'(3i4,5e24.16)') i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,iblock%variable(i,j,k)%res_3,iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
                end do
            end do
        end do
        !>
        !>
        !>
        !>
        !>
        deallocate(c2)
        deallocate(q2)
        deallocate(viscous_ave_l)
        deallocate(viscous_ave_t)
        deallocate(f2)
        deallocate(f3)
        deallocate(f4)
        deallocate(f5)
        deallocate(iwall)
        deallocate(jwall)
        deallocate(kwall)



        return
    end subroutine flux_viscous



