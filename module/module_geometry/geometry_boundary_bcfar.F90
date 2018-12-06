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
    subroutine bc_far(num_of_bc,nbl)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use bc_module
        use blocks_module
		!>
        !>
		implicit none
        !>
        type(bc_types),pointer     :: ibc
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !>
        integer :: num_of_bc,&
                   i,j,k,&
                   istep,jstep,kstep,&
                   nbl
        integer :: is,ie,&
                   js,je,&
                   ks,ke,&
                   id1,jd1,kd1,&
                   id2,jd2,kd2
        character(len =20) :: iface
        real(kind = dprec) :: btt,ctt,cx,cy,cz,ax,ay,az,chi
        real(kind = dprec) :: ut,vt,wt,ci,ce,qi,qe,rnt,rp,qn,c,t9,t12,t13,t16,t17,t18,t19,t20,ent0
        real(kind = dprec) :: logical_depending,uban1,uban0,uban,aface,uface,c1,c0,riemann1,riemann0,epson
        !>
        !>
        mesh   => grids(imesh)
        ibc    => mesh%bcs(num_of_bc)
        iblock => mesh%blocks(nbl)
        !>
        !>
        !>
        epson  = 1.e-15
        chi    = 1.341946
        !>
        !>
        is = ibc%istart
        ie = ibc%iend
        js = ibc%jstart
        je = ibc%jend
        ks = ibc%kstart
        ke = ibc%kend
        !>
        !>
        !> boundary conditions for opposed symmetry
        if (is .eq. ie) then
            iface = 'i'
        else if(js .eq. je)then
            iface = 'j'
        else if(ks .eq. ke)then
            iface = 'k'
        else
            write(*,*) 'ccfdv3.0:: have some error at the boundary face!!!'
        end if
        !>
        !>
        !>
        istep = sign(1,ie - is)
        jstep = sign(1,je - js)
        kstep = sign(1,ke - ks)
        !>
        !>
        !>
        !> fill in the variable for  the boundary cell
        !>
        do i = is,ie,istep
            do k = ks,ke,kstep
                do j = js,je,jstep
                    !>
                    !>
                    id1 = i
                    jd1 = j
                    kd1 = k
                    !>
                    !>
                    id2 = i
                    jd2 = j
                    kd2 = k
                    !>
                    !>
                    if (iface .eq. 'i' ) then
                        !>
                        !>
                        id1 = i+1*ibc%direction
                        id2 = i+2*ibc%direction
                        !>
                        !>
                        !>
                        ax = ibc%norm_index*iblock%metric(id1,jd1,kd1)%fi
                        ay = ibc%norm_index*iblock%metric(id1,jd1,kd1)%fj
                        az = ibc%norm_index*iblock%metric(id1,jd1,kd1)%fk
                    else if ( iface  .eq. 'j') then
                        !>
                        !>
                        jd1 = j+1*ibc%direction
                        jd2 = j+2*ibc%direction
                        !>
                        !>
                        !>
                        ax = ibc%norm_index*iblock%metric(id1,jd1,kd1)%gi
                        ay = ibc%norm_index*iblock%metric(id1,jd1,kd1)%gj
                        az = ibc%norm_index*iblock%metric(id1,jd1,kd1)%gk

                    else if ( iface  .eq. 'k') then
                        !>
                        !>
                        kd1 = k+1*ibc%direction
                        kd2 = k+2*ibc%direction
                        !>
                        !>
                        !>
                        ax = ibc%norm_index*iblock%metric(id1,jd1,kd1)%hi
                        ay = ibc%norm_index*iblock%metric(id1,jd1,kd1)%hj
                        az = ibc%norm_index*iblock%metric(id1,jd1,kd1)%hk

                    end if

                    ut  =    iblock%cell(i,j,k)%u
                    vt  =    iblock%cell(i,j,k)%v
                    wt  =    iblock%cell(i,j,k)%w

                    c1  = sqrt(gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r)
                    c0  = sqrt(gamma*p0/r0)
                    uban1  = ut*ax+vt*ay+wt*az
                    uban0  = u0*ax+v0*ay+w0*az


                    riemann1  = uban1  + 2.0*c1*gm
                    riemann0  = uban0  - 2.0*c0*gm

                    uface = 0.5*(riemann1 + riemann0)
                    aface = 0.25*gm1*(riemann1 - riemann0)


                    !>
                    t17 = uface - uban0
                    !> velocuty on the face
                    t18 = u0 + ax*t17
                    t19 = v0 + ay*t17
                    t20 = w0 + az*t17
!                    write(300+icyc,'(3i3,5e20.12)')i,j,k,t17,t19,ay,uface,uban0
                    !>
                    ent0 = gamma*p0/r0**gamma
                    t9   = ent0
                    t12  = 1.e0/(iblock%cell(i,j,k)%r**gm1)
                    !>
                    t13  = uface + 0.e0
!                    if(real(t13) .ge. 0.e0) write(*,'("t17-old",2e36.24)') uface-uban1,t17
                    t17  = logical_depending(uface-uban1,t17,(real(t13) .ge. 0.e0))
!                    if(real(t13) .ge. 0.e0) write(*,'("t17-new",2e36.24)') uface-uban1,t17
!                    if(real(t13) .ge. 0.e0) write(*,'("t18-old",2e36.24)') u(i,j,k)+ax*t17,t18
                    t18  = logical_depending(iblock%cell(i,j,k)%u+ax*t17,t18,(real(t13) .ge. 0.e0))
!                    if(real(t13) .ge. 0.e0) write(*,'("t18-new",2e36.24)') u(i,j,k)+ax*t17,t18
!                    if(real(t13) .ge. 0.e0) write(*,'("t19-old",2e36.24)') v(i,j,k)+ay*t17,t19
                    t19  = logical_depending(iblock%cell(i,j,k)%v+ay*t17,t19,(real(t13) .ge. 0.e0))
!                    if(real(t13) .ge. 0.e0) write(*,'("t19-new",2e36.24)') v(i,j,k)+ay*t17,t19
!                    if(real(t13) .ge. 0.e0) write(*,'("t20-old",2e36.24)') w(i,j,k)+az*t17,t20
                    t20  = logical_depending(iblock%cell(i,j,k)%w+az*t17,t20,(real(t13) .ge. 0.e0))
!                    if(real(t13) .ge. 0.e0) write(*,'("t20-new",2e36.24)') w(i,j,k)+az*t17,t20
!                    if(real(t13) .ge. 0.e0) write(*,'("t9-old",2e36.24)') c1*c1*t12,t9
                    t9   = logical_depending(c1*c1*t12,t9,(real(t13) .ge. 0.e0))
!                    if(real(t13) .ge. 0.e0) write(*,'("t9-new",2e36.24)') c1*c1*t12,t9

                    !>
                    !>
                    t16 = aface*aface
                    iblock%cell(id1,jd1,kd1)%r  =   (t16/t9)**(1.0/gm1)
                    iblock%cell(id1,jd1,kd1)%u  =   t18
                    iblock%cell(id1,jd1,kd1)%v  =   t19
                    iblock%cell(id1,jd1,kd1)%w  =   t20
                    iblock%cell(id1,jd1,kd1)%p  =   iblock%cell(id1,jd1,kd1)%r*t16/gamma
!                    write(112,'(3i5,8e24.16)') id1,jd1,kd1,iblock%cell(id1,jd1,kd1)%r,iblock%cell(id1,jd1,kd1)%u,&
!                    iblock%cell(id1,jd1,kd1)%v,iblock%cell(id1,jd1,kd1)%w,iblock%cell(id1,jd1,kd1)%p,ax,ay,az
                    !>
                    !>
                    !>
                    !> need to store the viscous of turbulent to computing the viscous coefficient
                    !>
                    if(nvisc .ge. 2)then
                        iblock%turbulent(id1,jd1,kd1)%viscous   =   iblock%turbulent(i,j,k)%viscous
                    end if
                    !> only need to do advanced model turbulence B.C.s on finest grid
                    if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
                        uban =  (ax*iblock%cell(id1,jd1,kd1)%u + &
                                 ay*iblock%cell(id1,jd1,kd1)%v + &
                                 az*iblock%cell(id1,jd1,kd1)%w)*ibc%norm_index
                        if(uban .lt. 0.0)then
                            iblock%turbulent(id1,jd1,kd1)%tur_save = chi
                        else
                            iblock%turbulent(id1,jd1,kd1)%tur_save = iblock%turbulent(i,j,k)%tur_save
                        end if
                    end if
                    !>
                    !>
                    iblock%cell(id2,jd2,kd2)%r  =   iblock%cell(id1,jd1,kd1)%r
                    iblock%cell(id2,jd2,kd2)%u  =   iblock%cell(id1,jd1,kd1)%u
                    iblock%cell(id2,jd2,kd2)%v  =   iblock%cell(id1,jd1,kd1)%v
                    iblock%cell(id2,jd2,kd2)%w  =   iblock%cell(id1,jd1,kd1)%w
                    iblock%cell(id2,jd2,kd2)%p  =   iblock%cell(id1,jd1,kd1)%p
                    !>
                    !>
                    if(iblock%cell(id2,jd2,kd2)%r < epson) iblock%cell(id2,jd2,kd2)%r = epson
                    if(iblock%cell(id2,jd2,kd2)%p < epson) iblock%cell(id2,jd2,kd2)%p = epson

                    !>
                    !> need to store the viscous of turbulent to computing the viscous coefficient
                    !>
                    if(nvisc .ge. 2)then
                         iblock%turbulent(id2,jd2,kd2)%viscous  =   0.0
                    end if

                end do
            end do
        end do
        !>
        !>
        return
    end subroutine bc_far

