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
    subroutine bc_symmetry(num_of_bc,nbl)
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use bc_module
		!>
		!>
		implicit none
        !>
        !>
        integer :: num_of_bc,nbl
        integer :: is,ie,&
                   js,je,&
                   ks,ke,&
                   istep,&
                   jstep,&
                   kstep

        integer :: id1,jd1,kd1,&
                   id2,jd2,kd2,&
                   ip,jp,kp,&
                   i,j,k
        character(len =20 ) :: iface
        !>
        type(blocks_type),pointer  :: iblock
        type(bc_types),pointer     :: ibc
        type(overlap_type),pointer :: mesh
        real(kind = dprec) :: agx,agy,agz
        real(kind = dprec) :: uu,vv,ww,uban,qq
        !>
        !>
        mesh   => grids(imesh)
        ibc    => mesh%bcs(num_of_bc)
        iblock => mesh%blocks(nbl)
        !>
        !>
        is = ibc%istart
        ie = ibc%iend
        js = ibc%jstart
        je = ibc%jend
        ks = ibc%kstart
        ke = ibc%kend

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
        istep = sign(1,ie-is)
        jstep = sign(1,je-js)
        kstep = sign(1,ke-ks)

        do i = is,ie,istep
            do k = ks,ke,kstep
                do j = js,je,jstep
                    !>
                    !>
                    id1 = i
                    jd1 = j
                    kd1 = k
                    id2 = i
                    jd2 = j
                    kd2 = k
                    !>
                    !>
                    ip  = i
                    jp  = j
                    kp  = k
                    !>
                    !>

                    if (iface .eq. 'i') then
                        !>
                        !>
                        id1 = i+1*ibc%direction
                        id2 = i+2*ibc%direction
                        !>
                        !>
                        ip  = i+sign(1,3-i)
                        !>
                        !>
                        agx = iblock%metric(id1,jd1,kd1)%fi
                        agy = iblock%metric(id1,jd1,kd1)%fj
                        agz = iblock%metric(id1,jd1,kd1)%fk
                    else if (iface .eq. 'j') then
                        !>
                        !>
                        jd1 = j+1*ibc%direction
                        jd2 = j+2*ibc%direction
                        !>
                        !>
                        jp  = j+sign(1,3-j)
                        !>
                        !>
                        agx = iblock%metric(id1,jd1,kd1)%gi
                        agy = iblock%metric(id1,jd1,kd1)%gj
                        agz = iblock%metric(id1,jd1,kd1)%gk
                    else if (iface .eq. 'k') then
                        !>
                        !>
                        kd1 = k+1*ibc%direction
                        kd2 = k+2*ibc%direction
                        !>
                        !>
                        kp  = k+sign(1,3-k)
                        !>
                        !>
                        agx = iblock%metric(id1,jd1,kd1)%hi
                        agy = iblock%metric(id1,jd1,kd1)%hj
                        agz = iblock%metric(id1,jd1,kd1)%hk
                    end if
                    !>
                    !>
                    uu   = iblock%cell(i,j,k)%u
                    vv   = iblock%cell(i,j,k)%v
                    ww   = iblock%cell(i,j,k)%w
                    !>
                    !>
                    uban = uu*agx + vv*agy + ww*agz
                    !>
                    !>
                    iblock%cell(id1,jd1,kd1)%r = iblock%cell(i,j,k)%r
                    iblock%cell(id1,jd1,kd1)%p = iblock%cell(i,j,k)%p
                    !>
                    !>
                    iblock%cell(id1,jd1,kd1)%u = ( uu-2.0*uban*agx)
                    iblock%cell(id1,jd1,kd1)%v = ( vv-2.0*uban*agy)
                    iblock%cell(id1,jd1,kd1)%w = ( ww-2.0*uban*agz)
                    !>
                    !>
                    !>
                    !> need to store the viscous of turbulent to computing the viscous coefficient
                    !>
!                    if(nbl ==1 )write(110,'("num_of_bc=",4i5,8e24.16)') num_of_bc,id1,jd1,kd1,iblock%cell(id1,jd1,kd1)%r,&
!                    iblock%cell(id1,jd1,kd1)%u,iblock%cell(id1,jd1,kd1)%v,iblock%cell(id1,jd1,kd1)%w,iblock%cell(id1,jd1,kd1)%p,agx,agy,agz
                    if(nvisc .ge. 2)then
                        iblock%turbulent(id1,jd1,kd1)%viscous  = iblock%turbulent(i,j,k)%viscous
                    end if
                    !>
                    !> only need to do advanced model turbulence B.C.s on finest grid
                    if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
                        iblock%turbulent(id1,jd1,kd1)%tur_save = iblock%turbulent(i,j,k)%tur_save
                    end if
                    !>
                    !>
                    if(iblock%cell(id1,jd1,kd1)%r < 1.e-10 .or. iblock%cell(id1,jd1,kd1)%p < 1.e-10)then
                        open(unit=89,file='ccfd.bc_error',status='unknown')
                        write(89,'("iter=",i5," at the bc_sym1 block=",i6,"  rho is too small r(",i3,",",i3,",",i3,")=",&
                        e16.10," p(id1,jd1,kd1)=",e16.10)') icyc,nbl,id1,jd1,kd1,iblock%cell(id1,jd1,kd1)%r,iblock%cell(id1,jd1,kd1)%p
                        write(89,'("iter=",i5," at the bc_sym1 block=",i6,"  rho is too small r(",i3,",",i3,",",i3,")=",&
                        e16.10," p(i,j,k)=",e16.10)') icyc,nbl,i,j,k,iblock%cell(i,j,k)%p,iblock%cell(i,j,k)%p
                        write(89,*) '-------------------------------------------------------------------------------------'

                    end if
                    !>0726
                    uu   = iblock%cell(ip,jp,kp)%u
                    vv   = iblock%cell(ip,jp,kp)%v
                    ww   = iblock%cell(ip,jp,kp)%w
                    !>
                    !>
                    uban = uu*agx + vv*agy + ww*agz
                    !>
                    !>
                    iblock%cell(id2,jd2,kd2)%r = iblock%cell(ip,jp,kp)%r
                    iblock%cell(id2,jd2,kd2)%p = iblock%cell(ip,jp,kp)%p
                    !>
                    !>
                    iblock%cell(id2,jd2,kd2)%u = (uu-2.*uban*agx)
                    iblock%cell(id2,jd2,kd2)%v = (vv-2.*uban*agy)
                    iblock%cell(id2,jd2,kd2)%w = (ww-2.*uban*agz)

!                    write(500+icyc,'(4i3,10e20.12)') num_of_bc,i,j,k,iblock%cell(id1,jd1,kd1)%r,iblock%cell(id1,jd1,kd1)%u,&
!                    iblock%cell(id1,jd1,kd1)%v,iblock%cell(id1,jd1,kd1)%w,iblock%cell(id1,jd1,kd1)%p,&
!                    iblock%cell(id2,jd2,kd2)%r,iblock%cell(id2,jd2,kd2)%u,iblock%cell(id2,jd2,kd2)%v,&
!                    iblock%cell(id2,jd2,kd2)%w,iblock%cell(id2,jd2,kd2)%p
                    !>
                    !> need to store the viscous of turbulent to computing the viscous coefficient
                    !>
                    if(nvisc .ge. 2)then
                        iblock%turbulent(id2,jd2,kd2)%viscous = iblock%turbulent(ip,jp,kp)%viscous
                    end if
                    !>0726
                    if(iblock%cell(id2,jd2,kd2)%r < 1.e-10 .or. iblock%cell(id2,jd2,kd2)%p < 1.e-10)then
                        open(unit=89,file='ccfd.bc_error',status='unknown')
                        write(89,'("iter=",i5," at the bc_sym2 block=",i5,"  rho is too small r(",i3,",",i3,",",i3,")=",&
                        e16.10," p(id1,jd1,kd1)=",e16.10)') icyc,nbl,id2,jd2,kd2,iblock%cell(id2,jd2,kd2)%r,iblock%cell(id2,jd2,kd2)%p
                        write(89,'("iter=",i5," at the bc_sym2 block=",i5,"  rho is too small r(",i3,",",i3,",",i3,")=",&
                        e16.10," p(i,j,k)=",e16.10)') icyc,nbl,ip,jp,kp,iblock%cell(ip,jp,kp)%r,iblock%cell(ip,jp,kp)%p
                        write(89,*) '-------------------------------------------------------------------------------------'
                    !close(69)
                    endif

                    !>
                end do
            end do
        end do
        return
    end subroutine bc_symmetry

