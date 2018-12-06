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
    subroutine bc_cut(num_of_bc)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use bc_module
        !>
        !>
        use nodes_paras
        use nodes_var_bc
        use nodes_var
        !>
        !>
		implicit none
		!>
		!>
        type(blocks_type),pointer  :: iblock,&
                                      iblock_cut
        type(bc_types),pointer     :: ibc
        type(overlap_type),pointer :: mesh
        !>
        integer :: num_of_bc,&
                   sblock,dblock,&
                   is,ie,&
                   js,je,&
                   ks,ke,&
                   is_cut,ie_cut,&
                   js_cut,je_cut,&
                   ks_cut,ke_cut,&
                   i,j,k,&
                   istep,&
                   jstep,&
                   kstep,&
                   istep_cut,&
                   jstep_cut,&
                   kstep_cut,&
                   iid,jjd,kkd,&
                   ip,jp,kp,&
                   ii,jj,kk,&
                   id1,jd1,kd1,&
                   id2,jd2,kd2
        character(len=20) :: iface
        !>
        !>
        !*************************************************************
        !>
        mesh => grids(imesh)
        ibc  => mesh%bcs(num_of_bc)
        !>
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
        is_cut = ibc%istart_cut
        ie_cut = ibc%iend_cut
        js_cut = ibc%jstart_cut
        je_cut = ibc%jend_cut
        ks_cut = ibc%kstart_cut
        ke_cut = ibc%kend_cut
        !>
        !>
        if( is .eq. ie)then
            iface = 'i'
        else if(js .eq. je)then
            iface = 'j'
        else if(ks .eq. ke)then
            iface = 'k'
        end if
        !>
        istep = sign(1,ie - is)
        jstep = sign(1,je - js)
        kstep = sign(1,ke - ks)
        !>
        !>
        istep_cut = sign(1,ie_cut -is_cut)
        jstep_cut = sign(1,je_cut -js_cut)
        kstep_cut = sign(1,ke_cut -ks_cut)
        !>
        !>
        dblock = ibc%iblock_index
        sblock = ibc%block_index
        !>
        !>
        iblock     => mesh%blocks(sblock)
        iblock_cut => mesh%blocks(dblock)
        !>
        !>
        do k= ks,ke,kstep
            do j=js,je,jstep
                do i=is,ie,istep
                    !>
                    !>
                    iid =   (i-is)*istep
                    jjd =   (j-js)*jstep
                    kkd =   (k-ks)*kstep
                    !>
                    !>
                    !>
                    if(ibc%irank  .eq. 1)then
                        ii = is_cut + iid * istep_cut
                        ip = ii
                    else if(ibc%irank  .eq. 2)then
                        jj = js_cut + iid * jstep_cut
                        jp = jj
                    else if(ibc%irank  .eq. 3)then
                        kk = ks_cut + iid * kstep_cut
                        kp = kk
                    end if
                    !>
                    !>
                    !>
                    if(ibc%jrank  .eq. 1)then
                        ii = is_cut + jjd * istep_cut
                        ip = ii
                    else if(ibc%jrank  .eq. 2)then
                        jj = js_cut + jjd * jstep_cut
                        jp = jj
                    else if(ibc%jrank  .eq. 3)then
                        kk = ks_cut + jjd * kstep_cut
                        kp = kk
                    end if
                    !>
                    !>
                    !>
                    if(ibc%krank  .eq. 1)then
                        ii = is_cut + kkd * istep_cut
                        ip = ii
                    else if(ibc%krank  .eq. 2)then
                        jj = js_cut + kkd * jstep_cut
                        jp = jj
                    else if(ibc%krank  .eq. 3)then
                        kk = ks_cut + kkd * kstep_cut
                        kp = kk
                    end if
                    !>
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
                    !>
                    !>
                    if(iface .eq. 'i')then
                        !>
                        !>
                        id1 = i + 1*ibc%direction
                        id2 = i + 2*ibc%direction
                        if(ibc%irank  .eq. 1)then
                            ip =ip + sign(1,3-ip)
                        else if(ibc%irank  .eq. 2)then
                            jp =jp + sign(1,3-jp)
                        else if(ibc%irank  .eq. 3)then
                            kp =kp + sign(1,3-kp)
                        end if
                        !>
                        !>
                    else if(iface .eq. 'j')then
                        !>
                        !>
                        jd1 = j + 1*ibc%direction
                        jd2 = j + 2*ibc%direction
                        if(ibc%jrank  .eq. 1)then
                            ip =ip + sign(1,3-ip)
                        else if(ibc%jrank  .eq. 2)then
                            jp =jp + sign(1,3-jp)
                        else if(ibc%jrank  .eq. 3)then
                            kp =kp + sign(1,3-kp)
                        end if
                        !>
                        !>
                    else if(iface .eq. 'k')then
                        !>
                        !>
                        kd1 = k + 1*ibc%direction
                        kd2 = k + 2*ibc%direction
                        if(ibc%krank  .eq. 1)then
                            ip =ip + sign(1,3-ip)
                        else if(ibc%krank  .eq. 2)then
                            jp =jp + sign(1,3-jp)
                        else if(ibc%krank  .eq. 3)then
                            kp =kp + sign(1,3-kp)
                        end if
                        !>
                        !>
                    end if
                    !>
                    !>
                    !>
!                    write(181,'("num_of_bc=",13i4)') num_of_bc,is,ie,js,je,ks,ke,is_cut,ie_cut,js_cut,je_cut,ks_cut,ke_cut
!                    write(181,'("num_of_bc=",i4,"(",3i4,")<=(",3i4,") and (",3i4,")<=(",3i4,")")') num_of_bc,id1,jd1,kd1,ii,jj,kk,id2,jd2,kd2,ip,jp,kp
                    iblock%cell(id1,jd1,kd1)%r = iblock_cut%cell(ii,jj,kk)%r
                    iblock%cell(id1,jd1,kd1)%p = iblock_cut%cell(ii,jj,kk)%p
                    iblock%cell(id1,jd1,kd1)%u = iblock_cut%cell(ii,jj,kk)%u
                    iblock%cell(id1,jd1,kd1)%v = iblock_cut%cell(ii,jj,kk)%v
                    iblock%cell(id1,jd1,kd1)%w = iblock_cut%cell(ii,jj,kk)%w
                    iblock%metric(id1,jd1,kd1)%volume = iblock_cut%metric(ii,jj,kk)%volume
                    iblock%metric(id1,jd1,kd1)%jacobian = iblock_cut%metric(ii,jj,kk)%jacobian
                    !>
                    !>
                    !>
                    if(iblock%cell(id1,jd1,kd1)%r < 1.e-10 .or. iblock%cell(id1,jd1,kd1)%p <1.e-10)then
                        open(unit=89,file='ccfd.bc_error',status='unknown')
                        write(89,'("iter=",i5," bc_num=",i8," at the bc_some1 block=",i6,"  rho is too small r1(",i3,",",i3,",",i3,")=",&
                        e16.10," p1(id1,jd1,kd1)=",e16.10)') icyc,num_of_bc,ibc%block_index,id1,jd1,kd1,iblock%cell(id1,jd1,kd1)%r ,iblock%cell(id1,jd1,kd1)%p
                        write(89,'("iter=",i5," bc_num=",i8," at the bc_some1 block=",i6,"  rho is too small r01(",i3,",",i3,",",i3,")=",&
                        e16.10," p01(i,j,k)=",e16.10)') icyc,num_of_bc,ibc%block_index,ii,jj,kk,iblock_cut%cell(ii,jj,kk)%r,iblock_cut%cell(ii,jj,kk)%p
                        write(89,*) '-------------------------------------------------------------------------------------'
                    endif
                    !>
                    !> need to store the viscous of turbulent to computing the viscous coefficient
                    !>
                    if(nvisc .ge. 2)then
                        iblock%turbulent(id1,jd1,kd1)%viscous = iblock_cut%turbulent(ii,jj,kk)%viscous
                        !iblock%turbulent(id1,jd1,kd1)%viscous = iblock%turbulent(i,j,k)%viscous
                    end if
                    !>
                    !>
                    !> only need to do advanced model turbulence B.C.s on finest grid
                    if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
                        iblock%turbulent(id1,jd1,kd1)%tur_save = iblock_cut%turbulent(ii,jj,kk)%tur_save
                        !iblock%turbulent(id1,jd1,kd1)%tur_save = iblock%turbulent(i,j,k)%tur_save
                    end if
                    !>
                    !>
                    iblock%cell(id2,jd2,kd2)%r = iblock_cut%cell(ip,jp,kp)%r
                    iblock%cell(id2,jd2,kd2)%p = iblock_cut%cell(ip,jp,kp)%p
                    iblock%cell(id2,jd2,kd2)%u = iblock_cut%cell(ip,jp,kp)%u
                    iblock%cell(id2,jd2,kd2)%v = iblock_cut%cell(ip,jp,kp)%v
                    iblock%cell(id2,jd2,kd2)%w = iblock_cut%cell(ip,jp,kp)%w
                    iblock%metric(id2,jd2,kd2)%volume = iblock_cut%metric(ip,jp,kp)%volume
                    iblock%metric(id2,jd2,kd2)%jacobian = iblock_cut%metric(ip,jp,kp)%jacobian
                    !write(*,'("num_bc=",i8,"  1-face:",3i4," <- ",3i4,"   2-face:",3i4," <- ",3i4)') num_of_bc,id1,jd1,kd1,ii,jj,kk,id2,jd2,kd2,ip,jp,kp
                    !>
                    !>
                    !>
                    !>
                    if(iblock%cell(id2,jd2,kd2)%r <1.e-10 .or. iblock%cell(id2,jd2,kd2)%p <1.e-10)then
                        open(unit=89,file='ccfd.bc_error',status='unknown')
                        write(89,'("iter=",i5," bc_num=",i8," at the bc_some2 block=",i6,"  rho is too small r2(",i3,",",i3,",",i3,")=",&
                        e16.10," p2(id1,jd1,kd1)=",e16.10)') icyc,num_of_bc,ibc%block_index,id2,jd2,kd2,iblock%cell(id2,jd2,kd2)%r,iblock%cell(id2,jd2,kd2)%p
                        write(89,'("iter=",i5," bc_num=",i8," at the bc_some2 block=",i6,"  rho is too small r02(",i3,",",i3,",",i3,")=",&
                        e16.10," p02(i,j,k)=",e16.10)') icyc,num_of_bc,ibc%iblock_index,ip,jp,kp,iblock_cut%cell(ip,jp,kp)%r,iblock_cut%cell(ip,jp,kp)%p
                        write(89,*) '-------------------------------------------------------------------------------------'
                    endif
                    !>

                    !>
                end do
            end do
        end do
        !>
        !>
        !>
        return
        !>
    end subroutine bc_cut

      !**********************************************
      !**********************************************
      !**********************************************
