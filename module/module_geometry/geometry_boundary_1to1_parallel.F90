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
	!> shared the delta_q with the processors
    !> must shared the 1-1 block boundarary information
    !> to other blocks
    !>
    subroutine bc_1to1_parallel(level)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use bc_module
        !>
        !>
        use nodes_paras
        use nodes_var
        use nodes_var_bc
        use nodes_mg
		!>
		!>
		implicit none
#if defined PMPI
        include "mpif.h"
#endif
        integer :: level
#if defined PMPI
        !>
        !>
        type(bc_types),pointer     :: ibc
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !>
        !>
        integer :: ib,ib_cut,&
                   i,j,k,&
                   dblocks,&
                   sblocks,&
                   id1,jd1,kd1,&
                   id2,jd2,kd2
        integer :: iireq1,&
                   iireq2,&
                   i_check,&
                   ii_check,&
                   ii_check_1,&
                   i_point,&
                   ii_point
        integer :: ierr
        integer :: istat2(mpi_status_size,num_bc)
        integer :: ndone,&
                  nrecvd,&
                  nnn,&
                  ilength ,&
                  iilength,&
                  mytag
        !>
        !>
        integer :: is,ie,&
                   js,je,&
                   ks,ke,&
                   is_cut,ie_cut,&
                   js_cut,je_cut,&
                   ks_cut,ke_cut,&
                   istep,&
                   jstep,&
                   kstep,&
                   istep_cut,&
                   jstep_cut,&
                   kstep_cut,&
                   iid,jjd,kkd
        integer::  ii,ip,&
                   jj,jp,&
                   kk,kp
        integer::  ii_send,&
                   ii_recv,&
                   ii_wait,&
                   ii_bc
        !>
        character(len=20) :: iface
        !>
        !>
#if defined single
		integer,parameter :: mpi_own_real = mpi_real
#else
		integer,parameter :: mpi_own_real = mpi_double_precision
#endif
        mesh => grids(imesh)
        !*************************************************************
        !>
        !>
        !> the blocks and boundary face at some processor
        !> can share the bc informations
        !>
        do  ii_bc = 1,num_bc
            !>
            ibc => mesh%bcs(ii_bc)
            if(ibc%bc_type .eq. 'cut')then
                ib     = ibc%block_index
                ib_cut = ibc%iblock_index
                !>
                !>
                if(level_mg(ib) .eq. level .and. level_mg(ib_cut) .eq. level)then
                    sblocks = n2p(ib)
                    dblocks = n2p(ib_cut)
                    if(myid .eq. sblocks .and. myid .eq. dblocks)then
                        call bc_cut(ii_bc)
                    end if
                end if
            end if
        end do

        !>
        !>
        !>

        !>*************************************************************
        !> send the 1-1 bc informations from source processor
        !> to dest processor
        !> and used isend api to fast send the info
        !>
        iireq1  = 0
        iireq2  = 0
        i_check = 0
        i_point = 1
        !>
        !>
        do ii_send=1,num_bc
            ibc => mesh%bcs(ii_send)
            if(ibc%bc_type .eq. 'cut')then
                !>
                is = ibc%istart
                ie = ibc%iend
                js = ibc%jstart
                je = ibc%jend
                ks = ibc%kstart
                ke = ibc%kend
                is_cut = ibc%istart_cut
                ie_cut = ibc%iend_cut
                js_cut = ibc%jstart_cut
                je_cut = ibc%jend_cut
                ks_cut = ibc%kstart_cut
                ke_cut = ibc%kend_cut
                !>
                !>
                ib      = ibc%block_index
                ib_cut  = ibc%iblock_index
                if(level_mg(ib) .eq. level .and. level_mg(ib_cut) .eq. level)then
                    sblocks = n2p(ib)
                    dblocks = n2p(ib_cut)
                    !>
                    !>
                    if(dblocks .eq. myid .and. sblocks .ne. myid)then
                        !>
                        !>
                        iblock =>mesh%blocks(ib_cut)
                        !>
                        !>
                        if(is == ie )then
                            iface = 'i'
                        else if(js == je)then
                            iface = 'j'
                        else if(ks == ke)then
                            iface = 'k'
                        end if
                        !>
                        !>
                        i_check = i_point
                        !>
                        istep     = sign(1,ie -is)
                        jstep     = sign(1,je -js)
                        kstep     = sign(1,ke -ks)
                        istep_cut = sign(1,ie_cut -is_cut)
                        jstep_cut = sign(1,je_cut -js_cut)
                        kstep_cut = sign(1,ke_cut -ks_cut)
                        !>***********************************************
                        !>***********************************************
                        !>***********************************************
!                        write(182,'("send=",7i4," send to ",i4," from=",i4)') ii_send,is_cut,ie_cut,js_cut,je_cut,ks_cut,ke_cut,sblocks,dblocks
!                        write(182,*) '*************'
                        do k=ks,ke,kstep
                            do j=js,je,jstep
                                do i=is,ie,istep
                                    iid = (i - is)*istep
                                    jjd = (j - js)*jstep
                                    kkd = (k - ks)*kstep
                                    !>
                                    if(ibc%irank .eq. 1)then
                                        ii = is_cut + iid*istep_cut
                                        ip = ii
                                    else if(ibc%irank .eq. 2)then
                                        jj = js_cut + iid*jstep_cut
                                        jp = jj
                                    else if(ibc%irank .eq. 3)then
                                        kk = ks_cut + iid*kstep_cut
                                        kp = kk
                                    end if
                                    !>
                                    if(ibc%jrank .eq. 1)then
                                        ii = is_cut + jjd*istep_cut
                                        ip = ii
                                    else if(ibc%jrank .eq. 2)then
                                        jj = js_cut + jjd*jstep_cut
                                        jp = jj
                                    else if(ibc%jrank .eq. 3)then
                                        kk = ks_cut + jjd*kstep_cut
                                        kp = kk
                                    end if
                                    !>
                                    if(ibc%krank .eq. 1)then
                                        ii = is_cut + kkd*istep_cut
                                        ip = ii
                                    else if(ibc%krank .eq. 2)then
                                        jj = js_cut + kkd*jstep_cut
                                        jp = jj
                                    else if(ibc%krank .eq. 3)then
                                        kk = ks_cut + kkd*kstep_cut
                                        kp = kk
                                    end if
                                    !>
                                    !>
                                    id1 = i
                                    jd1 = j
                                    kd1 = k
                                    id2 = i
                                    jd2 = j
                                    kd2 = k
                                    !>
                                    if(iface .eq. 'i')then
                                        id1 = i + 1*ibc%direction
                                        id2 = i + 2*ibc%direction
                                        if(ibc%irank  .eq. 1)then
                                            ip =ip + sign(1,3-ip)
                                        else if(ibc%irank  .eq. 2)then
                                            jp =jp + sign(1,3-jp)
                                        else if(ibc%irank  .eq. 3)then
                                            kp =kp + sign(1,3-kp)
                                        end if
                                    else if(iface .eq. 'j')then
                                        jd1 = j + 1*ibc%direction
                                        jd2 = j + 2*ibc%direction
                                        if(ibc%jrank  .eq. 1)then
                                            ip =ip + sign(1,3-ip)
                                        else if(ibc%jrank  .eq. 2)then
                                            jp =jp + sign(1,3-jp)
                                        else if(ibc%jrank  .eq. 3)then
                                            kp =kp + sign(1,3-kp)
                                        end if
                                    else if(iface .eq. 'k')then
                                        kd1 = k + 1*ibc%direction
                                        kd2 = k + 2*ibc%direction
                                        if(ibc%krank  .eq. 1)then
                                            ip =ip + sign(1,3-ip)
                                        else if(ibc%krank  .eq. 2)then
                                            jp =jp + sign(1,3-jp)
                                        else if(ibc%krank  .eq. 3)then
                                            kp =kp + sign(1,3-kp)
                                        end if
                                    end if
                                    !>
                                    !>

                                    !write(182,'("send=",7i4)') ii_send,ip,jp,kp,ii,jj,kk
                                    wk(i_point) = iblock%cell(ii,jj,kk)%r
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%cell(ii,jj,kk)%u
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%cell(ii,jj,kk)%v
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%cell(ii,jj,kk)%w
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%cell(ii,jj,kk)%p
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%metric(ii,jj,kk)%volume
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%metric(ii,jj,kk)%jacobian
                                    i_point     = i_point + 1
                                    !>
                                    !> need to store the viscous of turbulent to computing the viscous coefficient
                                    !>
                                    if(nvisc .ge. 2)then
                                        wk(i_point) = iblock%turbulent(ii,jj,kk)%viscous
                                        i_point     = i_point + 1
                                    end if
                                    !>
                                    !>
                                    !> only need to do advanced model turbulence B.C.s on finest grid
                                    if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
                                        wk(i_point) = iblock%turbulent(ii,jj,kk)%tur_save
                                        i_point     = i_point + 1
                                    end if
                                    !>
                                    !>
                                    wk(i_point) = iblock%cell(ip,jp,kp)%r
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%cell(ip,jp,kp)%u
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%cell(ip,jp,kp)%v
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%cell(ip,jp,kp)%w
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%cell(ip,jp,kp)%p
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%metric(ip,jp,kp)%volume
                                    i_point     = i_point + 1
                                    wk(i_point) = iblock%metric(ip,jp,kp)%jacobian
                                    i_point     = i_point + 1
                                    !if(icyc ==5 .and. ii_send==58 ) write(10151,'(6i3,7e24.16)') ii,jj,kk,ip,jp,kp,iblock%cell(ii,jj,kk)%r,iblock%cell(ii,jj,kk)%u,iblock%cell(ii,jj,kk)%v,&
                                    !iblock%cell(ii,jj,kk)%w,iblock%cell(ii,jj,kk)%p,iblock%turbulent(ii,jj,kk)%viscous,iblock%turbulent(ii,jj,kk)%tur_save
                                    !if(ii_send==81 .and. icyc==2) write(101511,'(6i3,7e24.16)') ii,jj,kk,ip,jp,kp,iblock%cell(ii,jj,kk)%r,iblock%cell(ii,jj,kk)%u,iblock%cell(ii,jj,kk)%v,&
                                    !iblock%cell(ii,jj,kk)%w,iblock%cell(ii,jj,kk)%p,iblock%turbulent(ii,jj,kk)%viscous,iblock%turbulent(ii,jj,kk)%tur_save
                                end do
                            end do
                        end do
                        !>***********************************************
                        !>***********************************************
                        !>***********************************************
                        ilength = i_point  - i_check
                        iireq2  = iireq2 + 1
                        mytag   = ii_send
                        !>
                        !>
!                        if(ib_cut ==7 .or. ib_cut==8) write(*,'("send myid=",i4," nbl=",i4," to =",i4," local=",i4, " buffer size=",i16, "and the ireq=",i6," tag=",i4)')myid,ib_cut,sblocks,dblocks,ilength,iireq2,mytag
                        call mpi_isend(wk(i_check),&
                                       ilength,&
                                       !mpi_own_real,&
                                       mpi_double_precision,&
                                       sblocks,&
                                       mytag,&
                                       mycomm,&
                                       ireq_s(iireq2),&
                                       ierr)
                    end if
                end if
            end if
        end do
        !>
        !>
        !>
        !****************************************************************\
        !> recive the info from the source processor
        !>
        !>
        do ii_recv =1 ,num_bc
            ibc => mesh%bcs(ii_recv)
            !>
            !>
            if(ibc%bc_type .eq. 'cut')then
                !>
                ib = ibc%block_index
                ib_cut = ibc%iblock_index
                if(level_mg(ib) .eq. level .and. level_mg(ib_cut) .eq. level)then
                    sblocks = n2p(ib)
                    dblocks = n2p(ib_cut)
                    if(myid .eq. sblocks  .and. myid .ne. dblocks)then
                        is = ibc%istart
                        ie = ibc%iend
                        js = ibc%jstart
                        je = ibc%jend
                        ks = ibc%kstart
                        ke = ibc%kend
                        istep = sign(1,ie-is)
                        jstep = sign(1,je-js)
                        kstep = sign(1,ke-ks)
                        ii_check = i_point
                        !>
                        do k =ks,ke,kstep
                            do j=js,je,jstep
                                do i=is,ie,istep
                                    i_point = i_point + 14
                                    if(nvisc .ge. 2)then
                                        i_point = i_point + 1
                                    end if
                                    !>
                                    !>
                                    !> only need to do advanced model turbulence B.C.s on finest grid
                                    if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
                                        i_point = i_point + 1
                                    end if
                                end do
                            end do
                        end do
                        !>
                        !>
                        iilength = i_point - ii_check
                        iireq1   = iireq1 + 1
                        mytag    = ii_recv
!                        if( ib==8) write(*,'("recv myid=",i4," nbl=",i4," from =",i4," local=",i4, " buffer size=",i16, "and the ireq=",i6," tag=",i4)')myid,ib,sblocks,dblocks,iilength,iireq1,mytag
                        !>
                        !if(ii_recv == 58) write(*,'(2i10)') ii_check,iilength
                        call mpi_irecv(wk(ii_check),&
                                       iilength,&
                                       !mpi_own_real,&
                                       mpi_double_precision,&
                                       dblocks,&
                                       mytag,&
                                       mycomm,&
                                       ireq_r(iireq1),&
                                       ierr)
                        !>
                        !>
                        keep1(iireq1) = ii_check
                        keep2(iireq1) = ii_recv
                        !if(myid ==10 )write(*,'("myid=",i4," send buffer=",2i10)') myid,ii_recv,iilength
                        !>
                        !>
                    end if
                end if
                !>
                !>
            end if
        end do
        !>
        !>
        !**************************************************************
        !> bound_sendrecv
        !> used the bc info before the info from the source processor
        !> have recv already
        !>
        ndone  = 0
        nrecvd = 0
        do while(ndone .lt. iireq1)
            !>
            !>
            !call  mpi_testsome(iireq1,ireq_r,nrecvd,index_r,istat2,ierr)
            call  mpi_waitsome(iireq1,ireq_r,nrecvd,index_r,istat2,ierr)
            !if(myid ==2) write(*,'("in waitsome=",i8)') nrecvd
            !>
            !>
            if(nrecvd > 0)then
                !>
                !>
                ndone = ndone + nrecvd
                !>
                do  nnn =1,nrecvd
                    ii_wait   = keep2(index_r(nnn))
                    ibc => mesh%bcs(ii_wait)
                    !>
                    ib      = ibc%block_index
                    ib_cut  = ibc%iblock_index
                    !>
                    iblock => mesh%blocks(ib)
                    sblocks = n2p(ib)
                    dblocks = n2p(ib_cut)
                    !if(ib==7 .or. ib==8 )write(*,'("waitsome=",2i8)') ii_wait,index_r(nnn)

                    if(myid == sblocks .and. myid /= dblocks)then
                        !>
                        !>
                        ii_check_1 = keep1(index_r(nnn))
                        ii_point   = ii_check_1
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
                        !write(182,'("recv=",7i4)') ii_wait,is,ie,js,je,ks,ke
                        !write(182,*) '*************'
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

                                    !write(182,'("recv=",7i4," recv to ",i4," from=",i4)') ii_wait,id1,jd1,kd1,id2,jd2,kd2,sblocks,dblocks
                                    !>
                                    iblock%cell(id1,jd1,kd1)%r  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%cell(id1,jd1,kd1)%u  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%cell(id1,jd1,kd1)%v  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%cell(id1,jd1,kd1)%w  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%cell(id1,jd1,kd1)%p  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%metric(id1,jd1,kd1)%volume  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%metric(id1,jd1,kd1)%jacobian= wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    !>
                                    !> need to store the viscous of turbulent to computing the viscous coefficient
                                    !>
                                    if(nvisc .ge. 2)then
                                        iblock%turbulent(id1,jd1,kd1)%viscous = wk(ii_point)
                                        ii_point                              = ii_point + 1
                                    end if
                                    !>
                                    !>
                                    !> only need to do advanced model turbulence B.C.s on finest grid
                                    if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
                                        iblock%turbulent(id1,jd1,kd1)%tur_save = wk(ii_point)
                                        ii_point                               = ii_point + 1
                                    end if
                                    !>
                                    iblock%cell(id2,jd2,kd2)%r  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%cell(id2,jd2,kd2)%u  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%cell(id2,jd2,kd2)%v  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%cell(id2,jd2,kd2)%w  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%cell(id2,jd2,kd2)%p  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%metric(id2,jd2,kd2)%volume  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    iblock%metric(id2,jd2,kd2)%jacobian  = wk(ii_point)
                                    ii_point                    = ii_point + 1
                                    if(id1 > iblock%idim+1 .or. jd1 > iblock%jdim+1 .or. kd1 > iblock%kdim+1 .or. ii_point-1 .gt. i_size1to1)then
                                        write(*,'("wait nbl=",i6," bcs=",i6," face1=(",3i4,") face2(",3i4,") and block axis=(",3i4,") size=",i16)') ib,ii_wait,&
                                        id1,jd1,kd1,id2,jd2,kd2,iblock%idim+2,iblock%jdim+2,iblock%kdim+2,i_size1to1-ii_point+1
                                    end if
                                    !>
                                    !>
                                end do
                            end do
                        end do
                        !>
                    end if
                end do
            end if
        end do
        !>
        !> make sure all the sends are completed before the subroutine exiting
        !>
        !>
        if(iireq2 .gt. 0 )then
           call mpi_waitall(iireq2,ireq_s,istat2,ierr)
        end if
#endif
        return
    end subroutine bc_1to1_parallel
    !**********************************************
    !**********************************************
