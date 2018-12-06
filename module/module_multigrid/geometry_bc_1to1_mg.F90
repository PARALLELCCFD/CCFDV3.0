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
    subroutine bc_1to1_mg(level,isize_level)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use bc_module
        use cell_module
        use variable_module
        use nodes_paras
        use nodes_var
        use nodes_mg
        use nodes_var_bc
		!>
		!>
		implicit none
#if defined PMPI
        include "mpif.h"
        integer::istat2(mpi_status_size)
#endif
        !>
        !>
        type(overlap_type),pointer   :: mesh
        type(bc_types),pointer       :: ibc
        type(blocks_type),pointer    :: iblock
        !>
        integer :: iireq1,iireq2,num_of_bc,iireq_p
        integer :: sblock,dblock,i_size1,level,isize_level
        integer :: i,j,k,&
                   ii,jj,kk,&
                   ii1,ii2,nd_src,nd_dest,&
                   nbl,&
                   i_point,i_check,ierr,&
                   ii_point,ii_check,ii_check_1,&
                   ndone,nrecvd,nnn,ii_wait,&
                   ilength,iilength,mytag
        integer :: is,js,ks,&
                   ie,je,ke,&
                   is_cut,js_cut,ks_cut,&
                   ie_cut,je_cut,ke_cut,&
                   istep,jstep,kstep,&
                   istep_cut,jstep_cut,kstep_cut,&
                   iid,jjd,kkd,&
                   id1,jd1,kd1,&
                   id2,jd2,kd2,&
                   ip,jp,kp
        !>
        character(len=20) :: iface
        !*************************************************************
        !>
        mesh => grids(imesh)
        !>
        !>if the block the iblock at the same nodes
        !>
        do num_of_bc=1,num_bc
            ibc => mesh%bcs(num_of_bc)
            if(ibc%bc_type .eq. 'cut')then
                ii1 = ibc%block_index
                ii2 = ibc%iblock_index
                if(level_mg(ii1) .eq. level .and. level_mg(ii2) .eq. level)then
                    nd_src  = n2p(ii1)
                    nd_dest = n2p(ii2)
                    if( nd_src .eq. myid .and.  nd_dest .eq. myid)then
                        call bc_cut(num_of_bc)
                    end if
                end if
            end if
        end do
        !>
        !>
#if defined PMPI
        iireq1  = num_bc*(global_level-level)+1
        iireq2  = num_bc*(global_level-level)+1
        iireq_p = 0
        i_point = isize_level
        bc_cut_send :  do num_of_bc=1,num_bc
            ibc => mesh%bcs(num_of_bc)
            if(ibc%bc_type .eq. 'cut')then
                ii1 = ibc%block_index
                ii2 = ibc%iblock_index
                if(level_mg(ii1) .eq. level .and. level_mg(ii2) .eq. level)then
                    nd_src = n2p(ii1)
                    nd_dest =n2p(ii2)
                    if(num_of_bc == 388)then
                        write(*,'("388 bc=",9i4)')num_of_bc,ii1,ii2,level_mg(ii1),level_mg(ii2),level,nd_src,nd_dest,myid
                    end if
                    if( nd_src .ne. myid .and. nd_dest .eq. myid)then
                        write(*,'("in the send proc and bc_num=",i4," myid=",i4)') num_of_bc,myid
                        !>
                        !>
                        is     = ibc%istart
                        ie     = ibc%iend
                        js     = ibc%jstart
                        je     = ibc%jend
                        ks     = ibc%kstart
                        ke     = ibc%kend
                        !>
                        is_cut = ibc%istart_cut
                        ie_cut = ibc%iend_cut
                        js_cut = ibc%jstart_cut
                        je_cut = ibc%jend_cut
                        ks_cut = ibc%kstart_cut
                        ke_cut = ibc%kend_cut
                        !>
                        !>
                        !>
                        i_check  = i_point
                        dblock   = ibc%iblock_index
                        iblock   => mesh%blocks(dblock)
!                        nd_dest  = n2p(ibc%block_index)
                        !>
                        !>
                        if(is .eq. ie )then
                            iface = 'i'
                        else if(js .eq. je )then
                            iface = 'j'
                        else if(ks .eq. ke )then
                            iface = 'k'
                        end if
                        !>
                        istep = sign(1,ie-is)
                        jstep = sign(1,je-js)
                        kstep = sign(1,ke-ks)
                        !>
                        !>
                        istep_cut = sign(1,ie_cut-is_cut)
                        jstep_cut = sign(1,je_cut-js_cut)
                        kstep_cut = sign(1,ke_cut-ks_cut)
                        !>
                        !>
                        write(183,'("send=",7i4)') num_of_bc,is_cut,ie_cut,js_cut,je_cut,ks_cut,ke_cut
                        write(183,*) '*************'
                        do k=ks,ke,kstep
                            do j=js,je,jstep
                                do i=is,ie,istep
                                    !>
                                    !>
                                    iid = (i-is)*istep
                                    jjd = (j-js)*jstep
                                    kkd = (k-ks)*kstep
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
                                    id1 = i
                                    jd1 = j
                                    kd1 = k
                                    id2 = i
                                    jd2 = j
                                    kd2 = k
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

                                    write(183,'("send=",7i4," send to ",i4," from=",i4)') num_of_bc,ip,jp,kp,ii,jj,kk,nd_dest,nd_src

!                                    write(*,'("the axis for send=",7i4)') dblock,ii,jj,kk,ip,jp,kp
                                    !> send the conservation variable of 1-1 block bc to the dest blocks
                                    !> the face of send include idim+1,idim+2,idim+3
                                    !> and i0-2,i0-1,i0
                                    !> the conservation variable at the cell centre
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

                                    !> update the residual from the neighbor blocks
                                    !> and the residual is restrict from the finer mesh
                                    !> just update two cell of bc
                                    !> if is need can update more cell of faces
                                    !>
                                    wk(i_point) =  iblock%variable(ip,jp,kp)%res_1
                                    i_point     = i_point + 1
                                    wk(i_point) =  iblock%variable(ip,jp,kp)%res_2
                                    i_point     = i_point + 1
                                    wk(i_point) =  iblock%variable(ip,jp,kp)%res_3
                                    i_point     = i_point + 1
                                    wk(i_point) =  iblock%variable(ip,jp,kp)%res_4
                                    i_point     = i_point + 1
                                    wk(i_point) =  iblock%variable(ip,jp,kp)%res_5
                                    i_point     = i_point + 1
                                    !>
                                end do
                            end do
                        end do
                        !>
                        !>
                        !>
                        ilength = i_point - i_check
                        iireq2  = iireq2 + 1
                        mytag   = num_of_bc
!                        write(*,'("send myid=",2i4," length=",2i8)') myid,nd_dest,ilength,i_point
                        call mpi_isend(wk(i_check),&
                                       ilength,&
                                       mpi_double_precision,&
                                       nd_src,&
                                       mytag,&
                                       mycomm,&
                                       ireq_s(iireq2),&
                                       ierr)
                    end if
                end if
            end if
        enddo bc_cut_send
        call mpi_barrier(mycomm,ierr)
        !>
        !>
        !****************************************************************
        bc_cut_recv : do num_of_bc=1,num_bc
            ibc => mesh%bcs(num_of_bc)
            if(ibc%bc_type .eq. 'cut')then
                ii1 = ibc%block_index
                ii2 = ibc%iblock_index
                if(level_mg(ii1) .eq. level .and. level_mg(ii2) .eq. level)then
                    nd_src  = n2p(ii1)
                    nd_dest = n2p(ii2)
                    if( myid .eq. nd_src .and. myid .ne. nd_dest)then
                        !>
                        !>
                        ii_check = i_point
!                        nd_src = n2p(ibc%iblock_index)
                        is = ibc%istart
                        ie = ibc%iend
                        js = ibc%jstart
                        je = ibc%jend
                        ks = ibc%kstart
                        ke = ibc%kend
                        istep = sign(1,ie-is)
                        jstep = sign(1,je-js)
                        kstep = sign(1,ke-ks)
                        !>
                        do k =ks,ke,kstep
                            do j=js,je,jstep
                                do i=is,ie,istep
                                    i_point = i_point + 15
                                end do
                            end do
                        end do
                        !>
                        !>
                        iilength = i_point - ii_check
                        iireq1  = iireq1 + 1
                        iireq_p = iireq_p + 1
                        mytag   = num_of_bc
                        !>
                        !>
                        call mpi_irecv(wk(ii_check),&
                                       iilength,&
                                       mpi_double_precision,&
                                       nd_dest,&
                                       mytag,&
                                       mycomm,&
                                       ireq_r(iireq1),&
                                       ierr)
                        keep1(iireq1) = ii_check
                        keep2(iireq1) = num_of_bc
                    end if
                end if
            end if
        enddo bc_cut_recv
        !>
        !>
        !>
        !> bound_sendrecv
        ndone  = 0
	    nrecvd = 0
        do while(ndone < iireq_p)
            !>
            !>
            call mpi_testsome(iireq1,&
                              ireq_r,&
                              nrecvd,&
                              index_r,&
                              istat2,&
                              ierr)
            !>
            !>
            if(nrecvd > 0)then
                !>
                !>
                ndone = ndone + nrecvd
                nrecvddo: do nnn=1,nrecvd
                    !>
                    !>
                    ii_wait   =  keep2(index_r(nnn))
                    ibc       => mesh%bcs(ii_wait)
                    sblock    =  ibc%block_index
                    iblock    => mesh%blocks(sblock)
                    nd_src    = n2p(sblock)
                    nd_dest   = n2p(ibc%iblock_index)
                    if(myid == nd_src .and. myid /= nd_dest)then
                        ii_check_1 = keep1(index_r(nnn))
                        ii_point   = ii_check_1
                        dblock     = ibc%iblock_index
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
                        kstep = sign(1,je - ks)
                        !>
                        !>
                        istep_cut = sign(1,ie_cut -is_cut)
                        jstep_cut = sign(1,je_cut -js_cut)
                        kstep_cut = sign(1,ke_cut -ks_cut)
                        !>
                        !>
                        write(183,'("recv=",7i4)') ii_wait,is,ie,js,je,ks,ke
                        write(183,*) '*************'
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

                                    write(183,'("recv=",7i4," recv to ",i4," from=",i4)') ii_wait,id1,jd1,kd1,id2,jd2,kd2,nd_src,nd_dest
!                                    write(*,'("the axis for recv=",7i4)') dblock,id1,jd1,kd1,id2,jd2,kd2
                                    !> send the conservation variable of 1-1 block bc to the dest blocks
                                    !> the face of send include idim+1,idim+2,idim+3
                                    !> and i0-2,i0-1,i0
                                    !> the conservation variable at the cell centre
                                    !>
                                    iblock%cell(id1,jd1,kd1)%r =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id1,jd1,kd1)%u =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id1,jd1,kd1)%v =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id1,jd1,kd1)%w =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id1,jd1,kd1)%p =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id2,jd2,kd2)%r =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id2,jd2,kd2)%u =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id2,jd2,kd2)%v =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id2,jd2,kd2)%w =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    iblock%cell(id2,jd2,kd2)%p =  wk(ii_point)
                                    ii_point            = ii_point+1
                                    !>
                                    !>
                                    iblock%variable(id1,jd1,kd1)%res_1 = wk(ii_point)
                                    ii_point                    = ii_point+1
                                    !>
                                    iblock%variable(id1,jd1,kd1)%res_2 = wk(ii_point)
                                    ii_point                    = ii_point+1
                                    !>
                                    iblock%variable(id1,jd1,kd1)%res_3 = wk(ii_point)
                                    ii_point                    = ii_point+1
                                    !>
                                    iblock%variable(id1,jd1,kd1)%res_4 = wk(ii_point)
                                    ii_point                    = ii_point+1
                                    !>
                                    iblock%variable(id1,jd1,kd1)%res_5 = wk(ii_point)
                                    ii_point                    = ii_point+1
                                end do
                            end do
                        end do
                    endif
                end do nrecvddo
            endif
        enddo
        !>
        !>
#endif
        return
    end subroutine bc_1to1_mg

