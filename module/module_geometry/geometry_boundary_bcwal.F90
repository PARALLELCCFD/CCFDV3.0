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
    subroutine bc_wall(num_of_bc,nbl)
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use bc_module
        !>
        !>
		implicit none
        !>
        type(bc_types),pointer     :: ibc
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !>
        !>
        integer :: num_of_bc ,nbl
        integer :: is,ie,&
                   js,je,&
                   ks,ke,&
                   istep,&
                   jstep,&
                   kstep
        integer :: ip,jp,kp,&
                   i,j,k,&
                   id1,jd1,kd1,&
                   id2,jd2,kd2,&
                   idir
        !>
        !>
        real(kind = dprec) :: ubb,vbb,wbb,&
                              uu,vv,ww,uban,&
                              epson,twtype,tminnumber,&
                              agx,agy,agz,cq,m,&
                              pw,dpw,c2,xm2
        !>
        !>
        character(len = 20 ) :: iface
        !>
        !>
        mesh   => grids(imesh)
        ibc    => mesh%bcs(num_of_bc)
        iblock => mesh%blocks(nbl)
        !> boundary conditions for opposed symmetry (on solid wall)
        epson      = 1.e-15

        !> set twtype
        !> twtype > 0 ... fixed wall temperature twall/tinf = twtype
        !> twtype = 0 ... adiabatic wall
        !> twtype < 0 ... fixed wall temperature = stagnation temp
        twall  = 0.0
        cq     = 0.0
        twtype = twall/tinf
        !> set suction/blowing boundary conditions
        !>via mass-flow coefficient
        !>
        !>                                 (rho * u)
        !> mass-flow coefficient cq = ------------------
        !>                                (rho0 * u0)
        !> cq = 0...standard solid wall b.c, no flow thru wall
        !> cq < 0...suction (mass flow out of the zone)
        !> cq > 0...blowing (mass flow into the zone)
        !>
        tminnumber = 1.0e-30
        !>
        is = ibc%istart
        ie = ibc%iend
        js = ibc%jstart
        je = ibc%jend
        ks = ibc%kstart
        ke = ibc%kend
        if(is .eq. ie )then
            iface = 'i'
        else if(js .eq. je )then
            iface = 'j'
        else if(ks .eq. ke )then
            iface = 'k'
        else
            write(*,*) 'ccfdv3.0:: have some error at the wall boundary face,check the grid and the boundary !!!'
        end if
        !>
!        dir   = bc_norm(mm)
!        nmco  = bc_dir(mm)
        !>
        istep = sign(1,ie-is)
        jstep = sign(1,je-js)
        kstep = sign(1,ke-ks)
        !>
        !>
        !>
        !>
        if(nvisc .ne. 0 .and. ibc%bc_type .eq. 'wallinviscid')then
            write(*,'("error::at wall bondary conditions bc_wall subroutine ")')
            write(*,'("error::the boundary conditions is invicous wall,and must computing by euler not n-s")')
            stop
        end if
        !>
        !>
        !>
        !>
        if(nvisc .ne. 0 .and.  ibc%bc_type .eq. 'wallviscid') then

            !> viscous boundary condition
            !> no slip boundary condition
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
                        id2 = i
                        jd2 = j
                        kd2 = k
                        !>
                        !>
                        ip  = i
                        jp  = j
                        kp  = k
                        if (iface .eq. 'i') then
                            !>
                            !>
                            id1 = i+1*ibc%direction
                            id2 = i+2*ibc%direction
                            ip  = i+sign(1,3-i)
                            idir= sign(1,3-i)
                            !>
                            !>
                            agx = iblock%metric(id1,jd1,kd1)%fi
                            agy = iblock%metric(id1,jd1,kd1)%fj
                            agz = iblock%metric(id1,jd1,kd1)%fk
                            !>
                            !>
                        else if (iface .eq. 'j') then
                            !>
                            !>
                            jd1 = j+1*ibc%direction
                            jd2 = j+2*ibc%direction
                            jp  = j+sign(1,3-j)
                            idir= sign(1,3-j)
                            !>
                            !>
                            agx = iblock%metric(id1,jd1,kd1)%gi
                            agy = iblock%metric(id1,jd1,kd1)%gj
                            agz = iblock%metric(id1,jd1,kd1)%gk
                            !>
                            !>
                        else if (iface .eq. 'k') then
                            !>
                            !>
                            kd1 = k+1*ibc%direction
                            kd2 = k+2*ibc%direction
                            kp  = k+sign(1,3-k)
                            idir= sign(1,3-k)
                            !>
                            !>
                            agx = iblock%metric(id1,jd1,kd1)%hi
                            agy = iblock%metric(id1,jd1,kd1)%hj
                            agz = iblock%metric(id1,jd1,kd1)%hk
                            !>
                            !>
                        end if
                        !>
                        !>
                        !> surface velocities
                        uu = 0.0
                        vv = 0.0
                        ww = 0.0
                        !>pressure gradient
                        pw  = iblock%cell(i,j,k)%p
                        dpw = real(idir)*(iblock%cell(ip,jp,kp)%p - iblock%cell(i,j,k)%p)
!                        write(*,'(" bc info =",6i6)') is,ie,js,je,ks,ke
!                        write(*,'("nbl=",i6," bc_num=",i6," dp =(ip,jp,kp)-",3i5,"(i,j,k)=",4i5)') nbl,num_of_bc,ip,jp,kp,i,j,k,idir
                        pw  = pw - real(idir)*dpw/2.0
                        if(pw .le. 0.e0) pw= iblock%cell(i,j,k)%p
                        !>
                        c2 = gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r
                        !>
                        !>
                        if(twtype .gt. 0.0)then
                            c2 = twtype
                        else if(twtype .lt. 0.0)then
                            c2 = 1.0 + 0.5e0*gm1*xmach*xmach
                        else !twtype == 0.0
                            xm2 = iblock%cell(i,j,k)%u**2 + &
                                  iblock%cell(i,j,k)%v**2 + &
                                  iblock%cell(i,j,k)%w**2
                            xm2 = xm2 /c2
                            c2  = c2*(1.e0 + 0.5*gm1*xm2)
                        end if
                        uu = uu + xmach*cq*agx*c2/(gamma*pw)
                        vv = vv + xmach*cq*agy*c2/(gamma*pw)
                        ww = ww + xmach*cq*agz*c2/(gamma*pw)
                        !>
                        iblock%cell(id1,jd1,kd1)%r  = gamma*pw/c2
                        iblock%cell(id1,jd1,kd1)%p  = pw
                        !>
                        iblock%cell(id1,jd1,kd1)%u = uu
                        iblock%cell(id1,jd1,kd1)%v = vv
                        iblock%cell(id1,jd1,kd1)%w = ww
!                        write(111,'(3i5,5e24.16)') id1,jd1,kd1,iblock%cell(id1,jd1,kd1)%r,iblock%cell(id1,jd1,kd1)%p,agx,agy,agz
                        !>
                        !>
                        !>
                        iblock%cell(id1,jd1,kd1)%wall_blank = 2.0
                        !>
                        !> need to store the viscous of turbulent to computing the viscous coefficient
                        !>
                        if(ivisc_i .ge. 2 .or. ivisc_j .gt. 2 .or. ivisc_k .gt. 2)then
                            iblock%turbulent(id1,jd1,kd1)%viscous = 0.0
                        end if
                        !>
                        !>
                        !> only need to do advanced model turbulence B.C.s on finest grid
                        if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
                            iblock%turbulent(id1,jd1,kd1)%tur_save = -iblock%turbulent(i,j,k)%tur_save
!                            write(*,'("i=",i4," idim=",i4," j=",i4," jdim=",i4," k=",i4," kdim=",i4," boundary axis=",3i4)') i,iblock%idim,j,iblock%jdim,k,iblock%kdim,id1,jd1,kd1
!                            write(*,'("axis=",3i4," viscous=",2e24.16)')id1,jd1,kd1,iblock%turbulent(id1,jd1,kd1)%viscous,iblock%turbulent(id1,jd1,kd1)%tur_save
                        end if
                        !>
                        !>
                        !>
                        if(iblock%cell(id1,jd1,kd1)%r <1.e-10 .or. iblock%cell(id1,jd1,kd1)%p < 1.e-10)then
                           open(unit=89,file='ccfd.bc_error',status='unknown')
                           write(89,'("iter=",i6," at the bc_wall_11 block=",i3,"  rho is too small r1(",i3,",",i3,",",i3,")=",&
                           e20.12," p1(id1,jd1,kd1)=",e20.12)') icyc,ibc%block_index,id1,jd1,kd1,iblock%cell(id1,jd1,kd1)%r,iblock%cell(id1,jd1,kd1)%p
                           write(89,'("iter=",i6," at the bc_wall_11 block=",i3,"  rho is too small r0(",i3,",",i3,",",i3,")=",&
                           e20.12," p0(i,j,k)=",e20.12)') icyc,ibc%block_index,i,j,k,iblock%cell(i,j,k)%r,iblock%cell(i,j,k)%p
                           write(89,*) '-------------------------------------------------------------------------------------'
                           !close(69)
                        endif
                        !>
                        iblock%cell(id2,jd2,kd2)%r = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%r - iblock%cell(id1,jd1,kd1)%r)
                        iblock%cell(id2,jd2,kd2)%p = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%p - iblock%cell(id1,jd1,kd1)%p)
                        !>
                        !>
                        !>
                        iblock%cell(id2,jd2,kd2)%u = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%u - iblock%cell(id1,jd1,kd1)%u)
                        iblock%cell(id2,jd2,kd2)%v = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%v - iblock%cell(id1,jd1,kd1)%v)
                        iblock%cell(id2,jd2,kd2)%w = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%w - iblock%cell(id1,jd1,kd1)%w)
!                        write(500+icyc,'(4i3,10e18.12)') num_of_bc,i,j,k,iblock%cell(id1,jd1,kd1)%r,iblock%cell(id1,jd1,kd1)%u,&
!                        iblock%cell(id1,jd1,kd1)%v,iblock%cell(id1,jd1,kd1)%w,iblock%cell(id1,jd1,kd1)%p,&
!                        iblock%cell(id2,jd2,kd2)%r,iblock%cell(id2,jd2,kd2)%u,iblock%cell(id2,jd2,kd2)%v,&
!                        iblock%cell(id2,jd2,kd2)%w,iblock%cell(id2,jd2,kd2)%p
                        !>
                        !>
                        iblock%cell(id2,jd2,kd2)%wall_blank = 1.0
                        !>
                        !> need to store the viscous of turbulent to computing the viscous coefficient
                        !>
                        if(nvisc .ge. 2)then
                            iblock%turbulent(id2,jd2,kd2)%viscous = 0.0
                        end if

                    end do
                end do
            end do

        else
            !>
            !> used the euler computing the viscous wall bc
            !>
            !>
            if(nvisc .eq. 0 .and. ibc%bc_type .eq. 'wallviscid')then
               write(*,'("warnning::the wall boundary condition is viscous bc ")')
               write(*,'("warnning::and you will used the euler computing the mesh")')
               write(*,'("warnning::please check the condition or the mesh!")')
            end if
            !> inviscid wall boundary condition
            !> slip b.c
            do k = ks,ke,kstep
                do j = js,je,jstep
                    do i = is,ie,istep
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

                        if (iface .eq. 'i') then
                            id1 = i + 1*ibc%direction
                            id2 = i + 2*ibc%direction
                            !>
                            !>
                            agx = iblock%metric(id1,jd1,kd1)%fi
                            agy = iblock%metric(id1,jd1,kd1)%fj
                            agz = iblock%metric(id1,jd1,kd1)%fk
                            !>
                            !>

                        else if(iface .eq. 'j') then

                            jd1 = j + 1*ibc%direction
                            jd2 = j + 2*ibc%direction
                            !>
                            !>
                            agx = iblock%metric(id1,jd1,kd1)%gi
                            agy = iblock%metric(id1,jd1,kd1)%gj
                            agz = iblock%metric(id1,jd1,kd1)%gk
                            !>
                            !>
                        else  if(iface .eq. 'k') then

                            kd1 = k + 1*ibc%direction
                            kd2 = k + 2*ibc%direction
                            !>
                            !>
                            agx = iblock%metric(id1,jd1,kd1)%hi
                            agy = iblock%metric(id1,jd1,kd1)%hj
                            agz = iblock%metric(id1,jd1,kd1)%hk
                            !>
                            !>
                        end if
!                        if(iblock%metric(id1,jd1,kd1)%hi .ne. iblock%metric(i,j,k)%hi) write(*,'(" at axis=",6i4," different the metric=",2e24.16)') i,j,k,id1,jd1,kd1,iblock%metric(i,j,k)%hi,iblock%metric(id1,jd1,kd1)%hi
!                        write(*,'(" bc info =",6i6)') is,ie,js,je,ks,ke
!                        write(*,'("nbl=",i6," bc_num=",i6," dp =(id1,jd1,kd1)",3i5," (i,j,k)=",3i5)') nbl,num_of_bc,id1,jd1,kd1,i,j,k
                        !>
                        !>
                        uu = iblock%cell(i,j,k)%u
                        vv = iblock%cell(i,j,k)%v
                        ww = iblock%cell(i,j,k)%w
                        !>
                        !>
                        uban = uu*agx + vv*agy + ww*agz

                        ubb = uu - uban*agx
                        vbb = vv - uban*agy
                        wbb = ww - uban*agz
                        !>
                        !>
                        iblock%cell(id1,jd1,kd1)%r = iblock%cell(i,j,k)%r
                        iblock%cell(id1,jd1,kd1)%p = iblock%cell(i,j,k)%p
                        !>
                        !>
                        !>
                        iblock%cell(id1,jd1,kd1)%u  =  ubb
                        iblock%cell(id1,jd1,kd1)%v  =  vbb
                        iblock%cell(id1,jd1,kd1)%w  =  wbb
                        !>
                        iblock%cell(id1,jd1,kd1)%wall_blank = 2.0
                        !>
                        !> need to store the viscous of turbulent to computing the viscous coefficient
                        !>
                        if(nvisc .ge. 2)then
                            iblock%turbulent(id1,jd1,kd1)%viscous = iblock%turbulent(i,j,k)%viscous
                        end if
                        !>
                        !>
                        !> only need to do advanced model turbulence B.C.s on finest grid
                        if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
                            iblock%turbulent(id1,jd1,kd1)%tur_save = iblock%turbulent(i,j,k)%tur_save
                        end if
                        !>
                        !>
                        !>
                        iblock%cell(id2,jd2,kd2)%r  = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%r - iblock%cell(id1,jd1,kd1)%r )
                        iblock%cell(id2,jd2,kd2)%p  = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%p - iblock%cell(id1,jd1,kd1)%p )
                        !>
                        !>
                        iblock%cell(id2,jd2,kd2)%u = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%u -  iblock%cell(id1,jd1,kd1)%u )
                        iblock%cell(id2,jd2,kd2)%v = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%v -  iblock%cell(id1,jd1,kd1)%v )
                        iblock%cell(id2,jd2,kd2)%w = 2.0*ibc%norm_index*(iblock%cell(i,j,k)%w -  iblock%cell(id1,jd1,kd1)%w )
                        !>
                        !>
                        iblock%cell(id2,jd2,kd2)%wall_blank = 1.0
                        !>
                        !>
                        !>
                        !> need to store the viscous of turbulent to computing the viscous coefficient
                        !>
                        if(nvisc .ge. 2)then
                            iblock%turbulent(id2,jd2,kd2)%viscous = iblock%turbulent(i,j,k)%viscous
                        end if
                        !?
                        !>
                        !>
                    enddo
                enddo
            enddo

        end if


        return
    end subroutine bc_wall


