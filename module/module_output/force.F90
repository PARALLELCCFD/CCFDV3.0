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
    !*************************************************************************!
    !                                                                         !
    !  subroutine:   force.F90                                                !
    !                                                                         !
    !  programmer:   liuxz                                                    !
    !                                                                         !
    !  update:       2019-1-7 23:04:10                                        !
    !*************************************************************************!
!----------------------------------------------------------------------------------------------
    subroutine force
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use bc_module
        use cell_module
        use coordinate_module
        use metric_module
        use turbulent_module
        use nodes_mg
        use nodes_var
        use nodes_paras
        !>
        !>
		implicit none
#if defined PMPI
        include "mpif.h"
#endif
        type(bc_types),pointer ::ibc
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !>
        !>
        real(kind=dprec) cpc,const,&
                         cosa,sina,cosb,sinb
        real(kind=dprec) fxyz(6,3),fxyz_out(6,3)
        real(kind=dprec) clw,cdw,&
                         xc,yc,zc,xmc,ymc,zmc,&
                         fx,fy,fz,cx,cy,cz,&
                         pfx,pfy,pfz,pcx,pcy,pcz,&
                         dfx,dfy,dfz,dcx,dcy,dcz,&
                         dpfx,dpfy,dpfz,dpcx,dpcy,dpcz,&
                         mex,mey,mez,mee,&
                         viscous,vnorm,dpf,&
                         tau
        integer :: i,j,k,ii,ijk,&
                   idim,jdim,kdim,&
                   is,ie,js,je,ks,ke,&
                   istep,jstep,kstep,&
                   nbl,icyc1,iseq,ierr
        !>
        mesh => grids(imesh)
        !>
        !>
#if defined PMPI
		if(myid .eq. myhost)then
#endif
            !>
            !>
            icyc1 = icyc
            !>
            do iseq=1,imeshseque-1
                icyc1 = icyc1 + ncycle(iseq)
            end do
#if defined PMPI
		end if
#endif
        !>
        !>
        do ii =1,6
            fxyz(ii,1) = 0.0
            fxyz(ii,2) = 0.0
            fxyz(ii,3) = 0.0
            fxyz_out(ii,1) = 0.0
            fxyz_out(ii,2) = 0.0
            fxyz_out(ii,3) = 0.0
        end do
        !>
        !>
        !>
        xmc = x_ref
        ymc = y_ref
        zmc = z_ref
        !>
        !>
        !>
        cpc   = 2.0/(gamma*xmach*xmach)
        const = 4.0/(reue*xmach)
        !>
        !>
        cosa  = cos(alpha)
        sina  = sin(alpha)
        cosb  = cos(beta)
        sinb  = sin(beta)
        !>
        !>
        !>
        do ii =1,num_bc
            !>
            !>
            ibc => mesh%bcs(ii)
            !>
            !>
            if(ibc%bc_type.ne.'wallviscid') cycle
            !>
            !>
            nbl = ibc%block_index
            !>
            !>
#if defined PMPI
            if(n2p(nbl) .eq. myid)then
#endif

				if((level_mg(nbl) .eq. global_level .and. meshseqflags .eq. 0) &
				&.or. (level_mg(nbl) .eq. imeshseque .and. meshseqflags .ne. 0))then
                    !>
                    !>
					iblock  =>  mesh%blocks(nbl)
					!>
					idim = iblock%idim
					jdim = iblock%jdim
					kdim = iblock%kdim
					!>
					!>
					is  = ibc%istart
					ie  = ibc%iend
					js  = ibc%jstart
					je  = ibc%jend
					ks  = ibc%kstart
					ke  = ibc%kend
					!>
					!>
					!>
					istep = 1
					jstep = 1
					kstep = 1
					!>
					!>
					if(is .gt. ie)then
						istep = -1
					end if
					if(js .gt. je)then
						jstep = -1
					end if
					if(ks .gt. ke)then
						kstep = -1
					end if
					!>
					!>
					fx   = 0.0
					fy   = 0.0
					fz   = 0.0
					cx   = 0.0
					cy   = 0.0
					cz   = 0.0
					pfx  = 0.0
					pfy  = 0.0
					pfz  = 0.0
					pcx  = 0.0
					pcy  = 0.0
					pcz  = 0.0
					!>
					!>
					!******************************************************************
					!     forces on i=constant surfaces
					!******************************************************************
                    if(is .eq. ie) then
                        !>
                        !>
                        i  = is
                        if(is .le. 2)       ijk = is -1
                        if(ie .ge. idim)    ijk = is + 1
                        !>
                        !>
                        do k  = ks,ke,kstep
                            do j  = js,je,jstep
                                !>
                                !>
                                mex = ibc%norm_index *iblock%metric(i,j,k)%fi
								mey = ibc%norm_index *iblock%metric(i,j,k)%fj
								mez = ibc%norm_index *iblock%metric(i,j,k)%fk
								mee = iblock%metric(i,j,k)%ff
								!>
								!>
								xc= (iblock%coordinate(i,j  ,k  )%x&
								   &+iblock%coordinate(i,j+1,k  )%x&
								   &+iblock%coordinate(i,j  ,k+1)%x&
								   &+iblock%coordinate(i,j+1,k+1)%x)*0.25
                                !>
                                !>
                                yc=  (iblock%coordinate(i,j  ,k  )%y&
                                    &+iblock%coordinate(i,j+1,k  )%y&
                                    &+iblock%coordinate(i,j  ,k+1)%y&
                                    &+iblock%coordinate(i,j+1,k+1)%y)*0.25
                                !>
                                !>
                                zc=  (iblock%coordinate(i,j  ,k  )%z&
                                    &+iblock%coordinate(i,j+1,k  )%z&
                                    &+iblock%coordinate(i,j  ,k+1)%z&
                                    &+iblock%coordinate(i,j+1,k+1)%z)*0.25
                                !>
                                !>
                                dpf  = -(gamma*iblock%cell(ijk,j,k)%p-1.0) *cpc *mee
                                dpfx = dpf *mex
                                dpfy = dpf *mey
                                dpfz = dpf *mez
                                !>
                                !>
                                pfx   = pfx  +dpfx
                                pfy   = pfy  +dpfy
                                pfz   = pfz  +dpfz
                                !>
                                !>
                                pcx  =  pcx  +dpfz*(yc-ymc) -dpfy*(zc-zmc)
                                pcy  =  pcy  -dpfz*(xc-xmc) +dpfx*(zc-zmc)
                                pcz  =  pcz  +dpfy*(xc-xmc) -dpfx*(yc-ymc)
                                !>
                                !>
                                if((nvisc.eq.0).or.(ibc%bc_type.ne.'wallviscid')) cycle
                                !>
                                !>
                                vnorm =  mex*iblock%cell(i,j,k)%u + mey*iblock%cell(i,j,k)%v + mez*iblock%cell(i,j,k)%w
                                !>
                                !>
                                viscous = iblock%cell(ijk,j,k)%viscous
                                !>
                                !>
                                if(nvisc.ge.2) then
                                    viscous = viscous + iblock%turbulent(ijk,j,k)%viscous
                                end if
                                !>
                                !>
                                tau = viscous *const *iblock%metric(i,j,k)%jacobian *mee**2
                                !>
                                dfx  = tau*(iblock%cell(i,j,k)%u - vnorm*mex)
                                dfy  = tau*(iblock%cell(i,j,k)%v - vnorm*mey)
                                dfz  = tau*(iblock%cell(i,j,k)%w - vnorm*mez)
                                !>
                                !>
                                fx  = fx +dfx
                                fy  = fy +dfy
                                fz  = fz +dfz
                                cx  = cx +dfz*(yc-ymc) -dfy*(zc-zmc)
                                cy  = cy -dfz*(xc-xmc) +dfx*(zc-zmc)
                                cz  = cz +dfy*(xc-xmc) -dfx*(yc-ymc)
                            end do
                        end do
                    end if
                    !******************************************************************
                    !     forces on j=constant surfaces
                    !******************************************************************
                    if(js .eq. je) then
                        !>
                        !>
                        j  = js
                        if(js .le. 2)       ijk = js -1
                        if(je .ge. jdim)    ijk = je + 1
                        !>
                        !>
                        do k  = ks,ke,kstep
                            do i  = is,ie,istep
                                mex = ibc%norm_index *iblock%metric(i,j,k)%gi
                                mey = ibc%norm_index *iblock%metric(i,j,k)%gj
                                mez = ibc%norm_index *iblock%metric(i,j,k)%gk
                                mee = iblock%metric(i,j,k)%gg
                                !>
                                !>
                                xc= (iblock%coordinate(i,  j,k  )%x&
                                   &+iblock%coordinate(i+1,j,k  )%x&
                                   &+iblock%coordinate(i,  j,k+1)%x&
                                   &+iblock%coordinate(i+1,j,k+1)%x)*0.25
                                !>
                                yc= (iblock%coordinate(i,  j,k  )%y&
                                   &+iblock%coordinate(i+1,j,k  )%y&
                                   &+iblock%coordinate(i,  j,k+1)%y&
                                   &+iblock%coordinate(i+1,j,k+1)%y)*0.25
                                !>
                                zc= (iblock%coordinate(i,  j,k  )%z&
                                   &+iblock%coordinate(i+1,j,k  )%z&
                                   &+iblock%coordinate(i,  j,k+1)%z&
                                   &+iblock%coordinate(i+1,j,k+1)%z)*0.25
                                !>
                                !>
                                !>
                                dpf  = -(gamma*iblock%cell(i,ijk,k)%p-1.0) *cpc * mee
                                dpfx = dpf *mex
                                dpfy = dpf *mey
                                dpfz = dpf *mez
                                !>
                                !>
                                pfx   = pfx  +dpfx
                                pfy   = pfy  +dpfy
                                pfz   = pfz  +dpfz
                                pcx  =  pcx  +dpfz*(yc-ymc) -dpfy*(zc-zmc)
                                pcy  =  pcy  -dpfz*(xc-xmc) +dpfx*(zc-zmc)
                                pcz  =  pcz  +dpfy*(xc-xmc) -dpfx*(yc-ymc)
                                !>
                                !>
                                !>
                                if((nvisc.eq.0).or.(ibc%bc_type.ne.'wallviscid')) cycle
                                !>
                                !>
                                vnorm = mex*iblock%cell(i,j,k)%u + mey*iblock%cell(i,j,k)%v + mez*iblock%cell(i,j,k)%w
                                !>
                                !>
                                viscous = iblock%cell(i,ijk,k)%viscous
                                !>
                                !>
                                if(nvisc.ge.2) then
                                    viscous = viscous + iblock%turbulent(i,ijk,k)%viscous
                                end if
                                !>
                                !>
                                tau = viscous *const *iblock%metric(i,j,k)%jacobian *mee**2
                                !>
                                !>
                                dfx  = tau*(iblock%cell(i,j,k)%u - vnorm*mex)
                                dfy  = tau*(iblock%cell(i,j,k)%v - vnorm*mey)
                                dfz  = tau*(iblock%cell(i,j,k)%w - vnorm*mez)
!                                write(1028,'(3i4,5e24.16)') i,j,k,vnorm,tau,iblock%cell(i,j,k)%u,iblock%cell(i,j,k)%v,iblock%cell(i,j,k)%w
                                !>
                                !>
                                fx  = fx +dfx
                                fy  = fy +dfy
                                fz  = fz +dfz
                                !>
                                !>
                                cx  = cx +dfz*(yc-ymc) -dfy*(zc-zmc)
                                cy  = cy -dfz*(xc-xmc) +dfx*(zc-zmc)
                                cz  = cz +dfy*(xc-xmc) -dfx*(yc-ymc)
                                !>
                                !>
                            end do
                        end do
                    end if
                    !******************************************************************
                    !     forces on k=constant surfaces
                    !******************************************************************
                    if(ks .eq. ke) then
                        !>
                        !>
                        k  = ks
                        if(ks .le. 2)       ijk = ks -1
                        if(ke .ge. kdim)    ijk = ks + 1
                        !>
                        !>
                        do j  = js,je,jstep
                            do i  = is,ie,istep
                                !>
                                !>
                                mex = ibc%norm_index *iblock%metric(i,j,k)%hi
                                mey = ibc%norm_index *iblock%metric(i,j,k)%hj
                                mez = ibc%norm_index *iblock%metric(i,j,k)%hk
                                mee = iblock%metric(i,j,k)%hh
                                !>
                                !>
                                xc= (iblock%coordinate(i  ,j  ,k)%x&
                                   &+iblock%coordinate(i+1,j  ,k)%x&
                                   &+iblock%coordinate(i,  j+1,k)%x&
                                   &+iblock%coordinate(i+1,j+1,k)%x)*0.25
                                !>
                                !>
                                yc= (iblock%coordinate(i  ,j  ,k)%y&
                                   &+iblock%coordinate(i+1,j  ,k)%y&
                                   &+iblock%coordinate(i,  j+1,k)%y&
                                   &+iblock%coordinate(i+1,j+1,k)%y)*0.25
                                !>
                                !>
                                zc= (iblock%coordinate(i  ,j  ,k)%z&
                                   &+iblock%coordinate(i+1,j  ,k)%z&
                                   &+iblock%coordinate(i,  j+1,k)%z&
                                   &+iblock%coordinate(i+1,j+1,k)%z)*0.25
                                !>
                                !>
                                dpf  = -(gamma*iblock%cell(i,j,ijk)%p-1.0) *cpc * mee
                                dpfx = dpf *mex
                                dpfy = dpf *mey
                                dpfz = dpf *mez
                                !>
                                !>
                                pfx   = pfx  +dpfx
                                pfy   = pfy  +dpfy
                                pfz   = pfz  +dpfz
                                !>
                                !>
                                pcx  =  pcx  +dpfz*(yc-ymc) -dpfy*(zc-zmc)
                                pcy  =  pcy  -dpfz*(xc-xmc) +dpfx*(zc-zmc)
                                pcz  =  pcz  +dpfy*(xc-xmc) -dpfx*(yc-ymc)
                                !>
                                !>
                                if((nvisc.eq.0).or.(ibc%bc_type.ne.'wallviscid')) cycle
                                !>
                                !>
                                vnorm = mex *iblock%cell(i,j,k)%u + mey*iblock%cell(i,j,k)%v + mez*iblock%cell(i,j,k)%w
                                !>
                                !>
                                viscous = iblock%cell(i,j,ijk)%viscous
                                if(nvisc.ge.2) then
                                    viscous = viscous + iblock%turbulent(i,j,ijk)%viscous
                                end if
                                !>
                                !>
                                tau = viscous *const *iblock%metric(i,j,k)%jacobian *mee**2
                                !>
                                !>
                                dfx  = tau*(iblock%cell(i,j,k)%u - vnorm*mex)
                                dfy  = tau*(iblock%cell(i,j,k)%v - vnorm*mey)
                                dfz  = tau*(iblock%cell(i,j,k)%w - vnorm*mez)
                                !>
                                !>
                                fx  = fx +dfx
                                fy  = fy +dfy
                                fz  = fz +dfz
                                !>
                                !>
                                cx  = cx +dfz*(yc-ymc) -dfy*(zc-zmc)
                                cy  = cy -dfz*(xc-xmc) +dfx*(zc-zmc)
                                cz  = cz +dfy*(xc-xmc) -dfx*(yc-ymc)
                                !>
                                !>
                            end do
                        end do
                    end if
                    !>
                    !>
                    fxyz(1,1) = fxyz(1,1) + pfx
                    fxyz(2,1) = fxyz(2,1) + pfy
                    fxyz(3,1) = fxyz(3,1) + pfz
                    fxyz(4,1) = fxyz(4,1) + pcx
                    fxyz(5,1) = fxyz(5,1) + pcy
                    fxyz(6,1) = fxyz(6,1) + pcz

                    fxyz(1,2) = fxyz(1,2) + fx
                    fxyz(2,2) = fxyz(2,2) + fy
                    fxyz(3,2) = fxyz(3,2) + fz
                    fxyz(4,2) = fxyz(4,2) + cx
                    fxyz(5,2) = fxyz(5,2) + cy
                    fxyz(6,2) = fxyz(6,2) + cz

                    fxyz(1,3) = fxyz(1,3) + pfx + fx
                    fxyz(2,3) = fxyz(2,3) + pfy + fy
                    fxyz(3,3) = fxyz(3,3) + pfz + fz
                    fxyz(4,3) = fxyz(4,3) + pcx + cx
                    fxyz(5,3) = fxyz(5,3) + pcy + cy
                    fxyz(6,3) = fxyz(6,3) + pcz + cz
                !>
                !>
                !>
                endif
                !>
                !>
#if defined PMPI
            endif
#endif
            !>
            !>
        end do !do ii =1,num_bc
        !>
        !>
        !>
#if defined PMPI
        !>
        !>
        !> reduce all the values
        call mpi_allreduce(fxyz,fxyz_out,18,mpi_double_precision,mpi_sum,mycomm,ierr)
        !>
        !>
        do ii =1,6
            fxyz(ii,1) = fxyz_out(ii,1)
            fxyz(ii,2) = fxyz_out(ii,2)
            fxyz(ii,3) = fxyz_out(ii,3)
        end do
        !>
        !>
        if(myid .eq. myhost)then
#endif

            cd =    fxyz(1,3)*cosa*cosb + fxyz(3,3)*sina*cosb - fxyz(2,3)*sinb
            cl =   -fxyz(1,3)*sina + fxyz(3,3)*cosa

#if defined PMPI
        end if
#endif
        !>
        !>
        !>
        return
    end subroutine force

