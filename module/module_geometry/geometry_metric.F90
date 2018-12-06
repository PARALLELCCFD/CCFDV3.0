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
    subroutine metric()
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
        use coordinate_module
        use bc_module
        use metric_module
        !> the mpi module
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
        !>
        !>
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !>
        real(kind = dprec):: x1,y1,z1
        real(kind = dprec):: x2,y2,z2
        real(kind = dprec):: x3,y3,z3
        real(kind = dprec):: x4,y4,z4
        real(kind = dprec):: x5,y5,z5
        real(kind = dprec):: x6,y6,z6
        real(kind = dprec):: x7,y7,z7
        real(kind = dprec):: x8,y8,z8
        !>
        !>
        real(kind = dprec):: t1,t2,t3
        real(kind = dprec):: t4,t5,t6
        real(kind = dprec):: kap
        !> the volume of six
        real(kind = dprec):: v1,&
                             v2,&
                             v3,&
                             v4,&
                             v5,&
                             v6
        !>
        real(kind = dprec),dimension(:),allocatable :: xfc
        real(kind = dprec),dimension(:),allocatable :: yfc
        real(kind = dprec),dimension(:),allocatable :: zfc
!        real*8,dimension(:),allocatable :: xfc
!        real*8,dimension(:),allocatable :: yfc
!        real*8,dimension(:),allocatable :: zfc
        real(kind = dprec),dimension(:,:,:),allocatable :: volume
        !>
        real(kind = dprec):: dxfc
        real(kind = dprec):: dyfc
        real(kind = dprec):: dzfc
        !>
        !>
        real(kind = dprec):: vx
        real(kind = dprec):: vy
        real(kind = dprec):: vz
        !>
        !>
        !>
        integer :: ii,jj,kk,ilevel
        integer :: nbl,i,j,k
        integer :: idim,jdim,kdim
        !>
        !>
        allocate(xfc(ijkmax))
        allocate(yfc(ijkmax))
        allocate(zfc(ijkmax))
        !> computation of cell metric
        !>
        !>
        mesh => grids(imesh)
        kap = 1.e0/24.e0


        do nbl =1,nblocks,global_level
#if defined PMPI
            if(myid .eq. n2p(nbl) .or. myid .eq. myhost)then
#endif
                !>
                !>
                iblock => mesh%blocks(nbl)
                !>
                idim = iblock%idim
                jdim = iblock%jdim
                kdim = iblock%kdim
                !> the corser mesh's volume restict by finer mesh
                !> so the corser mesh's volume must not computing by metic
                !>


!                if(nbl==1)then
!                open(unit=110,file='coordinate.dat',status='unknown')
!                do k=2,kdim+1
!                do j=2,jdim+1
!                do i=2,idim+1
!                write(110,'("axis=",3i4,"  and coordinate=",3e24.16)')i-1,j-1,k-1,&
!                iblock%coordinate(i,j,k)%x,iblock%coordinate(i,j,k)%y,iblock%coordinate(i,j,k)%z
!                end do
!                end do
!                end do
!                end if


                do k  =2,kdim
                    do j  =2,jdim
                        do i  =2,idim
                            !>
                            !> the (i,j,k) coordinate
                            x1  = iblock%coordinate(i,j,k)%x
                            y1  = iblock%coordinate(i,j,k)%y
                            z1  = iblock%coordinate(i,j,k)%z
                            !> the (i,j+1,k) coordinate
                            x2  = iblock%coordinate(i,j+1,k)%x
                            y2  = iblock%coordinate(i,j+1,k)%y
                            z2  = iblock%coordinate(i,j+1,k)%z
                            !> the (i,j,k+1) coordinate
                            x3  = iblock%coordinate(i,j,k+1)%x
                            y3  = iblock%coordinate(i,j,k+1)%y
                            z3  = iblock%coordinate(i,j,k+1)%z
                            !> the (i+1,j,k) coordinate
                            x4  = iblock%coordinate(i+1,j,k)%x
                            y4  = iblock%coordinate(i+1,j,k)%y
                            z4  = iblock%coordinate(i+1,j,k)%z
                            !> the (i,j+1,k+1) coordinate
                            x5  = iblock%coordinate(i,j+1,k+1)%x
                            y5  = iblock%coordinate(i,j+1,k+1)%y
                            z5  = iblock%coordinate(i,j+1,k+1)%z
                            !> the (i+1,j,k+1) coordinate
                            x6  = iblock%coordinate(i+1,j,k+1)%x
                            y6  = iblock%coordinate(i+1,j,k+1)%y
                            z6  = iblock%coordinate(i+1,j,k+1)%z
                            !> the (i+1,j+1,k) coordinate
                            x7  = iblock%coordinate(i+1,j+1,k)%x
                            y7  = iblock%coordinate(i+1,j+1,k)%y
                            z7  = iblock%coordinate(i+1,j+1,k)%z
                            !> the (i+1,j+1,k+1) coordinate
                            x8  = iblock%coordinate(i+1,j+1,k+1)%x
                            y8  = iblock%coordinate(i+1,j+1,k+1)%y
                            z8  = iblock%coordinate(i+1,j+1,k+1)%z
                            !>
                            !>
                            !>
                            v1  =   x1*(y2*(z5-z8)-y5*(z2-z8)+y8*(z2-z5))-&
                                    x2*(y1*(z5-z8)-y5*(z1-z8)+y8*(z1-z5))+&
                                    x5*(y1*(z2-z8)-y2*(z1-z8)+y8*(z1-z2))-&
                                    x8*(y1*(z2-z5)-y2*(z1-z5)+y5*(z1-z2))

                            v2  =   x1*(y2*(z8-z7)-y8*(z2-z7)+y7*(z2-z8))-&
                                    x2*(y1*(z8-z7)-y8*(z1-z7)+y7*(z1-z8))+&
                                    x8*(y1*(z2-z7)-y2*(z1-z7)+y7*(z1-z2))-&
                                    x7*(y1*(z2-z8)-y2*(z1-z8)+y8*(z1-z2))

                            v3  =   x1*(y3*(z8-z5)-y8*(z3-z5)+y5*(z3-z8))-&
                                    x3*(y1*(z8-z5)-y8*(z1-z5)+y5*(z1-z8))+&
                                    x8*(y1*(z3-z5)-y3*(z1-z5)+y5*(z1-z3))-&
                                    x5*(y1*(z3-z8)-y3*(z1-z8)+y8*(z1-z3))

                            v4  =   x1*(y4*(z7-z8)-y7*(z4-z8)+y8*(z4-z7))-&
                                    x4*(y1*(z7-z8)-y7*(z1-z8)+y8*(z1-z7))+&
                                    x7*(y1*(z4-z8)-y4*(z1-z8)+y8*(z1-z4))-&
                                    x8*(y1*(z4-z7)-y4*(z1-z7)+y7*(z1-z4))

                            v5  =   x1*(y4*(z8-z6)-y8*(z4-z6)+y6*(z4-z8))-&
                                    x4*(y1*(z8-z6)-y8*(z1-z6)+y6*(z1-z8))+&
                                    x8*(y1*(z4-z6)-y4*(z1-z6)+y6*(z1-z4))-&
                                    x6*(y1*(z4-z8)-y4*(z1-z8)+y8*(z1-z4))

                            v6  =   x1*(y3*(z6-z8)-y6*(z3-z8)+y8*(z3-z6))-&
                                    x3*(y1*(z6-z8)-y6*(z1-z8)+y8*(z1-z6))+&
                                    x6*(y1*(z3-z8)-y3*(z1-z8)+y8*(z1-z3))-&
                                    x8*(y1*(z3-z6)-y3*(z1-z6)+y6*(z1-z3))

                            iblock%metric(i,j,k)%volume =(abs(v1) + &
                                                          abs(v2) + &
                                                          abs(v3) + &
                                                          abs(v4) + &
                                                          abs(v5) + &
                                                          abs(v6))/6.0

                            if(iblock%metric(i,j,k)%volume .le. 0.0) then
                                open(unit=71,file='ccfd.error',status='unknown')
                                write(71,'("CCFDV3.0::error:in block=",i4," i=",i4," j=",i4," k=",i4,"  the volume is zero, and the volume=",e24.12)') nbl,i,j,k,iblock%metric(i,j,k)%volume
                                write(71,'("CCFDV3.0::6-v=",6e24.12)') v1,v2,v3,v4,v5,v6
                                stop
                            end if
                            !>
                            !>
                            !> the jacobian value is 1/volume
                            iblock%metric(i,j,k)%jacobian= 1.0 / iblock%metric(i,j,k)%volume
!                            write(1221,'("axis=",3i4," metric=",e24.16)')i,j,k,iblock%metric(i,j,k)%volume
!                            write(113,'("iblock=",i4," the axis=",3i4," volume=",e20.12," the jacobian=",e20.12)') nbl,i,j,k,iblock%metric(i,j,k)%volume,iblock%metric(i,j,k)%jacobian
                            !>
                            !>
                        end do
                    end do
                end do
                !>
!                iblock => null
#if defined PMPI
            end if
#endif
        end do
        !>
        !>
        !>
        !>
        !>  computing the surface vector and the volume , jocabian
        !>  |fi,fj,fk|
        !>  |gi,gj,gk|
        !>  |hi,hj,hk|
        !>
        !>
        do nbl =1,nblocks
#if defined PMPI
            if(myid .eq. n2p(nbl) .or. myid .eq. myhost)then
#endif
                iblock => mesh%blocks(nbl)
                idim = iblock%idim
                jdim = iblock%jdim
                kdim = iblock%kdim
                allocate(volume(idim+1,jdim+1,kdim+1))
                !> volume and area calculation
                !> x-direciion
                do j  = 2,jdim
                    do k  = 2,kdim
                        do i  = 2,idim + 1
                            !>
                            !>
                            t4 =iblock%coordinate(i,j+1,k+1)%x -iblock%coordinate(i,j,k)%x
                            !>
                            !>
                            t5 =iblock%coordinate(i,j+1,k+1)%y -iblock%coordinate(i,j,k)%y
                            !>
                            !>
                            t6 =iblock%coordinate(i,j+1,k+1)%z -iblock%coordinate(i,j,k)%z
                            !>
                            !>
                            !>
                            t1 =-iblock%coordinate(i,j,k+1)%x +iblock%coordinate(i,j+1,k)%x
                            !>
                            !>
                            t2 =-iblock%coordinate(i,j,k+1)%y +iblock%coordinate(i,j+1,k)%y
                            !>
                            !>
                            t3 =-iblock%coordinate(i,j,k+1)%z +iblock%coordinate(i,j+1,k)%z
                            !>
                            !> the face center coodinate
                            !>
                            !>
                            xfc(i) = iblock%coordinate(i,j  ,k  )%x + &
                                     iblock%coordinate(i,j  ,k+1)%x + &
                                     iblock%coordinate(i,j+1,k  )%x + &
                                     iblock%coordinate(i,j+1,k+1)%x
                            !>
                            !>
                            yfc(i) = iblock%coordinate(i,j  ,k  )%y + &
                                     iblock%coordinate(i,j  ,k+1)%y + &
                                     iblock%coordinate(i,j+1,k  )%y + &
                                     iblock%coordinate(i,j+1,k+1)%y
                            !>
                            !>
                            zfc(i) = iblock%coordinate(i,j  ,k  )%z + &
                                     iblock%coordinate(i,j  ,k+1)%z + &
                                     iblock%coordinate(i,j+1,k  )%z + &
                                     iblock%coordinate(i,j+1,k+1)%z
                            !>
                            !>
                            iblock%metric(i,j,k)%fi =0.5*(t2*t6-t3*t5)
                            !>
                            !>
                            iblock%metric(i,j,k)%fj =0.5*(t3*t4-t1*t6)
                            !>
                            !>
                            iblock%metric(i,j,k)%fk =0.5*(t1*t5-t2*t4)
                            !>
                            !>
                            if(i .gt. 2) then
                                !>
                                !>
                                dxfc =xfc(i) - xfc(i-1)
                                dyfc =yfc(i) - yfc(i-1)
                                dzfc =zfc(i) - zfc(i-1)
                                !>
                                !>
                                vx  =iblock%metric(i  ,j,k)%fi + &
                                     iblock%metric(i-1,j,k)%fi
                                !>
                                !>
                                vy  =iblock%metric(i  ,j,k)%fj + &
                                     iblock%metric(i-1,j,k)%fj
                                !>
                                !>
                                vz  =iblock%metric(i  ,j,k)%fk + &
                                     iblock%metric(i-1,j,k)%fk

                                volume(i-1,j,k) = kap*abs((vx*dxfc+vy*dyfc+vz*dzfc))
                                !>
                                !>
                                !>
                            end if
                        end do
                    end do
                end do

                !> y-direciion
                do i  = 2,idim
                    do k  = 2,kdim
                        do j  = 2,jdim + 1
                            !>
                            !>
                            t4 =iblock%coordinate(i+1,j,k+1)%x - iblock%coordinate(i,j,k)%x
                            !>
                            !>
                            t5 =iblock%coordinate(i+1,j,k+1)%y - iblock%coordinate(i,j,k)%y
                            !>
                            !>
                            t6 =iblock%coordinate(i+1,j,k+1)%z - iblock%coordinate(i,j,k)%z
                            !>
                            !>
                            !>
                            t1 =-iblock%coordinate(i+1,j,k)%x + iblock%coordinate(i,j,k+1)%x
                            !>
                            !>
                            t2 =-iblock%coordinate(i+1,j,k)%y + iblock%coordinate(i,j,k+1)%y
                            !>
                            !>
                            t3 =-iblock%coordinate(i+1,j,k)%z + iblock%coordinate(i,j,k+1)%z
                            !>
                            !> the face center coodinate
                            !>
                            !>
                            xfc(j) = iblock%coordinate(i  ,j,k  )%x + &
                                     iblock%coordinate(i  ,j,k+1)%x + &
                                     iblock%coordinate(i+1,j,k  )%x + &
                                     iblock%coordinate(i+1,j,k+1)%x
                            !>
                            !>
                            yfc(j) = iblock%coordinate(i  ,j,k  )%y + &
                                     iblock%coordinate(i  ,j,k+1)%y + &
                                     iblock%coordinate(i+1,j,k  )%y + &
                                     iblock%coordinate(i+1,j,k+1)%y
                            !>
                            !>
                            zfc(j) = iblock%coordinate(i  ,j,k  )%z + &
                                     iblock%coordinate(i  ,j,k+1)%z + &
                                     iblock%coordinate(i+1,j,k  )%z + &
                                     iblock%coordinate(i+1,j,k+1)%z
                            !>
                            !>
                            iblock%metric(i,j,k)%gi =0.5*(t2*t6-t3*t5)
                            !>
                            !>
                            iblock%metric(i,j,k)%gj =0.5*(t3*t4-t1*t6)
                            !>
                            !>
                            iblock%metric(i,j,k)%gk =0.5*(t1*t5-t2*t4)
                            !>
                            !>
                            !>
                            t1 =sqrt(iblock%metric(i,j,k)%gi**2+iblock%metric(i,j,k)%gj**2+iblock%metric(i,j,k)%gk**2)
!                            write(1223,'("axis=",3i4," metric=",3e24.16)') i,j,k,iblock%metric(i,j,k)%gi/t1,iblock%metric(i,j,k)%gj/t1,iblock%metric(i,j,k)%gk/t1
                            if(j.gt.2) then
                                !>
                                !>
                                dxfc =xfc(j)-xfc(j-1)
                                dyfc =yfc(j)-yfc(j-1)
                                dzfc =zfc(j)-zfc(j-1)
                                !>
                                !>
                                vx  =iblock%metric(i,j-1,k)%gi + &
                                     iblock%metric(i,j  ,k)%gi
                                !>
                                !>
                                vy  =iblock%metric(i,j-1,k)%gj + &
                                     iblock%metric(i,j  ,k)%gj
                                !>
                                !>
                                vz  =iblock%metric(i,j-1,k)%gk + &
                                     iblock%metric(i,j  ,k)%gk

                                volume(i,j-1,k) = volume(i,j-1,k)+ kap*abs((vx*dxfc+vy*dyfc+vz*dzfc))
!                                write(1223,'("axis=",3i4," metric=",e24.16)')i,j,k,volume(i,j-1,k)
                            end if
                        end do
                    end do
                end do
                !>
                !>
                !>z-direciion
                do i  = 2,idim
                    do j  = 2,jdim
                        do k  = 2,kdim + 1
                            !>
                            !>
                            t4 =-iblock%coordinate(i+1,j+1,k)%x + iblock%coordinate(i,j,k)%x
                            !>
                            !>
                            t5 =-iblock%coordinate(i+1,j+1,k)%y + iblock%coordinate(i,j,k)%y
                            !>
                            !>
                            t6 =-iblock%coordinate(i+1,j+1,k)%z + iblock%coordinate(i,j,k)%z
                            !>
                            !>
                            !>
                            t1 =iblock%coordinate(i,j+1,k)%x - iblock%coordinate(i+1,j,k)%x
                            !>
                            !>
                            t2 =iblock%coordinate(i,j+1,k)%y - iblock%coordinate(i+1,j,k)%y
                            !>
                            !>
                            t3 =iblock%coordinate(i,j+1,k)%z - iblock%coordinate(i+1,j,k)%z
                            !>
                            !> the face center coodinate
                            !>
                            !>
                            xfc(k) = iblock%coordinate(i  ,j  ,k)%x + &
                                     iblock%coordinate(i+1,j  ,k)%x + &
                                     iblock%coordinate(i  ,j+1,k)%x + &
                                     iblock%coordinate(i+1,j+1,k)%x
                            !>
                            !>
                            yfc(k) = iblock%coordinate(i  ,j  ,k)%y + &
                                     iblock%coordinate(i+1,j  ,k)%y + &
                                     iblock%coordinate(i  ,j+1,k)%y + &
                                     iblock%coordinate(i+1,j+1,k)%y
                            !>
                            !>
                            zfc(k) = iblock%coordinate(i  ,j  ,k)%z + &
                                     iblock%coordinate(i+1,j  ,k)%z + &
                                     iblock%coordinate(i  ,j+1,k)%z + &
                                     iblock%coordinate(i+1,j+1,k)%z
                            !>
                            !>
                            iblock%metric(i,j,k)%hi =0.5*(t2*t6-t3*t5)
                            !>
                            !>
                            iblock%metric(i,j,k)%hj =0.5*(t3*t4-t1*t6)
                            !>
                            !>DAT(NM)%GG(I,J,K)
                            iblock%metric(i,j,k)%hk =0.5*(t1*t5-t2*t4)
                            !>
!                            write(1224,'("axis=",3i4," metric=",3e24.16)') i,j,k,iblock%metric(i,j,k)%hi,iblock%metric(i,j,k)%hj,iblock%metric(i,j,k)%hk
                            if(k.gt.2) then
                                !>
                                !>
                                dxfc =xfc(k)-xfc(k-1)
                                dyfc =yfc(k)-yfc(k-1)
                                dzfc =zfc(k)-zfc(k-1)
                                !>
                                !>
                                vx  =iblock%metric(i,j,k-1)%hi + &
                                     iblock%metric(i,j,k  )%hi
                                !>
                                !>
                                vy  =iblock%metric(i,j,k-1)%hj + &
                                     iblock%metric(i,j,k  )%hj
                                !>
                                !>
                                vz  =iblock%metric(i,j,k-1)%hk + &
                                     iblock%metric(i,j,k  )%hk

                                volume(i,j,k-1) = volume(i,j,k-1)+kap*abs((vx*dxfc+vy*dyfc+vz*dzfc))
!                                write(1224,'("axis=",3i4," metric=",e24.16)')i,j,k,volume(i,j,k-1)
                            end if
                        end do
                    end do
                end do
                !> boundary metric conditions
                !>
                !>
                !>
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
                            !>
                            !>
                            if(volume(i,j,k) .le. 0.0)then
                                open(unit=72,file='ccfd.error',status='unknown')
                                write(72,'("CCFDV3.0::stopping ... negative volume at block# ",i6," axis=",3i5,": the volume is negative,volume=",e24.16)') nbl,i,k,k,volume(i,j,k)
                                stop
                            end if
                            iblock%metric(i,j,k)%volume     = volume(i,j,k)
                            iblock%metric(i,j,k)%jacobian   = 1.0/volume(i,j,k)
!                            write(666,'(4i6,2e24.16)') nbl,i-1,j-1,k-1,volume(i,j,k),1.0/volume(i,j,k)
                        end do
                    end do
                end do
                !> in  x direction
                !> fill in extra values for the surface for safety
                !>
                !>
                do k = 1,kdim+2
                    do j = 1,jdim+2
                        !>
                        !>
                        !> f term
                        iblock%metric(0,j,k)%fi         = iblock%metric(2,j,k)%fi
                        iblock%metric(0,j,k)%fj         = iblock%metric(2,j,k)%fj
                        iblock%metric(0,j,k)%fk         = iblock%metric(2,j,k)%fk
!                        iblock%metric(0,j,k)%ff         = iblock%metric(2,j,k)%ff
                        !>
                        !>
                        !>
                        iblock%metric(1,j,k)%fi         = iblock%metric(2,j,k)%fi
                        iblock%metric(1,j,k)%fj         = iblock%metric(2,j,k)%fj
                        iblock%metric(1,j,k)%fk         = iblock%metric(2,j,k)%fk
!                        iblock%metric(1,j,k)%ff         = iblock%metric(2,j,k)%ff
                        !>
                        !>
                        !>
                        !>
                        iblock%metric(idim+2,j,k)%fi         = iblock%metric(idim+1,j,k)%fi
                        iblock%metric(idim+2,j,k)%fj         = iblock%metric(idim+1,j,k)%fj
                        iblock%metric(idim+2,j,k)%fk         = iblock%metric(idim+1,j,k)%fk
!                        iblock%metric(idim+2,j,k)%ff         = iblock%metric(idim+1,j,k)%ff

                        !>
                        !>
                        !>
                        iblock%metric(idim+3,j,k)%fi         = iblock%metric(idim+1,j,k)%fi
                        iblock%metric(idim+3,j,k)%fj         = iblock%metric(idim+1,j,k)%fj
                        iblock%metric(idim+3,j,k)%fk         = iblock%metric(idim+1,j,k)%fk
!                        iblock%metric(idim+3,j,k)%ff         = iblock%metric(idim+1,j,k)%ff

                        !>
                        !> g term
                        iblock%metric(0,j,k)%gi  = iblock%metric(2,j,k)%gi
                        iblock%metric(0,j,k)%gj  = iblock%metric(2,j,k)%gj
                        iblock%metric(0,j,k)%gk  = iblock%metric(2,j,k)%gk
!                        iblock%metric(0,j,k)%gg  = iblock%metric(2,j,k)%gg
                        !>
                        !>
                        !>
                        iblock%metric(1,j,k)%gi  = iblock%metric(2,j,k)%gi
                        iblock%metric(1,j,k)%gj  = iblock%metric(2,j,k)%gj
                        iblock%metric(1,j,k)%gk  = iblock%metric(2,j,k)%gk
!                        iblock%metric(1,j,k)%gg  = iblock%metric(2,j,k)%gg
                        !>
                        iblock%metric(idim+1,j,k)%gi  = iblock%metric(idim,j,k)%gi
                        iblock%metric(idim+1,j,k)%gj  = iblock%metric(idim,j,k)%gj
                        iblock%metric(idim+1,j,k)%gk  = iblock%metric(idim,j,k)%gk
                        !>
                        !>
                        iblock%metric(idim+2,j,k)%gi  = iblock%metric(idim,j,k)%gi
                        iblock%metric(idim+2,j,k)%gj  = iblock%metric(idim,j,k)%gj
                        iblock%metric(idim+2,j,k)%gk  = iblock%metric(idim,j,k)%gk
!                        iblock%metric(idim+2,j,k)%gg  = iblock%metric(idim+1,j,k)%gg
                        !>
                        !>
                        !>
                        iblock%metric(idim+3,j,k)%gi  = iblock%metric(idim,j,k)%gi
                        iblock%metric(idim+3,j,k)%gj  = iblock%metric(idim,j,k)%gj
                        iblock%metric(idim+3,j,k)%gk  = iblock%metric(idim,j,k)%gk
!                        iblock%metric(idim+3,j,k)%gg  = iblock%metric(idim+1,j,k)%gg
                        !>
                        !>
                        !>
                        !> h term
                        iblock%metric(0,j,k)%hi  = iblock%metric(2,j,k)%hi
                        iblock%metric(0,j,k)%hj  = iblock%metric(2,j,k)%hj
                        iblock%metric(0,j,k)%hk  = iblock%metric(2,j,k)%hk
!                        iblock%metric(0,j,k)%hh  = iblock%metric(2,j,k)%hh
                        !>
                        !>
                        !>
                        iblock%metric(1,j,k)%hi  = iblock%metric(2,j,k)%hi
                        iblock%metric(1,j,k)%hj  = iblock%metric(2,j,k)%hj
                        iblock%metric(1,j,k)%hk  = iblock%metric(2,j,k)%hk
!                        iblock%metric(1,j,k)%hh  = iblock%metric(2,j,k)%hh
                        !>
                        !>
                        iblock%metric(idim+1,j,k)%hi  = iblock%metric(idim,j,k)%hi
                        iblock%metric(idim+1,j,k)%hj  = iblock%metric(idim,j,k)%hj
                        iblock%metric(idim+1,j,k)%hk  = iblock%metric(idim,j,k)%hk
                        !>
                        !>
                        iblock%metric(idim+2,j,k)%hi  = iblock%metric(idim,j,k)%hi
                        iblock%metric(idim+2,j,k)%hj  = iblock%metric(idim,j,k)%hj
                        iblock%metric(idim+2,j,k)%hk  = iblock%metric(idim,j,k)%hk
!                        iblock%metric(idim+2,j,k)%ff  = iblock%metric(idim+1,j,k)%hh
                        !>
                        !>
                        iblock%metric(idim+3,j,k)%hi  = iblock%metric(idim,j,k)%hi
                        iblock%metric(idim+3,j,k)%hj  = iblock%metric(idim,j,k)%hj
                        iblock%metric(idim+3,j,k)%hk  = iblock%metric(idim,j,k)%hk
!                        iblock%metric(idim+3,j,k)%hh  = iblock%metric(idim+1,j,k)%hh
                        !>
                        !>
                        iblock%metric(0,j,k)%jacobian       = iblock%metric(2,j,k)%jacobian
                        iblock%metric(0,j,k)%volume         = iblock%metric(2,j,k)%volume
                        iblock%metric(1,j,k)%jacobian       = iblock%metric(2,j,k)%jacobian
                        iblock%metric(1,j,k)%volume         = iblock%metric(2,j,k)%volume
                        iblock%metric(idim+1,j,k)%jacobian  = iblock%metric(idim,j,k)%jacobian
                        iblock%metric(idim+1,j,k)%volume    = iblock%metric(idim,j,k)%volume
                        iblock%metric(idim+2,j,k)%jacobian  = iblock%metric(idim,j,k)%jacobian
                        iblock%metric(idim+2,j,k)%volume    = iblock%metric(idim,j,k)%volume
                        iblock%metric(idim+3,j,k)%jacobian  = iblock%metric(idim,j,k)%jacobian
                        iblock%metric(idim+3,j,k)%volume    = iblock%metric(idim,j,k)%volume
                        !>
                        !>
                    end do
                end do
                !>
                !>
                !> in  j-direction
                !> fill in extra values for the surface for safety
                !>
                do i = 1,idim+2
                    do k = 1,kdim+2
                        !>
                        !>
                        !> f term
                        iblock%metric(i,0,k)%fi          = iblock%metric(i,2,k)%fi
                        iblock%metric(i,0,k)%fj          = iblock%metric(i,2,k)%fj
                        iblock%metric(i,0,k)%fk          = iblock%metric(i,2,k)%fk
!                        iblock%metric(i,0,k)%ff          = iblock%metric(i,2,k)%ff
                        !>
                        !>
                        !>
                        iblock%metric(i,1,k)%fi          = iblock%metric(i,2,k)%fi
                        iblock%metric(i,1,k)%fj          = iblock%metric(i,2,k)%fj
                        iblock%metric(i,1,k)%fk          = iblock%metric(i,2,k)%fk
!                        iblock%metric(i,1,k)%ff          = iblock%metric(i,2,k)%ff
                        !>
                        !>
                        iblock%metric(i,jdim+1,k)%fi          = iblock%metric(i,jdim,k)%fi
                        iblock%metric(i,jdim+1,k)%fj          = iblock%metric(i,jdim,k)%fj
                        iblock%metric(i,jdim+1,k)%fk          = iblock%metric(i,jdim,k)%fk
                        !>
                        !>
                        iblock%metric(i,jdim+2,k)%fi          = iblock%metric(i,jdim,k)%fi
                        iblock%metric(i,jdim+2,k)%fj          = iblock%metric(i,jdim,k)%fj
                        iblock%metric(i,jdim+2,k)%fk          = iblock%metric(i,jdim,k)%fk
!                        iblock%metric(i,jdim+2,k)%ff          = iblock%metric(i,jdim+1,k)%ff
                        !>
                        !>
                        iblock%metric(i,jdim+3,k)%fi          = iblock%metric(i,jdim,k)%fi
                        iblock%metric(i,jdim+3,k)%fj          = iblock%metric(i,jdim,k)%fj
                        iblock%metric(i,jdim+3,k)%fk          = iblock%metric(i,jdim,k)%fk
!                        iblock%metric(i,jdim+3,k)%ff          = iblock%metric(i,jdim+1,k)%ff
                        !>
                        !>
                        !>
                        !> g term
                        iblock%metric(i,0,k)%gi  = iblock%metric(i,2,k)%gi
                        iblock%metric(i,0,k)%gj  = iblock%metric(i,2,k)%gj
                        iblock%metric(i,0,k)%gk  = iblock%metric(i,2,k)%gk
!                        iblock%metric(i,0,k)%gg  = iblock%metric(i,2,k)%gg
                        !>
                        !>
                        !>
                        iblock%metric(i,1,k)%gi  = iblock%metric(i,2,k)%gi
                        iblock%metric(i,1,k)%gj  = iblock%metric(i,2,k)%gj
                        iblock%metric(i,1,k)%gk  = iblock%metric(i,2,k)%gk
!                        iblock%metric(i,1,k)%gg  = iblock%metric(i,2,k)%gg
                        !>
                        !>
                        !>
                        iblock%metric(i,jdim+2,k)%gi  = iblock%metric(i,jdim+1,k)%gi
                        iblock%metric(i,jdim+2,k)%gj  = iblock%metric(i,jdim+1,k)%gj
                        iblock%metric(i,jdim+2,k)%gk  = iblock%metric(i,jdim+1,k)%gk
!                        iblock%metric(i,jdim+2,k)%gg  = iblock%metric(i,jdim+1,k)%gg
                        !>
                        !>
                        !>
                        iblock%metric(i,jdim+3,k)%gi  = iblock%metric(i,jdim+1,k)%gi
                        iblock%metric(i,jdim+3,k)%gj  = iblock%metric(i,jdim+1,k)%gj
                        iblock%metric(i,jdim+3,k)%gk  = iblock%metric(i,jdim+1,k)%gk
!                        iblock%metric(i,jdim+3,k)%gg  = iblock%metric(i,jdim+1,k)%gg
                                                                        !>
                        !> h term
                        iblock%metric(i,0,k)%hi  = iblock%metric(i,2,k)%hi
                        iblock%metric(i,0,k)%hj  = iblock%metric(i,2,k)%hj
                        iblock%metric(i,0,k)%hk  = iblock%metric(i,2,k)%hk
!                        iblock%metric(i,0,k)%hh  = iblock%metric(i,2,k)%hh
                        !>
                        !>
                        !>
                        iblock%metric(i,1,k)%hi  = iblock%metric(i,2,k)%hi
                        iblock%metric(i,1,k)%hj  = iblock%metric(i,2,k)%hj
                        iblock%metric(i,1,k)%hk  = iblock%metric(i,2,k)%hk
!                        iblock%metric(i,1,k)%hh  = iblock%metric(i,2,k)%hh
                        !>
                        !>
                        iblock%metric(i,jdim+1,k)%hi  = iblock%metric(i,jdim,k)%hi
                        iblock%metric(i,jdim+1,k)%hj  = iblock%metric(i,jdim,k)%hj
                        iblock%metric(i,jdim+1,k)%hk  = iblock%metric(i,jdim,k)%hk
                        !>
                        !>
                        iblock%metric(i,jdim+2,k)%hi  = iblock%metric(i,jdim,k)%hi
                        iblock%metric(i,jdim+2,k)%hj  = iblock%metric(i,jdim,k)%hj
                        iblock%metric(i,jdim+2,k)%hk  = iblock%metric(i,jdim,k)%hk
!                        iblock%metric(i,jdim+2,k)%hh  = iblock%metric(i,jdim+1,k)%hh
                        !>
                        !>
                        iblock%metric(i,jdim+3,k)%hi  = iblock%metric(i,jdim,k)%hi
                        iblock%metric(i,jdim+3,k)%hj  = iblock%metric(i,jdim,k)%hj
                        iblock%metric(i,jdim+3,k)%hk  = iblock%metric(i,jdim,k)%hk
!                        iblock%metric(i,jdim+3,k)%hh  = iblock%metric(i,jdim+1,k)%hh
                        !>
                        !>
                        iblock%metric(i,0,k)%jacobian       = iblock%metric(i,2,k)%jacobian
                        iblock%metric(i,0,k)%volume         = iblock%metric(i,2,k)%volume
                        iblock%metric(i,1,k)%jacobian       = iblock%metric(i,2,k)%jacobian
                        iblock%metric(i,1,k)%volume         = iblock%metric(i,2,k)%volume
                        iblock%metric(i,jdim+1,k)%jacobian  = iblock%metric(i,jdim,k)%jacobian
                        iblock%metric(i,jdim+1,k)%volume    = iblock%metric(i,jdim,k)%volume
                        iblock%metric(i,jdim+2,k)%jacobian  = iblock%metric(i,jdim,k)%jacobian
                        iblock%metric(i,jdim+2,k)%volume    = iblock%metric(i,jdim,k)%volume
                        iblock%metric(i,jdim+3,k)%jacobian  = iblock%metric(i,jdim,k)%jacobian
                        iblock%metric(i,jdim+3,k)%volume    = iblock%metric(i,jdim,k)%volume
                    end do
                end do
                !>
                !>
                !> in  j  direction
                !> fill in extra values for the surface for safety
                do i = 1,idim+2
                    do j = 1,jdim+2
                        !>
                        !>
                        !> f term
                        iblock%metric(i,j,0)%fi         = iblock%metric(i,j,2)%fi
                        iblock%metric(i,j,0)%fj         = iblock%metric(i,j,2)%fj
                        iblock%metric(i,j,0)%fk         = iblock%metric(i,j,2)%fk
!                        iblock%metric(i,j,0)%ff         = iblock%metric(i,j,2)%ff
                        !>
                        !>
                        !>
                        iblock%metric(i,j,1)%fi          = iblock%metric(i,j,2)%fi
                        iblock%metric(i,j,1)%fj          = iblock%metric(i,j,2)%fj
                        iblock%metric(i,j,1)%fk          = iblock%metric(i,j,2)%fk
!                        iblock%metric(i,j,1)%ff          = iblock%metric(i,j,2)%ff
                        !>
                        !>
                        iblock%metric(i,j,kdim+1)%fi         = iblock%metric(i,j,kdim)%fi
                        iblock%metric(i,j,kdim+1)%fj         = iblock%metric(i,j,kdim)%fj
                        iblock%metric(i,j,kdim+1)%fk         = iblock%metric(i,j,kdim)%fk
                        !>
                        !>
                        iblock%metric(i,j,kdim+2)%fi         = iblock%metric(i,j,kdim)%fi
                        iblock%metric(i,j,kdim+2)%fj         = iblock%metric(i,j,kdim)%fj
                        iblock%metric(i,j,kdim+2)%fk         = iblock%metric(i,j,kdim)%fk
!                        iblock%metric(i,j,kdim+2)%ff         = iblock%metric(i,j,kdim+1)%ff
                        !>
                        !>
                        iblock%metric(i,j,kdim+3)%fi         = iblock%metric(i,j,kdim)%fi
                        iblock%metric(i,j,kdim+3)%fj         = iblock%metric(i,j,kdim)%fj
                        iblock%metric(i,j,kdim+3)%fk         = iblock%metric(i,j,kdim)%fk
!                        iblock%metric(i,j,kdim+3)%ff         = iblock%metric(i,j,kdim+1)%ff
                        !>
                        !>
                        !> g term
                        iblock%metric(i,j,0)%gi  = iblock%metric(i,j,2)%gi
                        iblock%metric(i,j,0)%gj  = iblock%metric(i,j,2)%gj
                        iblock%metric(i,j,0)%gk  = iblock%metric(i,j,2)%gk
!                        iblock%metric(i,j,0)%gg  = iblock%metric(i,j,2)%gg
                        !>
                        !>
                        !>
                        iblock%metric(i,j,1)%gi  = iblock%metric(i,j,2)%gi
                        iblock%metric(i,j,1)%gj  = iblock%metric(i,j,2)%gj
                        iblock%metric(i,j,1)%gk  = iblock%metric(i,j,2)%gk
!                        iblock%metric(i,j,1)%gg  = iblock%metric(i,j,2)%gg
                        !>
                        !>
                        iblock%metric(i,j,kdim+1)%gi  = iblock%metric(i,j,kdim)%gi
                        iblock%metric(i,j,kdim+1)%gj  = iblock%metric(i,j,kdim)%gj
                        iblock%metric(i,j,kdim+1)%gk  = iblock%metric(i,j,kdim)%gk
                        !>
                        !>
                        iblock%metric(i,j,kdim+2)%gi  = iblock%metric(i,j,kdim)%gi
                        iblock%metric(i,j,kdim+2)%gj  = iblock%metric(i,j,kdim)%gj
                        iblock%metric(i,j,kdim+2)%gk  = iblock%metric(i,j,kdim)%gk
!                        iblock%metric(i,j,kdim+2)%gg  = iblock%metric(i,j,kdim+1)%gg
                        !>
                        !>
                        iblock%metric(i,j,kdim+3)%gi  = iblock%metric(i,j,kdim)%gi
                        iblock%metric(i,j,kdim+3)%gj  = iblock%metric(i,j,kdim)%gj
                        iblock%metric(i,j,kdim+3)%gk  = iblock%metric(i,j,kdim)%gk
!                        iblock%metric(i,j,kdim+3)%gg  = iblock%metric(i,j,kdim+1)%gg
                        !>
                        !>
                        !>
                                                                        !>
                        !> h term
                        iblock%metric(i,j,0)%hi  = iblock%metric(i,j,2)%hi
                        iblock%metric(i,j,0)%hj  = iblock%metric(i,j,2)%hj
                        iblock%metric(i,j,0)%hk  = iblock%metric(i,j,2)%hk
!                        iblock%metric(i,j,0)%hh  = iblock%metric(i,j,2)%hh
                        !>
                        !>
                        !>
                        iblock%metric(i,j,1)%hi  = iblock%metric(i,j,2)%hi
                        iblock%metric(i,j,1)%hj  = iblock%metric(i,j,2)%hj
                        iblock%metric(i,j,1)%hk  = iblock%metric(i,j,2)%hk
!                        iblock%metric(i,j,1)%hh  = iblock%metric(i,j,2)%hh
                        !>
                        !>
                        !>
                        iblock%metric(i,j,kdim+2)%hi  = iblock%metric(i,j,kdim+1)%hi
                        iblock%metric(i,j,kdim+2)%hj  = iblock%metric(i,j,kdim+1)%hj
                        iblock%metric(i,j,kdim+2)%hk  = iblock%metric(i,j,kdim+1)%hk
!                        iblock%metric(i,j,kdim+2)%hh  = iblock%metric(i,j,kdim+1)%hh
                        !>
                        !>
                        !>
                        iblock%metric(i,j,kdim+3)%hi  = iblock%metric(i,j,kdim+1)%hi
                        iblock%metric(i,j,kdim+3)%hj  = iblock%metric(i,j,kdim+1)%hj
                        iblock%metric(i,j,kdim+3)%hk  = iblock%metric(i,j,kdim+1)%hk
!                        iblock%metric(i,j,kdim+3)%hh  = iblock%metric(i,j,kdim+1)%hh
                        !>
                        !>
                        iblock%metric(i,j,0)%jacobian       = iblock%metric(i,j,2)%jacobian
                        iblock%metric(i,j,0)%volume         = iblock%metric(i,j,2)%volume
                        iblock%metric(i,j,1)%jacobian       = iblock%metric(i,j,2)%jacobian
                        iblock%metric(i,j,1)%volume         = iblock%metric(i,j,2)%volume
                        iblock%metric(i,j,kdim+1)%jacobian  = iblock%metric(i,j,kdim)%jacobian
                        iblock%metric(i,j,kdim+1)%volume    = iblock%metric(i,j,kdim)%volume
                        iblock%metric(i,j,kdim+2)%jacobian  = iblock%metric(i,j,kdim)%jacobian
                        iblock%metric(i,j,kdim+2)%volume    = iblock%metric(i,j,kdim)%volume
                        iblock%metric(i,j,kdim+3)%jacobian  = iblock%metric(i,j,kdim)%jacobian
                        iblock%metric(i,j,kdim+3)%volume    = iblock%metric(i,j,kdim)%volume
                    end do
                end do
                !>
                !>
                !>
                !>
                !>
                do k=0,kdim+3
                    do j=0,jdim+3
                        do i=0,idim+3
                            iblock%metric(i,j,k)%ff = sqrt(iblock%metric(i,j,k)%fi**2 + &
                                                           iblock%metric(i,j,k)%fj**2 + &
                                                           iblock%metric(i,j,k)%fk**2 )
                            iblock%metric(i,j,k)%gg = sqrt(iblock%metric(i,j,k)%gi**2 + &
                                                           iblock%metric(i,j,k)%gj**2 + &
                                                           iblock%metric(i,j,k)%gk**2 )
                            iblock%metric(i,j,k)%hh = sqrt(iblock%metric(i,j,k)%hi**2 + &
                                                           iblock%metric(i,j,k)%hj**2 + &
                                                           iblock%metric(i,j,k)%hk**2 )
                            !>
                            !>
                            !>
                            !>
                            iblock%metric(i,j,k)%fi = iblock%metric(i,j,k)%fi /iblock%metric(i,j,k)%ff
                            iblock%metric(i,j,k)%fj = iblock%metric(i,j,k)%fj /iblock%metric(i,j,k)%ff
                            iblock%metric(i,j,k)%fk = iblock%metric(i,j,k)%fk /iblock%metric(i,j,k)%ff
                            !>
                            !>
                            iblock%metric(i,j,k)%gi = iblock%metric(i,j,k)%gi /iblock%metric(i,j,k)%gg
                            iblock%metric(i,j,k)%gj = iblock%metric(i,j,k)%gj /iblock%metric(i,j,k)%gg
                            iblock%metric(i,j,k)%gk = iblock%metric(i,j,k)%gk /iblock%metric(i,j,k)%gg
                            !>
                            !>
                            iblock%metric(i,j,k)%hi = iblock%metric(i,j,k)%hi /iblock%metric(i,j,k)%hh
                            iblock%metric(i,j,k)%hj = iblock%metric(i,j,k)%hj /iblock%metric(i,j,k)%hh
                            iblock%metric(i,j,k)%hk = iblock%metric(i,j,k)%hk /iblock%metric(i,j,k)%hh

                        end do
                    end do
                end do
                !>
                !>
!                write(666,'("nbl=",i6," axis=",3i6)') nbl,idim-1,jdim-1,kdim-1

                do k = 2,kdim
                    do j = 2,jdim
!                        iblock%metric(0,j,k)%volume     = iblock%metric(2,j,k)%volume
!                        iblock%metric(0,j,k)%jacobian   = iblock%metric(2,j,k)%jacobian
                        !>
                        iblock%metric(1,j,k)%volume     = iblock%metric(2,j,k)%volume
                        iblock%metric(1,j,k)%jacobian   = iblock%metric(2,j,k)%jacobian
                        !>
                        iblock%metric(idim+1,j,k)%volume     = iblock%metric(idim,j,k)%volume
                        iblock%metric(idim+1,j,k)%jacobian   = iblock%metric(idim,j,k)%jacobian
                        !>
!                        iblock%metric(idim+2,j,k)%volume     = iblock%metric(idim,j,k)%volume
!                        iblock%metric(idim+2,j,k)%jacobian   = iblock%metric(idim,j,k)%jacobian
                    end do
                end do
                do i = 2,idim
                    do k = 2,kdim
!                        iblock%metric(i,0,k)%volume      = iblock%metric(i,2,k)%volume
!                        iblock%metric(i,0,k)%jacobian    = iblock%metric(i,2,k)%jacobian
                        !>
                        iblock%metric(i,1,k)%volume      = iblock%metric(i,2,k)%volume
                        iblock%metric(i,1,k)%jacobian    = iblock%metric(i,2,k)%jacobian
                        !>
                        iblock%metric(i,jdim+1,k)%volume      = iblock%metric(i,jdim,k)%volume
                        iblock%metric(i,jdim+1,k)%jacobian    = iblock%metric(i,jdim,k)%jacobian
                        !>
!                        iblock%metric(i,jdim+2,k)%volume      = iblock%metric(i,jdim,k)%volume
!                        iblock%metric(i,jdim+2,k)%jacobian    = iblock%metric(i,jdim,k)%jacobian
                    end do
                end do
                do i = 2,idim
                    do j = 2,jdim
!                        iblock%metric(i,j,0)%volume     = iblock%metric(i,j,2)%volume
!                        iblock%metric(i,j,0)%jacobian   = iblock%metric(i,j,2)%jacobian
                        !>
                        iblock%metric(i,j,1)%volume      = iblock%metric(i,j,2)%volume
                        iblock%metric(i,j,1)%jacobian    = iblock%metric(i,j,2)%jacobian
                        !>
                        iblock%metric(i,j,kdim+1)%volume     = iblock%metric(i,j,kdim)%volume
                        iblock%metric(i,j,kdim+1)%jacobian   = iblock%metric(i,j,kdim)%jacobian
                        !>
!                        iblock%metric(i,j,kdim+2)%volume     = iblock%metric(i,j,kdim)%volume
!                        iblock%metric(i,j,kdim+2)%jacobian   = iblock%metric(i,j,kdim)%jacobian
                    end do
                end do
                !>
                deallocate(volume)
#if defined PMPI
            endif
#endif
        end do
        !>
        !>
        if(mgflags .ne. 0)then
        !>
        !>
            do nbl =1,nblocks,global_level
#if defined PMPI
                if(myid .eq. n2p(nbl) .or. myid .eq. myhost)then
#endif
                    !>
                    !> restrict the volumes to coarser meshes form finer meshes
                    !>
                    do ilevel=1,global_level-1
                        iblock => mesh%blocks(nbl+ilevel-1)
                        kk = 1
                        do k=2,iblock%kdim+1,2
                            jj = 1
                            kk = kk + 1
                            do j=2,iblock%jdim+1,2
                                ii = 1
                                jj = jj + 1
                                do i=2,iblock%idim+1,2
                                    ii = ii + 1
                                    mesh%blocks(nbl+ilevel)%metric(ii,jj,kk)%volume =   iblock%metric(i  ,j  ,k  )%volume +  &
                                                                                        iblock%metric(i+1,j  ,k  )%volume +  &
                                                                                        iblock%metric(i  ,j+1,k  )%volume +  &
                                                                                        iblock%metric(i  ,j  ,k+1)%volume +  &
                                                                                        iblock%metric(i+1,j+1,k  )%volume +  &
                                                                                        iblock%metric(i+1,j  ,k+1)%volume +  &
                                                                                        iblock%metric(i  ,j+1,k+1)%volume +  &
                                                                                        iblock%metric(i+1,j+1,k+1)%volume
                                    mesh%blocks(nbl+ilevel)%metric(ii,jj,kk)%jacobian = 1.0 / mesh%blocks(nbl+ilevel)%metric(ii,jj,kk)%volume
                                end do
                            end do
                        end do
                        !>
                        !> fill the ghost cell at the i-face
                        !>
                        i = mesh%blocks(nbl+ilevel)%idim
                        do k=2,mesh%blocks(nbl+ilevel)%kdim
                            do j=2,mesh%blocks(nbl+ilevel)%jdim
                                mesh%blocks(nbl+ilevel)%metric(1,j,k)%volume   =  mesh%blocks(nbl+ilevel)%metric(2,j,k)%volume
                                mesh%blocks(nbl+ilevel)%metric(1,j,k)%jacobian =  mesh%blocks(nbl+ilevel)%metric(2,j,k)%jacobian
                                !>
                                mesh%blocks(nbl+ilevel)%metric(i+1,j,k)%volume   =  mesh%blocks(nbl+ilevel)%metric(i,j,k)%volume
                                mesh%blocks(nbl+ilevel)%metric(i+1,j,k)%jacobian =  mesh%blocks(nbl+ilevel)%metric(i,j,k)%jacobian
                            end do
                        end do
                        !>
                        !> fill the ghost cell at the j-face
                        !>
                        j = mesh%blocks(nbl+ilevel)%jdim
                        do k=2,mesh%blocks(nbl+ilevel)%kdim
                            do i=2,mesh%blocks(nbl+ilevel)%idim
                                mesh%blocks(nbl+ilevel)%metric(i,1,k)%volume   =  mesh%blocks(nbl+ilevel)%metric(i,2,k)%volume
                                mesh%blocks(nbl+ilevel)%metric(i,1,k)%jacobian =  mesh%blocks(nbl+ilevel)%metric(i,2,k)%jacobian
                                !>
                                mesh%blocks(nbl+ilevel)%metric(i,j+1,k)%volume   =  mesh%blocks(nbl+ilevel)%metric(i,j,k)%volume
                                mesh%blocks(nbl+ilevel)%metric(i,j+1,k)%jacobian =  mesh%blocks(nbl+ilevel)%metric(i,j,k)%jacobian
                            end do
                        end do
                        !>
                        !> fill the ghost cell at the k-face
                        !>
                        k = mesh%blocks(nbl+ilevel)%kdim
                        do j=2,mesh%blocks(nbl+ilevel)%jdim
                            do i=2,mesh%blocks(nbl+ilevel)%idim
                                mesh%blocks(nbl+ilevel)%metric(i,j,1)%volume   =  mesh%blocks(nbl+ilevel)%metric(i,j,2)%volume
                                mesh%blocks(nbl+ilevel)%metric(i,j,1)%jacobian =  mesh%blocks(nbl+ilevel)%metric(i,j,2)%jacobian
                                !>
                                mesh%blocks(nbl+ilevel)%metric(i,j,k+1)%volume   =  mesh%blocks(nbl+ilevel)%metric(i,j,k)%volume
                                mesh%blocks(nbl+ilevel)%metric(i,j,k+1)%jacobian =  mesh%blocks(nbl+ilevel)%metric(i,j,k)%jacobian
                            end do
                        end do
                        !>
                        !>
                    end do
        !>
#if defined PMPI
                end if
#endif
            end do
    !>
        endif

        return
    end subroutine metric


