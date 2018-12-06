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
    subroutine lu_adi(nbl)

        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
        use coordinate_module
        use variable_module
        use metric_module
        use turbulent_module
        !>
        !>
		implicit none
		!>
		!>
        type(overlap_type),pointer :: mesh
        type(blocks_type),pointer  :: iblock
        !>
        !>
        real(kind = dprec),dimension(:,:,:,:),allocatable :: flmn,flmp,diag
        real(kind = dprec),dimension(:,:,:)  ,allocatable :: flam1
        real(kind = dprec),dimension(:,:,:)  ,allocatable :: abs_flam1
        real(kind = dprec),dimension(:,:,:)  ,allocatable :: flam2
        real(kind = dprec),dimension(:,:,:)  ,allocatable :: abs_flam2
        real(kind = dprec),dimension(:,:,:)  ,allocatable :: flam3
        real(kind = dprec),dimension(:,:,:)  ,allocatable :: abs_flam3
        real(kind = dprec),dimension(:,:,:)  ,allocatable :: implicit_viscous_term
        !> unit vectors in the plane
        real(kind = dprec),dimension(:,:,:)  ,allocatable :: kx,ky,kz,lx,ly,lz,mx,my,mz
        !>
        !>
        real(kind = dprec) :: avi0,const,rhs1,rhs2,rhs3,rhs4,rhs5
        real(kind = dprec) :: epsaa,epsbb,epscc,zero,jo,ca,cc,fai
        real(kind = dprec) :: dt,c1,c2,c3,c4,c5,t1,t2,t3,t4,t5
        real(kind = dprec) :: ra,u1,v1,w1,p1,&
                              xyz
        real(kind= dprec) :: uban,vban,wban
        !>
        integer :: i,j,k,&
                   nbl,&
                   idim,jdim,kdim,&
                   max_ijk
        !>
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)

        !>
        !>
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        !>
        !>
        max_ijk = 0
        max_ijk = max(idim,jdim)
        max_ijk = max(max_ijk,kdim)+1
        !>
        !>
        allocate(flmn(idim+1,jdim+1,kdim+1,3))
        allocate(flmp(idim+1,jdim+1,kdim+1,3))
        allocate(diag(idim+1,jdim+1,kdim+1,3))
        allocate(flam1(idim+1,jdim+1,kdim+1))
        allocate(abs_flam1(idim+1,jdim+1,kdim+1))
        allocate(flam2(idim+1,jdim+1,kdim+1))
        allocate(abs_flam2(idim+1,jdim+1,kdim+1))
        allocate(flam3(idim+1,jdim+1,kdim+1))
        allocate(abs_flam3(idim+1,jdim+1,kdim+1))
        allocate(implicit_viscous_term(idim+1,jdim+1,kdim+1))
        allocate(kx(idim+1,jdim+1,kdim+1))
        allocate(ky(idim+1,jdim+1,kdim+1))
        allocate(kz(idim+1,jdim+1,kdim+1))
        allocate(lx(idim+1,jdim+1,kdim+1))
        allocate(ly(idim+1,jdim+1,kdim+1))
        allocate(lz(idim+1,jdim+1,kdim+1))
        allocate(mx(idim+1,jdim+1,kdim+1))
        allocate(my(idim+1,jdim+1,kdim+1))
        allocate(mz(idim+1,jdim+1,kdim+1))

        zero   = 1.e-15
        implicit_viscous_term = 0.0
        !>
        !>
        !>
        if(nvisc .ne. 0) then
            avi0 = 2.0*rmre
        else
            avi0 =0.
        end if

        const =1./sqrt(2.0)
        !>
        !>
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    iblock%variable(i,j,k)%res_1 = -iblock%variable(i,j,k)%res_1
                    iblock%variable(i,j,k)%res_2 = -iblock%variable(i,j,k)%res_2
                    iblock%variable(i,j,k)%res_3 = -iblock%variable(i,j,k)%res_3
                    iblock%variable(i,j,k)%res_4 = -iblock%variable(i,j,k)%res_4
                    iblock%variable(i,j,k)%res_5 = -iblock%variable(i,j,k)%res_5
                end do
            end do
        end do
!       --------------------------------------------------------------
!       j-direction
!       --------------------------------------------------------------
        do k  = 2,kdim
            do i  = 2,idim
                do j = 2,jdim

                    ra     = iblock%cell(i,j,k)%r
                    u1     = iblock%cell(i,j,k)%u
                    v1     = iblock%cell(i,j,k)%v
                    w1     = iblock%cell(i,j,k)%w
                    !>
                    !>
                    kx(i,j,k)  = iblock%metric(i,j,k)%gi+iblock%metric(i,j+1,k)%gi
                    ky(i,j,k)  = iblock%metric(i,j,k)%gj+iblock%metric(i,j+1,k)%gj
                    kz(i,j,k)  = iblock%metric(i,j,k)%gk+iblock%metric(i,j+1,k)%gk
                    !>
                    !>
                    xyz = kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k) + kz(i,j,k)*kz(i,j,k)
                    xyz = 1.e0/sqrt(xyz)
                    kx(i,j,k)   = xyz*kx(i,j,k)
                    ky(i,j,k)   = xyz*ky(i,j,k)
                    kz(i,j,k)   = xyz*kz(i,j,k)
                    !>
                    !>
                    lx(i,j,k)  = ky(i,j,k) - kz(i,j,k)
                    ly(i,j,k)  = kz(i,j,k) - kx(i,j,k)
                    lz(i,j,k)  = kx(i,j,k) - ky(i,j,k)
                    xyz = lx(i,j,k)*lx(i,j,k) + ly(i,j,k)*ly(i,j,k) + lz(i,j,k)*lz(i,j,k)
                    xyz = 1.e0/sqrt(xyz)
                    lx(i,j,k) = xyz*lx(i,j,k)
                    ly(i,j,k) = xyz*ly(i,j,k)
                    lz(i,j,k) = xyz*lz(i,j,k)
                    !>
                    !>
                    mx(i,j,k) = ky(i,j,k)*lz(i,j,k) - ly(i,j,k)*kz(i,j,k)
                    my(i,j,k) = kz(i,j,k)*lx(i,j,k) - lz(i,j,k)*kx(i,j,k)
                    mz(i,j,k) = kx(i,j,k)*ly(i,j,k) - lx(i,j,k)*ky(i,j,k)
                    !>
                    !>
                    ca  =sqrt(gamma*iblock%cell(i,j,k)%p/ra)
                    c2  =ca*ca
                    uban =kx(i,j,k)*u1+ky(i,j,k)*v1+kz(i,j,k)*w1
                    !>
                    !>
                    flam1(i,j,k)     =uban
                    flam2(i,j,k)     =uban+ca
                    flam3(i,j,k)     =uban-ca
                    abs_flam1(i,j,k) = abs(uban)
                    abs_flam2(i,j,k) = abs(uban+ca)
                    abs_flam3(i,j,k) = abs(uban-ca)
                    !>
                    !>
                    !> used the limiter eigenvalues a la-harten and gnoffo
                    !> implicit tvd schemes for hyperbolic conservation laws in curvilinear coordinates
                    !> the shock-caputring schemes
                    !>
                    !> 0.01  < epsa_r < 0.4
                    !> 0.1   < epsaa  < 0.3     approximate range for blunt bodies
                    !> 0.005 < epsaa  < 0.05    approximate range for slender bodies
                    if(epsa_r .gt. 0)then
!                        write(*,*) 'run the epsa'
                        epsaa = 2.0*epsa_r*(ca + abs(u1) + abs(v1) + abs(w1))
                        epsbb = 0.25/max(epsaa,zero)
                        epscc = 2.00*epsaa
                        if(real(abs_flam1(i,j,k)) .lt. real(epscc)) abs_flam1(i,j,k)= abs_flam1(i,j,k)*abs_flam1(i,j,k)*epsbb+epsaa
                        if(real(abs_flam2(i,j,k)) .lt. real(epscc)) abs_flam2(i,j,k)= abs_flam2(i,j,k)*abs_flam2(i,j,k)*epsbb+epsaa
                        if(real(abs_flam3(i,j,k)) .lt. real(epscc)) abs_flam3(i,j,k)= abs_flam3(i,j,k)*abs_flam3(i,j,k)*epsbb+epsaa
                    end if

                end do
            end do
        end do
        !>
        !>
        do k=2,kdim
            do i=2,idim
                do j=2,jdim+1
                    !>
                    !>
                    !>
                    !>implicit treatment of viscous terms
                    !>
                    if(nvisc == 0)then
                        jo = 0.0
                    elseif( nvisc == 1)then
                        if(j == 2 )then
                            jo = iblock%metric(i,j,k)%jacobian/iblock%cell(i,j,k)%r
                        elseif(j == jdim+1)then
                            jo = iblock%metric(i,j-1,k)%jacobian/iblock%cell(i,j-1,k)%r
                        else
                            jo = (1.0/iblock%cell(i,j,k)%r+ 1.0/iblock%cell(i,j-1,k)%r)/(iblock%metric(i,j,k)%volume+iblock%metric(i,j-1,k)%volume)
                        end if
                    elseif( nvisc .gt. 1)then
                        if(j==2 )then
                            jo =(1.0+iblock%turbulent(i,j,k)%viscous)*iblock%metric(i,j,k)%jacobian/iblock%cell(i,j,k)%r
                        elseif(j == jdim +1)then
                            jo =(1.0+iblock%turbulent(i,j-1,k)%viscous)*iblock%metric(i,j-1,k)%jacobian/iblock%cell(i,j-1,k)%r
                        else
                            jo =((1.0+iblock%turbulent(i,j,k)%viscous)/iblock%cell(i,j,k)%r + (1.0+iblock%turbulent(i,j-1,k)%viscous)/iblock%cell(i,j-1,k)%r)/(iblock%metric(i,j,k)%volume+iblock%metric(i,j-1,k)%volume)
                        end if
                    else
                        write(*,*) 'have a wrong turbulent model!please check the mode used and modify!'
                    end if
                    !>
                    !>
                    implicit_viscous_term(i,j,k)      = iblock%metric(i,j,k)%gg*iblock%metric(i,j,k)%gg*jo*avi0
                    !>
                end do
                !>
                !>
                do j =3,jdim
!                   !>
!                   !> low-matrix
                    flmp(i,j,k,1) =-0.5*iblock%metric(i,j,k)%gg*(flam1(i,j-1,k)+abs_flam1(i,j-1,k))-implicit_viscous_term(i,j,k)
                    flmp(i,j,k,2) =-0.5*iblock%metric(i,j,k)%gg*(flam2(i,j-1,k)+abs_flam2(i,j-1,k))-implicit_viscous_term(i,j,k)
                    flmp(i,j,k,3) =-0.5*iblock%metric(i,j,k)%gg*(flam3(i,j-1,k)+abs_flam3(i,j-1,k))-implicit_viscous_term(i,j,k)
                    !> uper-matrix
                    flmn(i,j-1,k,1) = 0.5*iblock%metric(i,j,k)%gg*(flam1(i,j,k)-abs_flam1(i,j,k))-implicit_viscous_term(i,j,k)
                    flmn(i,j-1,k,2) = 0.5*iblock%metric(i,j,k)%gg*(flam2(i,j,k)-abs_flam2(i,j,k))-implicit_viscous_term(i,j,k)
                    flmn(i,j-1,k,3) = 0.5*iblock%metric(i,j,k)%gg*(flam3(i,j,k)-abs_flam3(i,j,k))-implicit_viscous_term(i,j,k)
                end do
                do j=2,jdim
                    !> diag-matrix
                    dt = iblock%cell(i,j,k)%dt
                    diag(i,j,k,1)   = iblock%metric(i,j,k)%volume/dt +  0.5 * iblock%metric(i,j+1,k)%gg * (flam1(i,j,k)+abs_flam1(i,j,k)) - &
                                                                        0.5 * iblock%metric(i,j  ,k)%gg * (flam1(i,j,k)-abs_flam1(i,j,k)) + &
                                                                        implicit_viscous_term(i,j,k) + implicit_viscous_term(i,j+1,k)
                    !>
                    !>
                    diag(i,j,k,2)   = iblock%metric(i,j,k)%volume/dt +  0.5 * iblock%metric(i,j+1,k)%gg * (flam2(i,j,k)+abs_flam2(i,j,k)) - &
                                                                        0.5 * iblock%metric(i,j  ,k)%gg * (flam2(i,j,k)-abs_flam2(i,j,k)) + &
                                                                        implicit_viscous_term(i,j,k) + implicit_viscous_term(i,j+1,k)
                    !>
                    !>
                    diag(i,j,k,3)   = iblock%metric(i,j,k)%volume/dt +  0.5 * iblock%metric(i,j+1,k)%gg * (flam3(i,j,k)+abs_flam3(i,j,k)) - &
                                                                        0.5 * iblock%metric(i,j  ,k)%gg * (flam3(i,j,k)-abs_flam3(i,j,k)) + &
                                                                        implicit_viscous_term(i,j,k) + implicit_viscous_term(i,j+1,k)
                    !>
                end do
            end do
        end do
        !>
        !>
        !> Multiply the inverse of the diagonalizing matrix
        !> T times the residual contribution.
        !> T(inverse)*RES
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    ra     = iblock%cell(i,j,k)%r
                    u1     = iblock%cell(i,j,k)%u
                    v1     = iblock%cell(i,j,k)%v
                    w1     = iblock%cell(i,j,k)%w
                    uban      = u1*kx(i,j,k) + v1*ky(i,j,k) + w1*kz(i,j,k)
                    vban      = u1*lx(i,j,k) + v1*ly(i,j,k) + w1*lz(i,j,k)
                    wban      = u1*mx(i,j,k) + v1*my(i,j,k) + w1*mz(i,j,k)
                    fai       = 0.5*(u1**2 + v1**2 + w1**2)
                    ca   = sqrt(gamma*iblock%cell(i,j,k)%p/ra)
                    c2   = ca*ca
                    !>
                    !>
                    t3    =  1.e0/ca
                    t1    = -fai*iblock%variable(i,j,k)%res_1 - &
                                 iblock%variable(i,j,k)%res_5 + &
                                 u1*iblock%variable(i,j,k)%res_2 + &
                                 v1*iblock%variable(i,j,k)%res_3 + &
                                 w1*iblock%variable(i,j,k)%res_4
                    t1    =  gm1*t1*t3*t3
                    t2    =  -uban*iblock%variable(i,j,k)%res_1  + &
                                   kx(i,j,k)*iblock%variable(i,j,k)%res_2  + &
                                   ky(i,j,k)*iblock%variable(i,j,k)%res_3  + &
                                   kz(i,j,k)*iblock%variable(i,j,k)%res_4
                    t2    =  t2*t3
                    t3    = -vban*iblock%variable(i,j,k)%res_1 + &
                                  lx(i,j,k)*iblock%variable(i,j,k)%res_2 + &
                                  ly(i,j,k)*iblock%variable(i,j,k)%res_3 + &
                                  lz(i,j,k)*iblock%variable(i,j,k)%res_4
                    !>
                    !>
                    iblock%variable(i,j,k)%res_3 =  -wban*iblock%variable(i,j,k)%res_1 + &
                                                    mx(i,j,k)*iblock%variable(i,j,k)%res_2 + &
                                                    my(i,j,k)*iblock%variable(i,j,k)%res_3 + &
                                                    mz(i,j,k)*iblock%variable(i,j,k)%res_4
                    iblock%variable(i,j,k)%res_1 =  iblock%variable(i,j,k)%res_1+t1
                    iblock%variable(i,j,k)%res_2 =  t3
                    iblock%variable(i,j,k)%res_4 =  0.5e0*(t2-t1)
                    iblock%variable(i,j,k)%res_5 =  iblock%variable(i,j,k)%res_4-t2
                    !>
                    !>
                end do
            end do
        end do
        !>
        !>
        do k=2,kdim
            do i=2,idim
                !>
                !>
                flmp(i,2,k,1)    =0.0
                flmp(i,2,k,2)    =0.0
                flmp(i,2,k,3)    =0.0

                flmn(i,jdim,k,1) =0.0
                flmn(i,jdim,k,2) =0.0
                flmn(i,jdim,k,3) =0.0
                !> Thomas algorithm
                !> perform the scalar tridiagonal LU decomposition
                !>
                diag(i,2,k,1) = diag(i,2,k,1)
                diag(i,2,k,2) = diag(i,2,k,2)
                diag(i,2,k,3) = diag(i,2,k,3)
                flmn(i,2,k,1) = flmn(i,2,k,1) /diag(i,2,k,1)
                flmn(i,2,k,2) = flmn(i,2,k,2) /diag(i,2,k,2)
                flmn(i,2,k,3) = flmn(i,2,k,3) /diag(i,2,k,3)
                do j=3,jdim
                    diag(i,j,k,1)   = diag(i,j,k,1) - flmp(i,j,k,1)*flmn(i,j-1,k,1)
                    diag(i,j,k,2)   = diag(i,j,k,2) - flmp(i,j,k,2)*flmn(i,j-1,k,2)
                    diag(i,j,k,3)   = diag(i,j,k,3) - flmp(i,j,k,3)*flmn(i,j-1,k,3)
                    flmn(i,j,k,1)   = flmn(i,j,k,1)/diag(i,j,k,1)
                    flmn(i,j,k,2)   = flmn(i,j,k,2)/diag(i,j,k,2)
                    flmn(i,j,k,3)   = flmn(i,j,k,3)/diag(i,j,k,3)
                end do
                !>
                !> forward sweep
                !>
                iblock%variable(i,2,k)%res_1  = iblock%variable(i,2,k)%res_1 / diag(i,2,k,1)
                iblock%variable(i,2,k)%res_2  = iblock%variable(i,2,k)%res_2 / diag(i,2,k,1)
                iblock%variable(i,2,k)%res_3  = iblock%variable(i,2,k)%res_3 / diag(i,2,k,1)
                iblock%variable(i,2,k)%res_4  = iblock%variable(i,2,k)%res_4 / diag(i,2,k,2)
                iblock%variable(i,2,k)%res_5  = iblock%variable(i,2,k)%res_5 / diag(i,2,k,3)
                do j=3,jdim
                    iblock%variable(i,j,k)%res_1 = (iblock%variable(i,j,k)%res_1  - flmp(i,j,k,1)* iblock%variable(i,j-1,k)%res_1)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_2 = (iblock%variable(i,j,k)%res_2  - flmp(i,j,k,1)* iblock%variable(i,j-1,k)%res_2)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_3 = (iblock%variable(i,j,k)%res_3  - flmp(i,j,k,1)* iblock%variable(i,j-1,k)%res_3)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_4 = (iblock%variable(i,j,k)%res_4  - flmp(i,j,k,2)* iblock%variable(i,j-1,k)%res_4)/diag(i,j,k,2)
                    iblock%variable(i,j,k)%res_5 = (iblock%variable(i,j,k)%res_5  - flmp(i,j,k,3)* iblock%variable(i,j-1,k)%res_5)/diag(i,j,k,3)
                end do
                !> back sweep
                !>
                iblock%variable(i,jdim,k)%res_1 = iblock%variable(i,jdim,k)%res_1
                iblock%variable(i,jdim,k)%res_2 = iblock%variable(i,jdim,k)%res_2
                iblock%variable(i,jdim,k)%res_3 = iblock%variable(i,jdim,k)%res_3
                iblock%variable(i,jdim,k)%res_4 = iblock%variable(i,jdim,k)%res_4
                iblock%variable(i,jdim,k)%res_5 = iblock%variable(i,jdim,k)%res_5

                do j=jdim-1,2,-1
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1 - flmn(i,j,k,1)*iblock%variable(i,j+1,k)%res_1
                    iblock%variable(i,j,k)%res_2 = iblock%variable(i,j,k)%res_2 - flmn(i,j,k,1)*iblock%variable(i,j+1,k)%res_2
                    iblock%variable(i,j,k)%res_3 = iblock%variable(i,j,k)%res_3 - flmn(i,j,k,1)*iblock%variable(i,j+1,k)%res_3
                    iblock%variable(i,j,k)%res_4 = iblock%variable(i,j,k)%res_4 - flmn(i,j,k,2)*iblock%variable(i,j+1,k)%res_4
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_5 - flmn(i,j,k,3)*iblock%variable(i,j+1,k)%res_5

                end do

            end do
        end do
        !>
        !>
        !> Multiply the inverse of the diagonalizing matrix
        !> T times the change in characteristic combination of variables.
        !> M(inverse)*T*R
        do k=2,kdim
            do i=2,idim
                do j=2,jdim
                    ca   = sqrt(gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r)
                    t1      = 1.0/iblock%cell(i,j,k)%r
                    t2      = t1*iblock%variable(i,j,k)%res_2
                    t3      = t1*iblock%variable(i,j,k)%res_3
                    t5      = t1* ca*(iblock%variable(i,j,k)%res_4-iblock%variable(i,j,k)%res_5)
                    !>
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_4+iblock%variable(i,j,k)%res_5
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1+iblock%variable(i,j,k)%res_5
                    iblock%variable(i,j,k)%res_5 = ca*ca*iblock%variable(i,j,k)%res_5
                    !>
                    iblock%variable(i,j,k)%res_2 = lx(i,j,k)*t2 + mx(i,j,k)*t3 + kx(i,j,k)*t5
                    iblock%variable(i,j,k)%res_3 = ly(i,j,k)*t2 + my(i,j,k)*t3 + ky(i,j,k)*t5
                    iblock%variable(i,j,k)%res_4 = lz(i,j,k)*t2 + mz(i,j,k)*t3 + kz(i,j,k)*t5
                end do
                !>
                !>
            end do
        end do
!       --------------------------------------------------------------
!       k-direction
!       --------------------------------------------------------------
		!>
		!> rhs = rhs*vol/dt
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                    iblock%variable(i,j,k)%res_2 = iblock%variable(i,j,k)%res_2*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                    iblock%variable(i,j,k)%res_3 = iblock%variable(i,j,k)%res_3*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                    iblock%variable(i,j,k)%res_4 = iblock%variable(i,j,k)%res_4*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_5*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                end do
            end do
        end do
        !>
        !>
        do j = 2,jdim
            do i = 2,idim
                do k = 2,kdim

                    ra     = iblock%cell(i,j,k)%r
                    u1     = iblock%cell(i,j,k)%u
                    v1     = iblock%cell(i,j,k)%v
                    w1     = iblock%cell(i,j,k)%w
                    !>
                    !>
                    kx(i,j,k)  = iblock%metric(i,j,k)%hi+iblock%metric(i,j,k+1)%hi
                    ky(i,j,k)  = iblock%metric(i,j,k)%hj+iblock%metric(i,j,k+1)%hj
                    kz(i,j,k)  = iblock%metric(i,j,k)%hk+iblock%metric(i,j,k+1)%hk
                    !>
                    !>
                    xyz = kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k) + kz(i,j,k)*kz(i,j,k)
                    xyz = 1.e0/sqrt(xyz)
                    kx(i,j,k)   = xyz*kx(i,j,k)
                    ky(i,j,k)   = xyz*ky(i,j,k)
                    kz(i,j,k)   = xyz*kz(i,j,k)
                    !>
                    !>
                    lx(i,j,k)  = ky(i,j,k) - kz(i,j,k)
                    ly(i,j,k)  = kz(i,j,k) - kx(i,j,k)
                    lz(i,j,k)  = kx(i,j,k) - ky(i,j,k)
                    xyz = lx(i,j,k)*lx(i,j,k) + ly(i,j,k)*ly(i,j,k) + lz(i,j,k)*lz(i,j,k)
                    xyz = 1.e0/sqrt(xyz)
                    lx(i,j,k) = xyz*lx(i,j,k)
                    ly(i,j,k) = xyz*ly(i,j,k)
                    lz(i,j,k) = xyz*lz(i,j,k)
                    !>
                    !>
                    mx(i,j,k) = ky(i,j,k)*lz(i,j,k) - ly(i,j,k)*kz(i,j,k)
                    my(i,j,k) = kz(i,j,k)*lx(i,j,k) - lz(i,j,k)*kx(i,j,k)
                    mz(i,j,k) = kx(i,j,k)*ly(i,j,k) - lx(i,j,k)*ky(i,j,k)
                    !>
                    !>
                    ca  =sqrt(gamma*iblock%cell(i,j,k)%p/ra)
                    c2  =ca*ca
                    uban =kx(i,j,k)*u1+ky(i,j,k)*v1+kz(i,j,k)*w1


                    flam1(i,j,k)     =uban
                    flam2(i,j,k)     =uban+ca
                    flam3(i,j,k)     =uban-ca
                    abs_flam1(i,j,k) = abs(uban)
                    abs_flam2(i,j,k) = abs(uban+ca)
                    abs_flam3(i,j,k) = abs(uban-ca)
                    !>
                    !>
                    !> used the limiter eigenvalues a la-harten and gnoffo
                    !> implicit tvd schemes for hyperbolic conservation laws in curvilinear coordinates
                    !> the shock-caputring schemes
                    !>
                    !> 0.01  < epsa_r < 0.4
                    !> 0.1   < epsaa  < 0.3     approximate range for blunt bodies
                    !> 0.005 < epsaa  < 0.05    approximate range for slender bodies
                    if(epsa_r .gt. 0)then
                        epsaa = 2.0*epsa_r*(ca + abs(u1) + abs(v1) + abs(w1))
                        epsbb = 0.25/max(epsaa,zero)
                        epscc = 2.00*epsaa
                        if(real(abs_flam1(i,j,k)) .lt. real(epscc)) abs_flam1(i,j,k)= abs_flam1(i,j,k)*abs_flam1(i,j,k)*epsbb+epsaa
                        if(real(abs_flam2(i,j,k)) .lt. real(epscc)) abs_flam2(i,j,k)= abs_flam2(i,j,k)*abs_flam2(i,j,k)*epsbb+epsaa
                        if(real(abs_flam3(i,j,k)) .lt. real(epscc)) abs_flam3(i,j,k)= abs_flam3(i,j,k)*abs_flam3(i,j,k)*epsbb+epsaa
                    end if
                end do
            end do
        end do
        !>
        !>
        do j=2,jdim
            do i=2,idim
                do k=2,kdim+1
                    !>
                    !>
                    !>
                    !>implicit treatment of viscous terms
                    !>
                    if(nvisc == 0)then
                        jo = 0.0
                    elseif( nvisc == 1)then
                        if(k==2 )then
                            jo = iblock%metric(i,j,k)%jacobian/iblock%cell(i,j,k)%r
                        elseif(k == kdim+1)then
                            jo = iblock%metric(i,j,k-1)%jacobian/iblock%cell(i,j,k-1)%r
                        else
                            jo = (1.0/iblock%cell(i,j,k)%r + 1.0/iblock%cell(i,j,k-1)%r)/(iblock%metric(i,j,k)%volume +iblock%metric(i,j,k-1)%volume)
                        end if
                    elseif( nvisc .gt. 1)then
                        if(k==2 )then
                            jo =(1.0+iblock%turbulent(i,j,k)%viscous)*iblock%metric(i,j,k)%jacobian/iblock%cell(i,j,k)%r
                        elseif(k == kdim +1)then
                            jo =(1.0+iblock%turbulent(i,j,k-1)%viscous)*iblock%metric(i,j,k-1)%jacobian/iblock%cell(i,j,k-1)%r
                        else
                            jo =((1.0+iblock%turbulent(i,j,k  )%viscous)/iblock%cell(i,j,k  )%r  + &
                                 (1.0+iblock%turbulent(i,j,k-1)%viscous)/iblock%cell(i,j,k-1)%r  ) / (iblock%metric(i,j,k)%volume+iblock%metric(i,j,k-1)%volume)
                        end if
                    else
                        write(*,*) 'have a wrong turbulent model!please check the mode used and modify!'
                    end if
                    !>
                    !>
                    implicit_viscous_term(i,j,k)      = iblock%metric(i,j,k)%hh*iblock%metric(i,j,k)%hh*jo*avi0

                end do
                do k =3,kdim
                    !>
                    !> low-matrix
                    flmp(i,j,k,1)   = -0.5*iblock%metric(i,j,k)%hh*(flam1(i,j,k-1)+abs_flam1(i,j,k-1))- implicit_viscous_term(i,j,k)
                    flmp(i,j,k,2)   = -0.5*iblock%metric(i,j,k)%hh*(flam2(i,j,k-1)+abs_flam2(i,j,k-1))- implicit_viscous_term(i,j,k)
                    flmp(i,j,k,3)   = -0.5*iblock%metric(i,j,k)%hh*(flam3(i,j,k-1)+abs_flam3(i,j,k-1))- implicit_viscous_term(i,j,k)
                    !> uper-matrix
                    flmn(i,j,k-1,1) = 0.5*iblock%metric(i,j,k)%hh*(flam1(i,j,k)-abs_flam1(i,j,k))- implicit_viscous_term(i,j,k)
                    flmn(i,j,k-1,2) = 0.5*iblock%metric(i,j,k)%hh*(flam2(i,j,k)-abs_flam2(i,j,k))- implicit_viscous_term(i,j,k)
                    flmn(i,j,k-1,3) = 0.5*iblock%metric(i,j,k)%hh*(flam3(i,j,k)-abs_flam3(i,j,k))- implicit_viscous_term(i,j,k)
                end do
                do k=2,kdim
                    !> diag-matrix
                    dt = iblock%cell(i,j,k)%dt
                    diag(i,j,k,1)   = iblock%metric(i,j,k)%volume/dt +  0.5 * iblock%metric(i,j,k+1)%hh * (flam1(i,j,k)+abs_flam1(i,j,k)) - &
                                                                        0.5 * iblock%metric(i,j,k  )%hh * (flam1(i,j,k)-abs_flam1(i,j,k)) + &
                                                                        implicit_viscous_term(i,j,k) + implicit_viscous_term(i,j,k+1)
                    !>
                    !>
                    diag(i,j,k,2)   = iblock%metric(i,j,k)%volume/dt +   0.5 * iblock%metric(i,j,k+1)%hh * (flam2(i,j,k)+abs_flam2(i,j,k)) - &
                                                                         0.5 * iblock%metric(i,j,k  )%hh * (flam2(i,j,k)-abs_flam2(i,j,k)) + &
                                                                         implicit_viscous_term(i,j,k) + implicit_viscous_term(i,j,k+1)
                    !>
                    !>
                    diag(i,j,k,3)   = iblock%metric(i,j,k)%volume/dt +  0.5 * iblock%metric(i,j,k+1)%hh * (flam3(i,j,k)+abs_flam3(i,j,k)) - &
                                                                        0.5 * iblock%metric(i,j,k  )%hh * (flam3(i,j,k)-abs_flam3(i,j,k)) + &
                                                                        implicit_viscous_term(i,j,k) + implicit_viscous_term(i,j,k+1)
                    !>
                    !>

                end do
            end do
        end do
        !>
        !> Multiply the inverse of the diagonalizing matrix
        !> T times the residual contribution.
        !> T(inverse)*M*R
        do j=2,jdim
            do i=2,idim
                do k=2,kdim
                    ca    = sqrt(gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r)
                    t1    =  1.e0/ca
                    iblock%variable(i,j,k)%res_5 =  iblock%variable(i,j,k)%res_5*t1*t1
                    iblock%variable(i,j,k)%res_1 =  iblock%variable(i,j,k)%res_1-iblock%variable(i,j,k)%res_5
                    t2    =  iblock%cell(i,j,k)%r*iblock%variable(i,j,k)%res_2
                    t3    =  iblock%cell(i,j,k)%r*iblock%variable(i,j,k)%res_3
                    t4    =  iblock%cell(i,j,k)%r*iblock%variable(i,j,k)%res_4
                    iblock%variable(i,j,k)%res_2 =  lx(i,j,k)*t2+ly(i,j,k)*t3+lz(i,j,k)*t4
                    iblock%variable(i,j,k)%res_3 =  mx(i,j,k)*t2+my(i,j,k)*t3+mz(i,j,k)*t4
                    iblock%variable(i,j,k)%res_4 =  0.5*(t1*(kx(i,j,k)*t2+ky(i,j,k)*t3+kz(i,j,k)*t4)+iblock%variable(i,j,k)%res_5)
                    iblock%variable(i,j,k)%res_5 =  -iblock%variable(i,j,k)%res_4+iblock%variable(i,j,k)%res_5

                end do
            end do
        end do
        !>
        !>
        !>
        do j=2,jdim
            do i=2,idim
                !>
                !>
                flmp(i,j,2,1)    =0.0
                flmp(i,j,2,2)    =0.0
                flmp(i,j,2,3)    =0.0

                flmn(i,j,kdim,1) =0.0
                flmn(i,j,kdim,2) =0.0
                flmn(i,j,kdim,3) =0.0

                !> thomas algorithm
                !> and the solution is exact
                !>perform the scalar tridiagonal lu decomposition
                !>
                diag(i,j,2,1) = diag(i,j,2,1)
                diag(i,j,2,2) = diag(i,j,2,2)
                diag(i,j,2,3) = diag(i,j,2,3)
                flmn(i,j,2,1) = flmn(i,j,2,1) /diag(i,j,2,1)
                flmn(i,j,2,2) = flmn(i,j,2,2) /diag(i,j,2,2)
                flmn(i,j,2,3) = flmn(i,j,2,3) /diag(i,j,2,3)
                do k=3,kdim
                    diag(i,j,k,1)   = diag(i,j,k,1) - flmp(i,j,k,1)*flmn(i,j,k-1,1)
                    diag(i,j,k,2)   = diag(i,j,k,2) - flmp(i,j,k,2)*flmn(i,j,k-1,2)
                    diag(i,j,k,3)   = diag(i,j,k,3) - flmp(i,j,k,3)*flmn(i,j,k-1,3)
                    flmn(i,j,k,1)   = flmn(i,j,k,1)/diag(i,j,k,1)
                    flmn(i,j,k,2)   = flmn(i,j,k,2)/diag(i,j,k,2)
                    flmn(i,j,k,3)   = flmn(i,j,k,3)/diag(i,j,k,3)
                end do
                !> forward sweep
                !>
                iblock%variable(i,j,2)%res_1  = iblock%variable(i,j,2)%res_1 / diag(i,j,2,1)
                iblock%variable(i,j,2)%res_2  = iblock%variable(i,j,2)%res_2 / diag(i,j,2,1)
                iblock%variable(i,j,2)%res_3  = iblock%variable(i,j,2)%res_3 / diag(i,j,2,1)
                iblock%variable(i,j,2)%res_4  = iblock%variable(i,j,2)%res_4 / diag(i,j,2,2)
                iblock%variable(i,j,2)%res_5  = iblock%variable(i,j,2)%res_5 / diag(i,j,2,3)

                do k=3,kdim
                    iblock%variable(i,j,k)%res_1 = (iblock%variable(i,j,k)%res_1 - flmp(i,j,k,1)*iblock%variable(i,j,k-1)%res_1)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_2 = (iblock%variable(i,j,k)%res_2 - flmp(i,j,k,1)*iblock%variable(i,j,k-1)%res_2)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_3 = (iblock%variable(i,j,k)%res_3 - flmp(i,j,k,1)*iblock%variable(i,j,k-1)%res_3)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_4 = (iblock%variable(i,j,k)%res_4 - flmp(i,j,k,2)*iblock%variable(i,j,k-1)%res_4)/diag(i,j,k,2)
                    iblock%variable(i,j,k)%res_5 = (iblock%variable(i,j,k)%res_5 - flmp(i,j,k,3)*iblock%variable(i,j,k-1)%res_5)/diag(i,j,k,3)

                end do
                !> backward sweep
                !>
                iblock%variable(i,j,kdim)%res_1 = iblock%variable(i,j,kdim)%res_1
                iblock%variable(i,j,kdim)%res_2 = iblock%variable(i,j,kdim)%res_2
                iblock%variable(i,j,kdim)%res_3 = iblock%variable(i,j,kdim)%res_3
                iblock%variable(i,j,kdim)%res_4 = iblock%variable(i,j,kdim)%res_4
                iblock%variable(i,j,kdim)%res_5 = iblock%variable(i,j,kdim)%res_5

                do k=kdim-1,2,-1
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1 - flmn(i,j,k,1)*iblock%variable(i,j,k+1)%res_1
                    iblock%variable(i,j,k)%res_2 = iblock%variable(i,j,k)%res_2 - flmn(i,j,k,1)*iblock%variable(i,j,k+1)%res_2
                    iblock%variable(i,j,k)%res_3 = iblock%variable(i,j,k)%res_3 - flmn(i,j,k,1)*iblock%variable(i,j,k+1)%res_3
                    iblock%variable(i,j,k)%res_4 = iblock%variable(i,j,k)%res_4 - flmn(i,j,k,2)*iblock%variable(i,j,k+1)%res_4
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_5 - flmn(i,j,k,3)*iblock%variable(i,j,k+1)%res_5
                end do
            end do
        end do
        !>
        !>
        !> Multiply the inverse of the diagonalizing matrix
        !> T times the change in characteristic combination of variables.
        !> M(inverse)*T*R
        do j=2,jdim
            do i=2,idim
                do k=2,kdim
                    ca   = sqrt(gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r)
                    t1      = 1.0/iblock%cell(i,j,k)%r
                    t2      = t1*iblock%variable(i,j,k)%res_2
                    t3      = t1*iblock%variable(i,j,k)%res_3
                    t5      = t1* ca*(iblock%variable(i,j,k)%res_4-iblock%variable(i,j,k)%res_5)
                    !>
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_4+iblock%variable(i,j,k)%res_5
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1+iblock%variable(i,j,k)%res_5
                    iblock%variable(i,j,k)%res_5 = ca*ca*iblock%variable(i,j,k)%res_5
                    !>
                    iblock%variable(i,j,k)%res_2 = lx(i,j,k)*t2 + mx(i,j,k)*t3 + kx(i,j,k)*t5
                    iblock%variable(i,j,k)%res_3 = ly(i,j,k)*t2 + my(i,j,k)*t3 + ky(i,j,k)*t5
                    iblock%variable(i,j,k)%res_4 = lz(i,j,k)*t2 + mz(i,j,k)*t3 + kz(i,j,k)*t5
                end do
            end do
        end do
!       --------------------------------------------------------------
!       i-direction
!       --------------------------------------------------------------
		!>
		!> rhs = rhs*vol/dt
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                    iblock%variable(i,j,k)%res_2 = iblock%variable(i,j,k)%res_2*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                    iblock%variable(i,j,k)%res_3 = iblock%variable(i,j,k)%res_3*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                    iblock%variable(i,j,k)%res_4 = iblock%variable(i,j,k)%res_4*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_5*iblock%metric(i,j,k)%volume/iblock%cell(i,j,k)%dt
                end do
            end do
        end do
        do k = 2,kdim
            do j = 2,jdim
                do i = 2,idim

                    ra     = iblock%cell(i,j,k)%r
                    u1     = iblock%cell(i,j,k)%u
                    v1     = iblock%cell(i,j,k)%v
                    w1     = iblock%cell(i,j,k)%w
                    !>
                    !>
                    kx(i,j,k)  = iblock%metric(i,j,k)%fi+iblock%metric(i+1,j,k)%fi
                    ky(i,j,k)  = iblock%metric(i,j,k)%fj+iblock%metric(i+1,j,k)%fj
                    kz(i,j,k)  = iblock%metric(i,j,k)%fk+iblock%metric(i+1,j,k)%fk
                    !>
                    !>
                    xyz = kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k) + kz(i,j,k)*kz(i,j,k)
                    xyz = 1.e0/sqrt(xyz)
                    kx(i,j,k)   = xyz*kx(i,j,k)
                    ky(i,j,k)   = xyz*ky(i,j,k)
                    kz(i,j,k)   = xyz*kz(i,j,k)
                    !>
                    !>
                    lx(i,j,k)  = ky(i,j,k) - kz(i,j,k)
                    ly(i,j,k)  = kz(i,j,k) - kx(i,j,k)
                    lz(i,j,k)  = kx(i,j,k) - ky(i,j,k)
                    xyz    = lx(i,j,k)*lx(i,j,k) + ly(i,j,k)*ly(i,j,k) + lz(i,j,k)*lz(i,j,k)
                    xyz    = 1.e0/sqrt(xyz)
                    lx(i,j,k) = xyz*lx(i,j,k)
                    ly(i,j,k) = xyz*ly(i,j,k)
                    lz(i,j,k) = xyz*lz(i,j,k)
                    !>
                    !>
                    mx(i,j,k) = ky(i,j,k)*lz(i,j,k) - ly(i,j,k)*kz(i,j,k)
                    my(i,j,k) = kz(i,j,k)*lx(i,j,k) - lz(i,j,k)*kx(i,j,k)
                    mz(i,j,k) = kx(i,j,k)*ly(i,j,k) - lx(i,j,k)*ky(i,j,k)
                    !>
                    !>
                    ca  =sqrt(gamma*iblock%cell(i,j,k)%p/ra)
                    c2  =ca*ca
                    uban =kx(i,j,k)*u1+ky(i,j,k)*v1+kz(i,j,k)*w1

                    flam1(i,j,k)     = uban
                    flam2(i,j,k)     = uban+ca
                    flam3(i,j,k)     = uban-ca
                    abs_flam1(i,j,k) = abs(uban)
                    abs_flam2(i,j,k) = abs(uban+ca)
                    abs_flam3(i,j,k) = abs(uban-ca)
                    !>
                    !>
                    !> used the limiter eigenvalues a La-Harten and Gnoffo
                    !> implicit tvd schemes for hyperbolic conservation laws in curvilinear coordinates
                    !> the shock-caputring schemes
                    !>
                    !> 0.01  < epsa_r < 0.4
                    !> 0.1   < epsaa  < 0.3     approximate range for blunt bodies
                    !> 0.005 < epsaa  < 0.05    approximate range for slender bodies
                    if(epsa_r .gt. 0)then
                        epsaa = 2.0*epsa_r*(ca + abs(u1) + abs(v1) + abs(w1))
                        epsbb = 0.25/max(epsaa,zero)
                        epscc = 2.0*epsaa
                        if(real(abs_flam1(i,j,k)) .lt. real(epscc)) abs_flam1(i,j,k)= abs_flam1(i,j,k)*abs_flam1(i,j,k)*epsbb+epsaa
                        if(real(abs_flam2(i,j,k)) .lt. real(epscc)) abs_flam2(i,j,k)= abs_flam2(i,j,k)*abs_flam2(i,j,k)*epsbb+epsaa
                        if(real(abs_flam3(i,j,k)) .lt. real(epscc)) abs_flam3(i,j,k)= abs_flam3(i,j,k)*abs_flam3(i,j,k)*epsbb+epsaa
                    end if
                end do
            end do
        end do
        !>
        !>
        do k=2,kdim
            do j=2,jdim
                do i=2,idim+1
                    !>
                    !>
                    !>
                    !>implicit treatment of viscous term
                    !>
                    if(nvisc == 0)then
                        jo = 0.0
                    elseif( nvisc == 1)then
                        if(i==2 )then
                            jo = iblock%metric(i,j,k)%jacobian/iblock%cell(i,j,k)%r
                        elseif(i == idim+1)then
                            jo = iblock%metric(i-1,j,k)%jacobian/iblock%cell(i-1,j,k)%r
                        else
                            jo = (1.0/iblock%cell(i,j,k)%r + 1.0/iblock%cell(i-1,j,k)%r)/(iblock%metric(i,j,k)%volume+iblock%metric(i-1,j,k)%volume)
                        end if
                    elseif( nvisc .gt. 1)then
                        if(i==2 )then
                            jo =(1.0+iblock%turbulent(i,j,k)%viscous)*iblock%metric(i,j,k)%jacobian/iblock%cell(i,j,k)%r
                        elseif(i == idim +1)then
                            jo =(1.0+iblock%turbulent(i-1,j,k)%viscous)*iblock%metric(i-1,j,k)%jacobian/iblock%cell(i-1,j,k)%r
                        else
                            jo =((1.0+iblock%turbulent(i,j,k)%viscous)/iblock%cell(i,j,k)%r+(1.0+iblock%turbulent(i-1,j,k)%viscous)/iblock%cell(i-1,j,k)%r)/(iblock%metric(i,j,k)%volume+iblock%metric(i-1,j,k)%volume)
                        end if
                    else
                        write(*,*) 'have a wrong turbulent model!please check the mode used and modify!'
                    end if
                    !>
                    !>
                    implicit_viscous_term(i,j,k) = iblock%metric(i,j,k)%ff*iblock%metric(i,j,k)%ff*jo*avi0
                    !>
                end do
                do i =3,idim
                    !>
                    !> low-matrix
                    flmp(i,j,k,1) = -0.5 * iblock%metric(i,j,k)%ff * (flam1(i-1,j,k)+abs_flam1(i-1,j,k))- implicit_viscous_term(i,j,k)
                    flmp(i,j,k,2) = -0.5 * iblock%metric(i,j,k)%ff * (flam2(i-1,j,k)+abs_flam2(i-1,j,k))- implicit_viscous_term(i,j,k)
                    flmp(i,j,k,3) = -0.5 * iblock%metric(i,j,k)%ff * (flam3(i-1,j,k)+abs_flam3(i-1,j,k))- implicit_viscous_term(i,j,k)
                    !> uper-matrix
                    flmn(i-1,j,k,1) = 0.5 * iblock%metric(i,j,k)%ff * (flam1(i,j,k)-abs_flam1(i,j,k))- implicit_viscous_term(i,j,k)
                    flmn(i-1,j,k,2) = 0.5 * iblock%metric(i,j,k)%ff * (flam2(i,j,k)-abs_flam2(i,j,k))- implicit_viscous_term(i,j,k)
                    flmn(i-1,j,k,3) = 0.5 * iblock%metric(i,j,k)%ff * (flam3(i,j,k)-abs_flam3(i,j,k))- implicit_viscous_term(i,j,k)
                end do
                do i=2,idim
                    !> diag-matrix
                    dt = iblock%cell(i,j,k)%dt
                    diag(i,j,k,1)   = iblock%metric(i,j,k)%volume/dt + 0.5 * iblock%metric(i+1,j,k)%ff * (flam1(i,j,k)+abs_flam1(i,j,k)) - &
                                                                       0.5 * iblock%metric(i  ,j,k)%ff * (flam1(i,j,k)-abs_flam1(i,j,k)) + &
                                                                       implicit_viscous_term(i,j,k) + implicit_viscous_term(i+1,j,k)
                    !>
                    !>
                    diag(i,j,k,2)   = iblock%metric(i,j,k)%volume/dt + 0.5 * iblock%metric(i+1,j,k)%ff * (flam2(i,j,k)+abs_flam2(i,j,k)) - &
                                                                       0.5 * iblock%metric(i  ,j,k)%ff * (flam2(i,j,k)-abs_flam2(i,j,k)) + &
                                                                       implicit_viscous_term(i,j,k) + implicit_viscous_term(i+1,j,k)
                    !>
                    !>
                    diag(i,j,k,3)   = iblock%metric(i,j,k)%volume/dt + 0.5 * iblock%metric(i+1,j,k)%ff * (flam3(i,j,k)+abs_flam3(i,j,k)) - &
                                                                       0.5 * iblock%metric(i  ,j,k)%ff * (flam3(i,j,k)-abs_flam3(i,j,k)) + &
                                                                       implicit_viscous_term(i,j,k) + implicit_viscous_term(i+1,j,k)
                end do
            end do
        end do
        !>
        !>
        !> Multiply the inverse of the diagonalizing matrix
        !> T times the residual contribution.
        !> T(inverse)*M*R
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    ca    = sqrt(gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r)
                    t1    =  1.e0/ca
                    iblock%variable(i,j,k)%res_5 =  iblock%variable(i,j,k)%res_5*t1*t1
                    iblock%variable(i,j,k)%res_1 =  iblock%variable(i,j,k)%res_1-iblock%variable(i,j,k)%res_5
                    t2    =  iblock%cell(i,j,k)%r*iblock%variable(i,j,k)%res_2
                    t3    =  iblock%cell(i,j,k)%r*iblock%variable(i,j,k)%res_3
                    t4    =  iblock%cell(i,j,k)%r*iblock%variable(i,j,k)%res_4
                    iblock%variable(i,j,k)%res_2 =  lx(i,j,k)*t2+ly(i,j,k)*t3+lz(i,j,k)*t4
                    iblock%variable(i,j,k)%res_3 =  mx(i,j,k)*t2+my(i,j,k)*t3+mz(i,j,k)*t4
                    iblock%variable(i,j,k)%res_4 =  0.5*(t1*(kx(i,j,k)*t2+ky(i,j,k)*t3+kz(i,j,k)*t4)+iblock%variable(i,j,k)%res_5)
                    iblock%variable(i,j,k)%res_5 =  -iblock%variable(i,j,k)%res_4+iblock%variable(i,j,k)%res_5
                end do
            end do
        end do
        !>
        do k=2,kdim
            do j=2,jdim
                flmp(2,j,k,1)    =0.0
                flmp(2,j,k,2)    =0.0
                flmp(2,j,k,3)    =0.0

                flmn(idim,j,k,1) =0.0
                flmn(idim,j,k,2) =0.0
                flmn(idim,j,k,3) =0.0
                !> Thomas algorithm
                !> perform the scalar tridiagonal LU decomposition
                !>
                diag(2,j,k,1) = diag(2,j,k,1)
                diag(2,j,k,2) = diag(2,j,k,2)
                diag(2,j,k,3) = diag(2,j,k,3)
                flmn(2,j,k,1) = flmn(2,j,k,1) /diag(2,j,k,1)
                flmn(2,j,k,2) = flmn(2,j,k,2) /diag(2,j,k,2)
                flmn(2,j,k,3) = flmn(2,j,k,3) /diag(2,j,k,3)
                do i=3,idim
                    diag(i,j,k,1)   = diag(i,j,k,1) - flmp(i,j,k,1)*flmn(i-1,j,k,1)
                    diag(i,j,k,2)   = diag(i,j,k,2) - flmp(i,j,k,2)*flmn(i-1,j,k,2)
                    diag(i,j,k,3)   = diag(i,j,k,3) - flmp(i,j,k,3)*flmn(i-1,j,k,3)
                    flmn(i,j,k,1)   = flmn(i,j,k,1)/diag(i,j,k,1)
                    flmn(i,j,k,2)   = flmn(i,j,k,2)/diag(i,j,k,2)
                    flmn(i,j,k,3)   = flmn(i,j,k,3)/diag(i,j,k,3)
                end do
                !> forward sweep
                !>
                iblock%variable(2,j,k)%res_1 = iblock%variable(2,j,k)%res_1 / diag(2,j,k,1)
                iblock%variable(2,j,k)%res_2 = iblock%variable(2,j,k)%res_2 / diag(2,j,k,1)
                iblock%variable(2,j,k)%res_3 = iblock%variable(2,j,k)%res_3 / diag(2,j,k,1)
                iblock%variable(2,j,k)%res_4 = iblock%variable(2,j,k)%res_4 / diag(2,j,k,2)
                iblock%variable(2,j,k)%res_5 = iblock%variable(2,j,k)%res_5 / diag(2,j,k,3)

                do i=3,idim
                    iblock%variable(i,j,k)%res_1 = (iblock%variable(i,j,k)%res_1 - flmp(i,j,k,1)*iblock%variable(i-1,j,k)%res_1)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_2 = (iblock%variable(i,j,k)%res_2 - flmp(i,j,k,1)*iblock%variable(i-1,j,k)%res_2)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_3 = (iblock%variable(i,j,k)%res_3 - flmp(i,j,k,1)*iblock%variable(i-1,j,k)%res_3)/diag(i,j,k,1)
                    iblock%variable(i,j,k)%res_4 = (iblock%variable(i,j,k)%res_4 - flmp(i,j,k,2)*iblock%variable(i-1,j,k)%res_4)/diag(i,j,k,2)
                    iblock%variable(i,j,k)%res_5 = (iblock%variable(i,j,k)%res_5 - flmp(i,j,k,3)*iblock%variable(i-1,j,k)%res_5)/diag(i,j,k,3)
                end do
                !> backward sweep
                !>
                iblock%variable(idim,j,k)%res_1 = iblock%variable(idim,j,k)%res_1
                iblock%variable(idim,j,k)%res_2 = iblock%variable(idim,j,k)%res_2
                iblock%variable(idim,j,k)%res_3 = iblock%variable(idim,j,k)%res_3
                iblock%variable(idim,j,k)%res_4 = iblock%variable(idim,j,k)%res_4
                iblock%variable(idim,j,k)%res_5 = iblock%variable(idim,j,k)%res_5

                do i=idim-1,2,-1
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1 - flmn(i,j,k,1)*iblock%variable(i+1,j,k)%res_1
                    iblock%variable(i,j,k)%res_2 = iblock%variable(i,j,k)%res_2 - flmn(i,j,k,1)*iblock%variable(i+1,j,k)%res_2
                    iblock%variable(i,j,k)%res_3 = iblock%variable(i,j,k)%res_3 - flmn(i,j,k,1)*iblock%variable(i+1,j,k)%res_3
                    iblock%variable(i,j,k)%res_4 = iblock%variable(i,j,k)%res_4 - flmn(i,j,k,2)*iblock%variable(i+1,j,k)%res_4
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_5 - flmn(i,j,k,3)*iblock%variable(i+1,j,k)%res_5
                    !>
                    !>
                end do
            end do
        end do
        !>
        !>
        !> Multiply the inverse of the diagonalizing matrix
        !> T times the change in characteristic combination of variables.
        !> M(inverse)*T*R
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    ca   = sqrt(gamma*iblock%cell(i,j,k)%p/iblock%cell(i,j,k)%r)
                    t1      = 1.0/iblock%cell(i,j,k)%r
                    t2      = t1*iblock%variable(i,j,k)%res_2
                    t3      = t1*iblock%variable(i,j,k)%res_3
                    t5      = t1* ca*(iblock%variable(i,j,k)%res_4-iblock%variable(i,j,k)%res_5)
                    !>
                    iblock%variable(i,j,k)%res_5 = iblock%variable(i,j,k)%res_4+iblock%variable(i,j,k)%res_5
                    iblock%variable(i,j,k)%res_1 = iblock%variable(i,j,k)%res_1+iblock%variable(i,j,k)%res_5
                    iblock%variable(i,j,k)%res_5 = ca*ca*iblock%variable(i,j,k)%res_5
                    !>
                    iblock%variable(i,j,k)%res_2 = lx(i,j,k)*t2 + mx(i,j,k)*t3 + kx(i,j,k)*t5
                    iblock%variable(i,j,k)%res_3 = ly(i,j,k)*t2 + my(i,j,k)*t3 + ky(i,j,k)*t5
                    iblock%variable(i,j,k)%res_4 = lz(i,j,k)*t2 + mz(i,j,k)*t3 + kz(i,j,k)*t5
                end do
            end do
        end do
		!>
		!>
        !>deallocate the dimension allocate
        deallocate(flmn)
        deallocate(flmp)
        deallocate(diag)
        deallocate(flam1)
        deallocate(flam2)
        deallocate(flam3)
        deallocate(abs_flam1)
        deallocate(abs_flam2)
        deallocate(abs_flam3)
        deallocate(implicit_viscous_term)
        deallocate(kx,ky,kz)
        deallocate(lx,ly,lz)
        deallocate(mx,my,mz)


        return
    end subroutine



