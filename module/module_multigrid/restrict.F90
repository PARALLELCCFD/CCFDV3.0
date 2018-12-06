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
    !> \brief
    !> purpose:restrict q(the primative variables )withe a volume
    !> weighted interpolation and residual to coarse meshes.also
    !> resrrict turbulent eddy viscosty in the case of turbulent flows
    !> to coarser meshes!
    !!
    subroutine restrict(nbl)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use turbulent_module
        implicit none
        !>
        !>
        integer i,j,k
        integer ii,jj,kk
        integer nbl,nbl_c
        integer idim,jdim,kdim,idim_c,jdim_c,kdim_c
        real(kind = dprec) :: qq,re
        real(kind = dprec),dimension(:,:,:,:),allocatable::q_v
        real(kind = dprec),dimension(:,:,:,:),allocatable::res_v
        !>
        !>
        type(overlap_type),pointer :: mesh
        type(blocks_type),pointer  :: iblock,iblock_c
        !>
        !>
        nbl_c   = nbl+1
        !>
        !>
        mesh     => grids(imesh)
        iblock   => mesh%blocks(nbl)
        iblock_c => mesh%blocks(nbl_c)
        !> nbl_c is coarser mesh blocks indices
        !> --------------------------------------------
        !>      v-cycle   nblocks  = ngrid*global_level
        !>     1       1    ---level 3 (global_level)
        !>       2   2      ---level 2
        !>         3        ---level 1
        !>   i(n,n-1)q(n) to n-1 level mesh
        !>   i(n,n-1)residual(n) to n-1 level mesh
        !>   the detail of methods in  the user manual files
        !>   by liuxz edit
        !> --------------------------------------------
        !> nbl is finer mesh blocks indices
        !> idim,jdim,kdim is finer mesh dimension indices
        !> idim_c,jdim_c,kdim_c is coarser meshs dimension indices
        idim    = iblock%idim
        jdim    = iblock%jdim
        kdim    = iblock%kdim
        !>
        !> the coordinate of the coarser meshes
        !>
        idim_c    = iblock_c%idim
        jdim_c    = iblock_c%jdim
        kdim_c    = iblock_c%kdim

        allocate(q_v(idim,jdim,kdim,9))
        allocate(res_v(idim,jdim,kdim,5))

        !> init the primitive variable
        do k=2,kdim_c!ii1
            do j=2,jdim_c!jj2
                do i=2,idim_c!kk2
                    iblock_c%variable(i,j,k)%res_1  = 0.0
                    iblock_c%variable(i,j,k)%res_2  = 0.0
                    iblock_c%variable(i,j,k)%res_3  = 0.0
                    iblock_c%variable(i,j,k)%res_4  = 0.0
                    iblock_c%variable(i,j,k)%res_5  = 0.0
                end do
            end do
        end do
        !> \brief
        !> the primitive variable with a volume-weighted interpolation
        !> and residual to coarser mesh.
        !> if turbulent module is useful is also restrict eddy viscosity in
        !> the case of turbulent flows to coarser meshs
        !!
        do k=2,kdim
            do j=2,jdim
                do i=2,idim
                    !>rho*vol  volume weighting of rho
                    q_v(i,j,k,1)  =  iblock%cell(i,j,k)%r * iblock%metric(i,j,k)%volume
                    !>u*vol of volume weighting of u
                    q_v(i,j,k,2)  =  iblock%cell(i,j,k)%u * iblock%metric(i,j,k)%volume
                    !>v*vol of volume weighting of v
                    q_v(i,j,k,3)  =  iblock%cell(i,j,k)%v * iblock%metric(i,j,k)%volume
                    !>w*vol of volume weighting of w
                    q_v(i,j,k,4)  =  iblock%cell(i,j,k)%w * iblock%metric(i,j,k)%volume
                    !>p*vol of volume weighting of pressuer
                    q_v(i,j,k,5)  =  iblock%cell(i,j,k)%p * iblock%metric(i,j,k)%volume
                end do
            end do
        end do
        !> \brief
        !> operating the two or three dimension meshs  by volume-weighted methods
        !!
        if(idim .gt. 2)then
            !> three-dimension mesh operating
            !> at i-j,j-k,i-k's face the r u v w p will be update by the idim+3-1's
            !> value,so do scale is kdim+2
            kk = 1
            do k=2,kdim,2
                kk = kk + 1
                jj = 1
                do j=2,jdim,2
                    jj = jj + 1
                    ii = 1
                    do i=2,idim,2
                        ii = ii + 1
                        !> restric rho
                        iblock_c%cell(ii,jj,kk)%r = ( q_v(i  ,j  ,k  ,1)   + &
                                                      q_v(i+1,j  ,k  ,1)   + &
                                                      q_v(i  ,j+1,k  ,1)   + &
                                                      q_v(i  ,j  ,k+1,1)   + &
                                                      q_v(i  ,j+1,k+1,1)   + &
                                                      q_v(i+1,j+1,k  ,1)   + &
                                                      q_v(i+1,j  ,k+1,1)   + &
                                                      q_v(i+1,j+1,k+1,1) ) * iblock_c%metric(ii,jj,kk)%jacobian
                        !> restric u
                        iblock_c%cell(ii,jj,kk)%u = ( q_v(i  ,j  ,k  ,2)   + &
                                                      q_v(i+1,j  ,k  ,2)   + &
                                                      q_v(i  ,j+1,k  ,2)   + &
                                                      q_v(i  ,j  ,k+1,2)   + &
                                                      q_v(i  ,j+1,k+1,2)   + &
                                                      q_v(i+1,j+1,k  ,2)   + &
                                                      q_v(i+1,j  ,k+1,2)   + &
                                                      q_v(i+1,j+1,k+1,2) ) * iblock_c%metric(ii,jj,kk)%jacobian !/ iblock_c%cell(ii,jj,kk)%r
                        !> restric v
                        iblock_c%cell(ii,jj,kk)%v = ( q_v(i  ,j  ,k  ,3)   + &
                                                      q_v(i+1,j  ,k  ,3)   + &
                                                      q_v(i  ,j+1,k  ,3)   + &
                                                      q_v(i  ,j  ,k+1,3)   + &
                                                      q_v(i  ,j+1,k+1,3)   + &
                                                      q_v(i+1,j+1,k  ,3)   + &
                                                      q_v(i+1,j  ,k+1,3)   + &
                                                      q_v(i+1,j+1,k+1,3) ) * iblock_c%metric(ii,jj,kk)%jacobian !/ iblock_c%cell(ii,jj,kk)%r
                        !> restric w
                        iblock_c%cell(ii,jj,kk)%w = ( q_v(i  ,j  ,k  ,4)   + &
                                                      q_v(i+1,j  ,k  ,4)   + &
                                                      q_v(i  ,j+1,k  ,4)   + &
                                                      q_v(i  ,j  ,k+1,4)   + &
                                                      q_v(i  ,j+1,k+1,4)   + &
                                                      q_v(i+1,j+1,k  ,4)   + &
                                                      q_v(i+1,j  ,k+1,4)   + &
                                                      q_v(i+1,j+1,k+1,4) ) * iblock_c%metric(ii,jj,kk)%jacobian !/ iblock_c%cell(ii,jj,kk)%r
                        !> restric the pressure
                        iblock_c%cell(ii,jj,kk)%p = (  q_v(i  ,j  ,k  ,5)   + &
                                                       q_v(i+1,j  ,k  ,5)   + &
                                                       q_v(i  ,j+1,k  ,5)   + &
                                                       q_v(i  ,j  ,k+1,5)   + &
                                                       q_v(i  ,j+1,k+1,5)   + &
                                                       q_v(i+1,j+1,k  ,5)   + &
                                                       q_v(i+1,j  ,k+1,5)   + &
                                                       q_v(i+1,j+1,k+1,5) ) * iblock_c%metric(ii,jj,kk)%jacobian
                        !>
!                        write(121,'("axis=",3i4,6e24.16)') ii,jj,kk,iblock_c%cell(ii,jj,kk)%r ,&
!                        iblock_c%cell(ii,jj,kk)%u ,iblock_c%cell(ii,jj,kk)%v ,&
!                        iblock_c%cell(ii,jj,kk)%w ,iblock_c%cell(ii,jj,kk)%p,iblock_c%metric(ii,jj,kk)%jacobian
                        !> \brief
                        if(nvisc .gt. 1 )then
                            !> restrict the turbulent viscous from the finer mesh to coarse mesh
                            !>
                            iblock_c%turbulent(ii,jj,kk)%viscous = 0.125*(iblock%turbulent(i  ,j  ,k  )%viscous + &
                                                                          iblock%turbulent(i+1,j  ,k  )%viscous + &
                                                                          iblock%turbulent(i  ,j+1,k  )%viscous + &
                                                                          iblock%turbulent(i  ,j  ,k+1)%viscous + &
                                                                          iblock%turbulent(i+1,j+1,k  )%viscous + &
                                                                          iblock%turbulent(i+1,j  ,k+1)%viscous + &
                                                                          iblock%turbulent(i  ,j+1,k+1)%viscous + &
                                                                          iblock%turbulent(i+1,j+1,k+1)%viscous)
                        endif
                        !>
                        !> the parameter for the SA model
                        !>
!                        if(nvisc .ge. 3 )then
!                            !> restrict the turbulent viscous from the finer mesh to coarse mesh
!                            !>
!                            iblock_c%turbulent(ii,jj,kk)%tur_save= 0.125*(iblock%turbulent(i  ,j  ,k  )%tur_save + &
!                                                                          iblock%turbulent(i+1,j  ,k  )%tur_save + &
!                                                                          iblock%turbulent(i  ,j+1,k  )%tur_save + &
!                                                                          iblock%turbulent(i  ,j  ,k+1)%tur_save + &
!                                                                          iblock%turbulent(i+1,j+1,k  )%tur_save + &
!                                                                          iblock%turbulent(i+1,j  ,k+1)%tur_save + &
!                                                                          iblock%turbulent(i  ,j+1,k+1)%tur_save + &
!                                                                          iblock%turbulent(i+1,j+1,k+1)%tur_save)
!                        endif
                        !>
                        !>
                    end do
                end do
            end do
            !> restrict the rn run rvn rwn ren to the coarser mesh from finer mesh
            !>
            !>
            kk = 1
            do k=2,kdim,2
                kk = kk + 1
                jj =1
                do j=2,jdim,2
                    jj = jj + 1
                    ii = 1
                    do i=2,idim,2
                        ii = ii + 1
                        !> restric rhs[1]
                        iblock_c%variable_mg(ii,jj,kk)%res1_mg =    iblock%variable(i  ,j  ,k  )%res_1  + &
                                                                    iblock%variable(i+1,j  ,k  )%res_1  + &
                                                                    iblock%variable(i  ,j+1,k  )%res_1  + &
                                                                    iblock%variable(i  ,j  ,k+1)%res_1  + &
                                                                    iblock%variable(i  ,j+1,k+1)%res_1  + &
                                                                    iblock%variable(i+1,j+1,k  )%res_1  + &
                                                                    iblock%variable(i+1,j  ,k+1)%res_1  + &
                                                                    iblock%variable(i+1,j+1,k+1)%res_1

                        !> restric rhs[2]
                        iblock_c%variable_mg(ii,jj,kk)%res2_mg =    iblock%variable(i  ,j  ,k  )%res_2  + &
                                                                    iblock%variable(i+1,j  ,k  )%res_2  + &
                                                                    iblock%variable(i  ,j+1,k  )%res_2  + &
                                                                    iblock%variable(i  ,j  ,k+1)%res_2  + &
                                                                    iblock%variable(i  ,j+1,k+1)%res_2  + &
                                                                    iblock%variable(i+1,j+1,k  )%res_2  + &
                                                                    iblock%variable(i+1,j  ,k+1)%res_2  + &
                                                                    iblock%variable(i+1,j+1,k+1)%res_2
                        !> restric rhs[3]
                        iblock_c%variable_mg(ii,jj,kk)%res3_mg =    iblock%variable(i  ,j  ,k  )%res_3  + &
                                                                    iblock%variable(i+1,j  ,k  )%res_3  + &
                                                                    iblock%variable(i  ,j+1,k  )%res_3  + &
                                                                    iblock%variable(i  ,j  ,k+1)%res_3  + &
                                                                    iblock%variable(i  ,j+1,k+1)%res_3  + &
                                                                    iblock%variable(i+1,j+1,k  )%res_3  + &
                                                                    iblock%variable(i+1,j  ,k+1)%res_3  + &
                                                                    iblock%variable(i+1,j+1,k+1)%res_3
                        !> restric rhs[4]
                        iblock_c%variable_mg(ii,jj,kk)%res4_mg =    iblock%variable(i  ,j  ,k  )%res_4  + &
                                                                    iblock%variable(i+1,j  ,k  )%res_4  + &
                                                                    iblock%variable(i  ,j+1,k  )%res_4  + &
                                                                    iblock%variable(i  ,j  ,k+1)%res_4  + &
                                                                    iblock%variable(i  ,j+1,k+1)%res_4  + &
                                                                    iblock%variable(i+1,j+1,k  )%res_4  + &
                                                                    iblock%variable(i+1,j  ,k+1)%res_4  + &
                                                                    iblock%variable(i+1,j+1,k+1)%res_4
                        !> restric rhs[5]
                        iblock_c%variable_mg(ii,jj,kk)%res5_mg =    iblock%variable(i  ,j  ,k  )%res_5  + &
                                                                    iblock%variable(i+1,j  ,k  )%res_5  + &
                                                                    iblock%variable(i  ,j+1,k  )%res_5  + &
                                                                    iblock%variable(i  ,j  ,k+1)%res_5  + &
                                                                    iblock%variable(i  ,j+1,k+1)%res_5  + &
                                                                    iblock%variable(i+1,j+1,k  )%res_5  + &
                                                                    iblock%variable(i+1,j  ,k+1)%res_5  + &
                                                                    iblock%variable(i+1,j+1,k+1)%res_5
                        !>
                        !>
                    end do
                end do
            end do
!            do k=2,kdim_c
!                do j=2,jdim_c
!                    do i=2,idim_c
!                        write(93001,'("flow=",3i4,5e24.16)') i-1,j-1,k-1,iblock_c%cell(i,j,k)%r,&
!                        iblock_c%cell(i,j,k)%u,iblock_c%cell(i,j,k)%v,iblock_c%cell(i,j,k)%w,&
!                        iblock_c%cell(i,j,k)%p
!                        write(93002,'("res=",3i4,6e24.16)') i-1,j-1,k-1,iblock_c%turbulent(i,j,k)%viscous,&
!                        iblock_c%variable_mg(i,j,k)%res1_mg,iblock_c%variable_mg(i,j,k)%res2_mg,&
!                        iblock_c%variable_mg(i,j,k)%res3_mg,iblock_c%variable_mg(i,j,k)%res4_mg,&
!                        iblock_c%variable_mg(i,j,k)%res5_mg
!                    end do
!                end do
!            end do
        else
            !> \brief
            !> two-dimension mesh operating
            !> i-direction dimension is 2
            !> restrict the r u v w to coarser mesh from finer mesh
            !>
            ii = 2
            i  = 2
            kk = 1
            do k=2,kdim,2
                kk = kk + 1
                jj = 1
                do j =2,jdim,2
                        jj = jj + 1
                        iblock_c%cell(ii,jj,kk)%r =(q_v(i,j  ,k  ,1) +&
                                                    q_v(i,j+1,k  ,1) +&
                                                    q_v(i,j  ,k+1,1) +&
                                                    q_v(i,j+1,k+1,1))*iblock_c%metric(ii,jj,kk)%jacobian
                        !> restric u
                        iblock_c%cell(ii,jj,kk)%u =(q_v(i,j  ,k  ,2) +&
                                                    q_v(i,j+1,k  ,2) +&
                                                    q_v(i,j  ,k+1,2) +&
                                                    q_v(i,j+1,k+1,2))*iblock_c%metric(ii,jj,kk)%jacobian !/ iblock_c%cell(ii,jj,kk)%r
                        !> restric v
                        iblock_c%cell(ii,jj,kk)%v =(q_v(i,j  ,k  ,3) +&
                                                    q_v(i,j+1,k  ,3) +&
                                                    q_v(i,j  ,k+1,3) +&
                                                    q_v(i,j+1,k+1,3))*iblock_c%metric(ii,jj,kk)%jacobian !/ iblock_c%cell(ii,jj,kk)%r
                        !> restric w
                        iblock_c%cell(ii,jj,kk)%w =(q_v(i,j  ,k  ,4) +&
                                                    q_v(i,j+1,k  ,4) +&
                                                    q_v(i,j  ,k+1,4) +&
                                                    q_v(i,j+1,k+1,4))*iblock_c%metric(ii,jj,kk)%jacobian !/ iblock_c%cell(ii,jj,kk)%r
                        !> restric the pressure
                        iblock_c%cell(ii,jj,kk)%p =(q_v(i,j  ,k  ,5) +&
                                                    q_v(i,j+1,k  ,5) +&
                                                    q_v(i,j  ,k+1,5) +&
                                                    q_v(i,j+1,k+1,5))*iblock_c%metric(ii,jj,kk)%jacobian
                        !> \brief
                        !>if viscous is not zero ,must restric from fine meshs to coarser meshs
                        !!
                        if(nvisc .gt. 1 )then
                            !> restrict the turbulent viscous from the finer mesh to coarse mesh
                            !>
                            iblock_c%turbulent(ii,jj,kk)%viscous = 0.25*( iblock%turbulent(i  ,j  ,k  )%viscous + &
                                                                          iblock%turbulent(i  ,j+1,k  )%viscous + &
                                                                          iblock%turbulent(i  ,j  ,k+1)%viscous + &
                                                                          iblock%turbulent(i  ,j+1,k+1)%viscous)
                        endif
                        !> the parameter for the SA model
                        if(nvisc .eq. 3)then
                            !> restrict the turbulent viscous from the finer mesh to coarse mesh
                            !>
                            iblock_c%turbulent(ii,jj,kk)%tur_save= 0.25*( iblock%turbulent(i  ,j  ,k  )%tur_save + &
                                                                          iblock%turbulent(i  ,j+1,k  )%tur_save + &
                                                                          iblock%turbulent(i  ,j  ,k+1)%tur_save + &
                                                                          iblock%turbulent(i  ,j+1,k+1)%tur_save)
                        endif
                end do
            end do
            !> \brief
            !>two-dimension mesh operating
            !!
            ii = 2
            i  = 2
            kk = 1
            do k=2,kdim,2
                kk = kk + 1
                jj = 1
                do j =2,jdim,2
                        jj = jj + 1
                        !> restric rhs[1]
                        iblock_c%variable_mg(ii,jj,kk)%res1_mg = iblock%variable(i,j  ,k  )%res_1 + &
                                                                 iblock%variable(i,j+1,k  )%res_1 + &
                                                                 iblock%variable(i,j  ,k+1)%res_1 + &
                                                                 iblock%variable(i,j+1,k+1)%res_1
                        !> restric rhs[2]
                        iblock_c%variable_mg(ii,jj,kk)%res2_mg = iblock%variable(i,j  ,k  )%res_2 + &
                                                                 iblock%variable(i,j+1,k  )%res_2 + &
                                                                 iblock%variable(i,j  ,k+1)%res_2 + &
                                                                 iblock%variable(i,j+1,k+1)%res_2
                        !> restric rhs[3]
                        iblock_c%variable_mg(ii,jj,kk)%res3_mg = iblock%variable(i,j  ,k  )%res_3 + &
                                                                 iblock%variable(i,j+1,k  )%res_3 + &
                                                                 iblock%variable(i,j  ,k+1)%res_3 + &
                                                                 iblock%variable(i,j+1,k+1)%res_3
                        !> restric rhs[4]
                        iblock_c%variable_mg(ii,jj,kk)%res4_mg = iblock%variable(i,j  ,k  )%res_4 + &
                                                                 iblock%variable(i,j+1,k  )%res_4 + &
                                                                 iblock%variable(i,j  ,k+1)%res_4 + &
                                                                 iblock%variable(i,j+1,k+1)%res_4
                        !> restric rhs[5]
                        iblock_c%variable_mg(ii,jj,kk)%res5_mg = iblock%variable(i,j  ,k  )%res_5 + &
                                                                 iblock%variable(i,j+1,k  )%res_5 + &
                                                                 iblock%variable(i,j  ,k+1)%res_5 + &
                                                                 iblock%variable(i,j+1,k+1)%res_5
                end do
            end do
        endif

        !> \brief
        !> fills the edges of the q array for safety  using
        !> multi-plane vectorization technique
        !!
!> the fill have some wrong!
!>
#if 0
        if(idim .gt. 1)then

            do  i=2,idim_c
                !> j-direction and i,k-face at jdim_c
                !!
                do k =2,kdim_c
                    iblock_c%cell(i,jdim_c+1,k)%r = iblock_c%cell(i,jdim_c,k)%r
                    iblock_c%cell(i,jdim_c+1,k)%u = iblock_c%cell(i,jdim_c,k)%u
                    iblock_c%cell(i,jdim_c+1,k)%v = iblock_c%cell(i,jdim_c,k)%v
                    iblock_c%cell(i,jdim_c+1,k)%w = iblock_c%cell(i,jdim_c,k)%w
                    iblock_c%cell(i,jdim_c+1,k)%p = iblock_c%cell(i,jdim_c,k)%p
                    !>
                    !>
                    iblock_c%cell(i,1,k)%r = iblock_c%cell(i,2,k)%r
                    iblock_c%cell(i,1,k)%u = iblock_c%cell(i,2,k)%u
                    iblock_c%cell(i,1,k)%v = iblock_c%cell(i,2,k)%v
                    iblock_c%cell(i,1,k)%w = iblock_c%cell(i,2,k)%w
                    iblock_c%cell(i,1,k)%p = iblock_c%cell(i,2,k)%p

                end do
                !> k-direction and i,j-face at kdim_c
                !!
                do j =2,jdim_c
                    iblock_c%cell(i,j,kdim_c+1)%r = iblock_c%cell(i,j,kdim_c)%r
                    iblock_c%cell(i,j,kdim_c+1)%u = iblock_c%cell(i,j,kdim_c)%v
                    iblock_c%cell(i,j,kdim_c+1)%v = iblock_c%cell(i,j,kdim_c)%u
                    iblock_c%cell(i,j,kdim_c+1)%w = iblock_c%cell(i,j,kdim_c)%w
                    iblock_c%cell(i,j,kdim_c+1)%p = iblock_c%cell(i,j,kdim_c)%p
                    !>
                    !>
                    iblock_c%cell(i,j,1)%r = iblock_c%cell(i,j,2)%r
                    iblock_c%cell(i,j,1)%u = iblock_c%cell(i,j,2)%u
                    iblock_c%cell(i,j,1)%v = iblock_c%cell(i,j,2)%v
                    iblock_c%cell(i,j,1)%w = iblock_c%cell(i,j,2)%w
                    iblock_c%cell(i,j,1)%p = iblock_c%cell(i,j,2)%p

               end do
            end do
            !> i-direction and j,k-face at idim_c
            ! !
            do j=2,jdim_c
                do k =2,kdim_c
                    iblock_c%cell(idim_c+1,j,k)%r = iblock_c%cell(idim_c,j,k)%r
                    iblock_c%cell(idim_c+1,j,k)%u = iblock_c%cell(idim_c,j,k)%u
                    iblock_c%cell(idim_c+1,j,k)%v = iblock_c%cell(idim_c,j,k)%v
                    iblock_c%cell(idim_c+1,j,k)%w = iblock_c%cell(idim_c,j,k)%w
                    iblock_c%cell(idim_c+1,j,k)%p = iblock_c%cell(idim_c,j,k)%p

                    !>
                    !>
                    iblock_c%cell(1,j,k)%r = iblock_c%cell(2,j,k)%r
                    iblock_c%cell(1,j,k)%u = iblock_c%cell(2,j,k)%u
                    iblock_c%cell(1,j,k)%v = iblock_c%cell(2,j,k)%v
                    iblock_c%cell(1,j,k)%w = iblock_c%cell(2,j,k)%w
                    iblock_c%cell(1,j,k)%p = iblock_c%cell(2,j,k)%p

                end do
            end do
            !>
        else
            !> j-direction and i,k-face at jdim_c
            !!
            i = 1
            do k =2,kdim_c
                iblock_c%cell(i,jdim_c+1,k)%r = iblock_c%cell(i,jdim_c,k)%r
                iblock_c%cell(i,jdim_c+1,k)%u = iblock_c%cell(i,jdim_c,k)%u
                iblock_c%cell(i,jdim_c+1,k)%v = iblock_c%cell(i,jdim_c,k)%v
                iblock_c%cell(i,jdim_c+1,k)%w = iblock_c%cell(i,jdim_c,k)%w
                iblock_c%cell(i,jdim_c+1,k)%p = iblock_c%cell(i,jdim_c,k)%p

                !>
                !>
                iblock_c%cell(i,1,k)%r = iblock_c%cell(i,2,k)%r
                iblock_c%cell(i,1,k)%u = iblock_c%cell(i,2,k)%u
                iblock_c%cell(i,1,k)%v = iblock_c%cell(i,2,k)%v
                iblock_c%cell(i,1,k)%w = iblock_c%cell(i,2,k)%w
                iblock_c%cell(i,1,k)%p = iblock_c%cell(i,2,k)%p

            end do
            !> k-direction and i,j-face at kdim_c
            !!
            do j =2,jdim_c

                iblock_c%cell(i,j,kdim_c+1)%r = iblock_c%cell(i,j,kdim_c)%r
                iblock_c%cell(i,j,kdim_c+1)%u = iblock_c%cell(i,j,kdim_c)%u
                iblock_c%cell(i,j,kdim_c+1)%v = iblock_c%cell(i,j,kdim_c)%v
                iblock_c%cell(i,j,kdim_c+1)%w = iblock_c%cell(i,j,kdim_c)%w
                iblock_c%cell(i,j,kdim_c+1)%p = iblock_c%cell(i,j,kdim_c)%p

                !>
                !>
                iblock_c%cell(i,j,1)%r = iblock_c%cell(i,j,2)%r
                iblock_c%cell(i,j,1)%u = iblock_c%cell(i,j,2)%u
                iblock_c%cell(i,j,1)%v = iblock_c%cell(i,j,2)%v
                iblock_c%cell(i,j,1)%w = iblock_c%cell(i,j,2)%w
                iblock_c%cell(i,j,1)%p = iblock_c%cell(i,j,2)%p

            end do

        end if
#endif
        !> \brief
        !> the q_v array must delete before return the subroutine
        !!
        deallocate(q_v)
        deallocate(res_v)
        !> end the subroutine
        !>
        return

    end subroutine restrict
