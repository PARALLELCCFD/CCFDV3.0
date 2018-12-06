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
    !> purpose:compute the tau term for correct and store the value i(n,n-1)q(n),
    !> the i(n,n-1)q(n) will used for correct the finer mesh at prolonged
    !> computing the tau term
    !> tau = l(i(n,n-1)q(n))-i(n,n-1)(l(q(n))) + i(n,n-1)tau(n)
    !> use in detemining correction to finer mesh in the multigrid iteration scheme
    !> by liuxz 2016/12/25

    subroutine tau_term(nbl,level)
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
		implicit none
		!>
		!>
        integer::idim,jdim,kdim
        integer:: i,j,k,n,level,nbl
        real(kind= dprec) :: qq
        type(blocks_type),pointer  :: iblock
        type(overlap_type),pointer :: mesh
        !> load the blocks is demension indiex correction
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        idim = iblock%idim
        jdim = iblock%jdim
        kdim = iblock%kdim
        !>
        !> store the primitive variables for later use in the prolonged
        !> determining delta q (when mode = 1,kode = 2)
        !> residual correction (when mode = 1,kode = 1)
        !> kode=1   mode=0  standard update
        !>                  rhs=roe(q)
        !>                  lu-sgs(delaq)=rhs(roe)
        !>
        !> kode=2   mode=1  first time on lower level multi-grid
        !>                  computing the tau term
        !>                  tau(n)=rhs(i(n,n-1)q(n))-i(n,n-1)(rhs(q(n)))
        !>                  detail of the term can see user's manual
        !>                  it used after comoputing the rhs(i(n,n-1)q(n))
        !> kode=1   mode=1   subsequnet time on lower level multi-grid
        !>                   add the tau(n) of n level
        !>                   rhs(q(n))+tau(n) is will used at n+1 level of mesh
        !>                   lu-sgs(delaq)=rhs(roe)

        if(kode .ge.2 .and. mgflags .ne.0)then
             !
            if(level .ge. global_level) goto 100! .and. mgflags .lt. 2)  goto 100
            !> store the variables
            do k=2,kdim
                do j=2,jdim
                    do i=2,idim
                        !> store the rho
                        iblock%variable_mg(i,j,k)%r_mg   = iblock%cell(i,j,k)%r
                        !> store the u
                        iblock%variable_mg(i,j,k)%u_mg   = iblock%cell(i,j,k)%u
                        !> store the v
                        iblock%variable_mg(i,j,k)%v_mg   = iblock%cell(i,j,k)%v
                        !> store the w
                        iblock%variable_mg(i,j,k)%w_mg   = iblock%cell(i,j,k)%w
                        !> store the p
                        iblock%variable_mg(i,j,k)%p_mg   = iblock%cell(i,j,k)%p
                    end do
                end do
            end do
            !>
            !>
        end if
100     continue
        !>  kode =2 and mode = 1 then residual restrict from finer meshs is rn_res etc
        !>  rn_res = rn_res - rn (rn is computing by the flux scheme)
        !>  kode =1 and mode = 1 then  rn_res = rn_res + rn  correction the residual
        !>
        if(mode .ne. 0)then
            if(level .lt. global_level)then !no-embedded grid
                !> rn_res is restrict from finer mesh ,so at the globla_level
                !> rn_res will not used
                !> rn_res = rn_res - rn (rn is computing by the flux scheme)
                !> residual correction from global grids
                !> res-i(n,n-1)(l(q(n)))
                if(kode .ge. 2)then
                    do k=2,kdim
                        do j=2,jdim
                            do i=2,idim
                                !> corrction the rho by restrict values
                                iblock%variable_mg(i,j,k)%res1_mg = iblock%variable_mg(i,j,k)%res1_mg - iblock%variable(i,j,k)%res_1
                                !> corrction the rho by restrict values
                                iblock%variable_mg(i,j,k)%res2_mg = iblock%variable_mg(i,j,k)%res2_mg - iblock%variable(i,j,k)%res_2
                                !> corrction the rho by restrict value
                                iblock%variable_mg(i,j,k)%res3_mg = iblock%variable_mg(i,j,k)%res3_mg - iblock%variable(i,j,k)%res_3
                                !> corrction the rho by restrict values
                                iblock%variable_mg(i,j,k)%res4_mg = iblock%variable_mg(i,j,k)%res4_mg - iblock%variable(i,j,k)%res_4
                                !> corrction the rho by restrict values
                                iblock%variable_mg(i,j,k)%res5_mg = iblock%variable_mg(i,j,k)%res5_mg - iblock%variable(i,j,k)%res_5

                            end do
                        end do
                    end do
                endif
                !elseif( kode .eq. 1)then
                !> kode = 1 and mode =1 is rn_res = rn_res + rn  correction the residual
                !>i(n,n-1)(l(q(n))) + res
                do k=2,kdim
                    do j=2,jdim
                        do i=2,idim
!                            write(200+nbl,'("tau=",4i4,2x,5e24.16)')nbl,i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,&
!                            iblock%variable(i,j,k)%res_3,&
!                            iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
                            !> rn_res = rn_res + rn
                            iblock%variable(i,j,k)%res_1   = iblock%variable_mg(i,j,k)%res1_mg   + iblock%variable(i,j,k)%res_1
                            !> run_res = run_res + run
                            iblock%variable(i,j,k)%res_2   = iblock%variable_mg(i,j,k)%res2_mg   + iblock%variable(i,j,k)%res_2
                            !> rvn_res = rvn_res + rvn
                            iblock%variable(i,j,k)%res_3   = iblock%variable_mg(i,j,k)%res3_mg   + iblock%variable(i,j,k)%res_3
                            !> rwn_res = rwn_res + rwn
                            iblock%variable(i,j,k)%res_4   = iblock%variable_mg(i,j,k)%res4_mg   + iblock%variable(i,j,k)%res_4
                            !> ren_res = ren_res + ren
                            iblock%variable(i,j,k)%res_5   = iblock%variable_mg(i,j,k)%res5_mg   + iblock%variable(i,j,k)%res_5
!                            write(210+nbl,'("tau=",4i4,2x,5e24.16)')nbl,i,j,k,iblock%variable_mg(i,j,k)%res1_mg,iblock%variable_mg(i,j,k)%res2_mg,&
!                            iblock%variable_mg(i,j,k)%res3_mg,&
!                            iblock%variable_mg(i,j,k)%res4_mg,iblock%variable_mg(i,j,k)%res5_mg
!                            write(220+nbl,'("tau=",4i4,2x,5e24.16)')nbl,i,j,k,iblock%variable(i,j,k)%res_1,iblock%variable(i,j,k)%res_2,&
!                            iblock%variable(i,j,k)%res_3,&
!                            iblock%variable(i,j,k)%res_4,iblock%variable(i,j,k)%res_5
                        end do
                    end do
                end do
                !end if
                !> end the  correction
            !> residual correction from the embedded grids
            !> if is need it can expend by user
            !>
            else
                !> null
                !> null
                !> null
            end if
        end if
        !> end the subroutine
        return
    end subroutine tau_term

