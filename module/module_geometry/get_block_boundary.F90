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
    subroutine get_cgns_blocks_boundary
        use global_parameter
        use blocks_module
        use mesh_overlap_module
        use bc_module
        use nodes_paras
        use nodes_var
        use nodes_var_bc
		!>
		!>
!		implicit none
        include "cgnslib_f.h"
!
#if defined PMPI
        include "mpif.h"
!        character(len=1),allocatable :: packbuff(:)
!        integer :: position1,buff_size
#endif
        integer stats,num1to1,n_bc,n_121,ier,nbl,nb,nbocos,ibocotype,iii,jjj
        integer isize(3,3),ibase,iccg,ierr,nccg
        integer,dimension(:,:),allocatable :: irange
        integer,dimension(:,:),allocatable :: idonor_range
        integer,dimension(:,:),allocatable :: itransform
!        integer,dimension(:),allocatable   :: bc_send
        type(overlap_type), pointer        :: mesh
!********************* imhoy  for  pointwise type time:2017-1-12 14:14:48******************
        type bc_family_type
        integer::bctype,userdefinebc
        character(len=50)::familyname,fambcname !familyname是family_t的名字，fambcname是familybc_t的名字
        end type bc_family_type
        type(bc_family_type),pointer,dimension(:)::bc_fam
        !>
        !>
        !>
        character(len=50)::familyname,fambcname
!********************* imhoy  for  pointwise type time:2017-1-12 14:14:48******************
        character*32,dimension(:),allocatable::connectname,znname,donorname,zonename
        character*32 boconame
        integer normalindex(3),ipnts(6),ipoint,ii,bc_type
        real(kind = dprec) data_double(6)
        allocate(zonename(nobl))
        mesh => grids(imesh)

#if defined PMPI
        if(myid == myhost)then
#endif
    !  Read CGNS file
            ibase  = 1
            ier    = 1
            !IER    = 0
            !write(*,'(" ibase= ",i3," ier=",i3)') ibase,ier
            open (11,file =grid_files,status='old')
            call cg_open_f(grid_files,MODE_READ,iccg,ier)
            if (ier .ne. 0) call cg_error_exit_f
            !***************************************
            do n=1,nobl
            !   Read general zone information
                call cg_zone_read_f(iccg, ibase, n, zonename(n), isize, ier)
                if (ier .ne. 0) call cg_error_exit_f
            enddo
          !******************************************
    !   Now get 1-to-1 connectivity info
            call cg_n1to1_global_f(iccg,ibase,num1to1,ier)
          !
            ALLOCATE(IRANGE(6,NUM1TO1),STAT=STATS)
            allocate( idonor_range(6,num1to1), stat=stats )
            allocate( itransform(3,num1to1), stat=stats )
            allocate(connectname(num1to1),STAT=STATS)
            allocate(znname(num1to1),STAT=STATS)
            allocate(donorname(num1to1),STAT=STATS)
            !
            if (ier .ne. 0) call cg_error_exit_f
            call cg_1to1_read_global_f(iccg,ibase,connectname,znname,&
            donorname,irange,idonor_range,itransform,ier)

            if (ier .ne. 0) call cg_error_exit_f

    !********************* imhoy  for  POINTWISE type Time:2017-1-12 14:14:48******************
            call cg_nfamilies_f(iccg, ibase, nfamilies, ier)
            write(6,*)"POINTWISE type BC FAMILY NUM=:",nfamilies
            allocate( bc_fam(nfamilies), stat=stats )
            do i=1,nfamilies
                call cg_family_read_f(iccg,ibase,i,FamilyName,nFamBC,nGeo,ier)
                if (nFamBC==1)then
                    bc_fam(i)%familyName = FamilyName
                    call cg_fambc_read_f(iccg,ibase,i,1,FamBCName,ifambctype,ier)
                    bc_fam(i)%FamBCName = FamBCName
                    bc_fam(i)%BCType = ifambctype
    !         call upcase(FamilyName)
                    write(*,*)"BCfamily: ",trim(familyName)," bctype:",ifambctype
                endif
            enddo
    !********************* imhoy  for  POINTWISE type Time:2017-1-12 14:14:48******************
            !>
            !>
            n_bc = 0
            do n=1,nobl
                call cg_nbocos_f(iccg, ibase, n, nbocos, ier)
                n_bc = n_bc + nbocos
!                write(61,'("the nbl=",i6," have the physical bc =",i6)') n,nbocos
            end do
            write(61,'("the total physical =",i6)') n_bc
            n_bc = n_bc + num1to1*2
#if defined DEBUG
            write(61,*)'**********************************************************'
            write(61,'('' 1-1blocks num= '',i10,", bc_num= ",i10)') num1to1*2,n_bc
            write(61,*)'** boundary condition and coordinate informations*********'
#endif

            !> set the multigrid method's bonudary conditions
            !> set the coarser meshs bc
            !>
            if(mgflags .ne. 0)then
                 allocate(mesh%bcs(n_bc*global_level))
            else
                 allocate(mesh%bcs(n_bc))
            end if

            !    ################################################################

            !>   see if there is any bc info
            n_bc  = 1
            !>
            do n=1,nobl
                !>
                call cg_nbocos_f(iccg, ibase, n, nbocos, ier)
                if (ier .ne. 0) call cg_error_exit_f
#if defined DEBUG
                write(61,*) "*********************************************************"
                write(61,'(" The original block=",i6," physics boundary conditions informations")') n
#endif
                !
                if (nbocos .gt. 0) then
                    do m=1,nbocos
                        call cg_boco_info_f(iccg, ibase, n, m, boconame, ibocotype,&
                        iptset_type, npnts, normalindex, normalflag,&
                        idatatype, ndataset, ier)
                        if (ier .ne. 0) call cg_error_exit_f
                        if (npnts .ne. 2) then
                            write(*,'('' cnoblannot use... npnts must be 2'')')
                            cycle!goto 333
                        end if
                        if (iptset_type .ne. PointRange) then
                            write(*,'('' cannot use... iptset_type must be'',&
                                '' PointRange'')')
                            cycle !goto 333
                        end if
                        !> \brief
                        !!
                        !
!********************* imhoy  for  pointwise type time:2017-1-12 14:14:48******************
                        if (ibocotype==FamilySpecified)then
                            call cg_goto_f(iccg,ibase,ier,"Zone_t",n,"ZoneBC_t",1,&
                                   "BC_t",m,"end")  !n:blocknumber ,m:bcnumber
                            call cg_famname_read_f(FamilyName,ier)
                            do i=1,nfamilies
                                if (trim(bc_fam(i)%familyName)==trim(FamilyName))then
                                    ibocotype= bc_fam(i)%BCType
                                    exit
                                endif
                            enddo
                        endif
!********************* imhoy  for  pointwise type time:2017-1-12 14:14:48******************
						!>
						!>BCFarfield=7
                        if( ibocotype .eq. BCFarfield )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'farfield'
                            !>
                            !>
						!>BCSymmetryPlane=16
                        else if(ibocotype .eq. BCSymmetryPlane )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'symmetryplane'
                            !>
                            !>
						!>BCWallInviscid=21
                        else if(ibocotype .eq. BCWallInviscid )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'wallinviscid'
                            !>
                            !>
                        !>BCWall=20
						!>BCWallViscous=22
						!>BCWallViscousHeatFlux=23
						!>BCWallViscousIsothermal=24
                        else if(ibocotype .eq. BCWall &
                               .or. ibocotype .eq. BCWallViscous&
                               .or. ibocotype .eq. BCWallViscousHeatFlux &
                               .or. ibocotype .eq. BCWallViscousIsothermal )then
                            !>
                            !>
                            if(nvisc .gt. 0 )then
                                mesh%bcs(n_bc)%bc_type =  'wallviscid'
                            else
                                !> set the wall viscid to the inviscid wall bc
                                mesh%bcs(n_bc)%bc_type =  'wallinviscid'
                            end if
                            !>
                            !>
						!>BCInflow=9
                        else if(ibocotype .eq. BCInflow )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'inflow'
                            !>
                            !>
						!>BCExtrapolate=6
                        else if(ibocotype .eq. BCExtrapolate )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'extrapolat'
                            !>
                            !>
						!>BCInflowSupersonic=11
                        else if(ibocotype .eq. BCInflowSupersonic )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'inflowsupersonic'
                            !>
                            !>
						!>BCTunnelInflow=18
                        else if(ibocotype .eq. BCTunnelInflow )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'tunnelinflow'
                            !>
                            !>
						!>BCSymmetryPolar=17
                        else if(ibocotype .eq. BCSymmetryPolar )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'symmetrypolar'
                            !>
                            !>
						!>BCDegenerateLine=3
                        else if(ibocotype .eq. BCDegenerateLine )then
                            !>
                            !>
                            mesh%bcs(n_bc)%bc_type =  'degenerateline'
                            !>
                            !>
                        else
                            write(*,*) '*********************************************'
                            write(*,*) 'ccfdv3.0 error: have some undefined bc_type  '
                            write(*,*) 'the bc_type is ',bctypename(ibocotype)
                            write(*,*) 'please check the bc_type before computing!'
                            write(*,*) '*********************************************'
                            cycle
                        endif
                        !>
                        call cg_boco_read_f(iccg, ibase, n, m, ipnts,data_double, ier)
                        if (ier .ne. 0) call cg_error_exit_f
                        if ((ipnts(1).eq.ipnts(4)).and.(ipnts(2).eq.ipnts(5)).or.&
                            (ipnts(1).eq.ipnts(4)).and.(ipnts(3).eq.ipnts(6)).or.&
                            (ipnts(2).eq.ipnts(5)).and.(ipnts(3).eq.ipnts(6))) then
                            cycle
                        end if
!
                        mesh%bcs(n_bc)%istart = ipnts(1)
                        mesh%bcs(n_bc)%jstart = ipnts(2)
                        mesh%bcs(n_bc)%kstart = ipnts(3)
                        mesh%bcs(n_bc)%iend   = ipnts(4)
                        mesh%bcs(n_bc)%jend   = ipnts(5)
                        mesh%bcs(n_bc)%kend   = ipnts(6)
                        mesh%bcs(n_bc)%iblock_index   = 0


                        !> multigrid methods and the cgns grid original num of blocks
                        !> so must computing the num of mg methods blocks
                        !>
                        if(mgflags .ne. 0 )then
                            mesh%bcs(n_bc)%block_index   = (n-1)*(global_level-1) + n
                        else
                            mesh%bcs(n_bc)%block_index   = n
                        end if

#if defined DEBUG
                        open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                        write(61,'(" block bc_num istart   iend jstart   jend kstart   kend  bc_type  ")')
                        write(61,'(i6,7i7,2x,a25)') &
                                    mesh%bcs(n_bc)%block_index,&
                                    n_bc,&
                                    mesh%bcs(n_bc)%istart,&
                                    mesh%bcs(n_bc)%iend,&
                                    mesh%bcs(n_bc)%jstart,&
                                    mesh%bcs(n_bc)%jend,&
                                    mesh%bcs(n_bc)%kstart,&
                                    mesh%bcs(n_bc)%kend,&
                                    mesh%bcs(n_bc)%bc_type
#endif
                        !> multigrid methods bc for coarser meshs
                        !> ijk-dir dimensions and bc type
                        !> bc at the nblocks
                        !>
                        !>
                        if(mgflags .ne. 0)then
                            do nccg=1,global_level-1
                                !>
                                !>
                                n_bc               = n_bc + 1
                                mesh%bcs(n_bc)%block_index = mesh%bcs(n_bc-1)%block_index + 1
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%istart .eq. 1)then
                                    mesh%bcs(n_bc)%istart   = 1
                                else
                                    mesh%bcs(n_bc)%istart   = mesh%bcs(n_bc-1)%istart/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jstart .eq. 1)then
                                    mesh%bcs(n_bc)%jstart   = 1
                                else
                                    mesh%bcs(n_bc)%jstart   = mesh%bcs(n_bc-1)%jstart/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kstart .eq. 1)then
                                    mesh%bcs(n_bc)%kstart   = 1
                                else
                                    mesh%bcs(n_bc)%kstart   = mesh%bcs(n_bc-1)%kstart/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%iend .eq. 1)then
                                    mesh%bcs(n_bc)%iend   = 1
                                else
                                    mesh%bcs(n_bc)%iend   = mesh%bcs(n_bc-1)%iend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jend .eq. 1)then
                                    mesh%bcs(n_bc)%jend   = 1
                                else
                                    mesh%bcs(n_bc)%jend   = mesh%bcs(n_bc-1)%jend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kend .eq. 1)then
                                    mesh%bcs(n_bc)%kend   = 1
                                else
                                    mesh%bcs(n_bc)%kend   = mesh%bcs(n_bc-1)%kend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                mesh%bcs(n_bc)%bc_type = mesh%bcs(n_bc-1)%bc_type
                                mesh%bcs(n_bc)%iblock_index   = 0
!                                mesh%bcs(n_bc)%surface=  mesh%bcs(n_bc-1)%surface
                                !>
                                !>
                                !>
#if defined DEBUG
                                open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                                write(61,'(" block_id bc_num istart   iend jstart   jend kstart   kend  bc_type  ")')
                                write(61,'(i9,7i7,2x,a25)')&
                                              mesh%bcs(n_bc)%block_index,&
                                              n_bc,&
                                              mesh%bcs(n_bc)%istart,&
                                              mesh%bcs(n_bc)%iend,&
                                              mesh%bcs(n_bc)%jstart,&
                                              mesh%bcs(n_bc)%jend,&
                                              mesh%bcs(n_bc)%kstart,&
                                              mesh%bcs(n_bc)%kend,&
                                              mesh%bcs(n_bc)%bc_type
#endif
                            end do
                            n_bc = n_bc + 1
                        !> no multigrid methods
                        !> one level meshs
                        else
                            n_bc = n_bc + 1
                        end if
                    enddo
                end if

                if(mgflags .ne. 0)then
                    nb = (n-1)*(global_level-1) + n
                else
                    nb = n
                end if
#if defined DEBUG
                open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                write(61,'(" The original block=",i6,"    1-1 boundary conditions informations")') n
                write(61,*) "*****************************************************************"
#endif
!                write(61,'("original block=",i4," have ",i4," num of bcs")') n,num1to1
                do n_121 =1,num1to1
                    if(znname(n_121) .eq. zonename(n))then
                        !>
                        !>
                        !>
                        mesh%bcs(n_bc)%bc_type = 'cut'
                        !>
                        !>
                        mesh%bcs(n_bc)%block_index   =  nb
                        !>
                        !>
!                        nn = donorname(n_121)
!                        if(mgflags .ne. 0)then
!                            mesh%bcs(n_bc)%iblock_index = (nn-1)*(global_level-1) + nn
!                        else
!                            mesh%bcs(n_bc)%iblock_index = nn
!                        end if
                        do nn=1,nobl
                            if(donorname(n_121).eq. zonename(nn) )then
                                if(mgflags .ne. 0)then
                                    mesh%bcs(n_bc)%iblock_index = (nn-1)*(global_level-1) + nn
                                else
                                    mesh%bcs(n_bc)%iblock_index = nn
                                end if
                            end if
                        end do

!                        nn = mesh%bcs(n_bc)%iblock_index
                        !> \brief
                        !!
                        !> the start and the end dimension of the cut bc
                        !>
                        mesh%bcs(n_bc)%istart  =  irange(1,n_121)
                        mesh%bcs(n_bc)%jstart  =  irange(2,n_121)
                        mesh%bcs(n_bc)%kstart  =  irange(3,n_121)
                        mesh%bcs(n_bc)%iend    =  irange(4,n_121)
                        mesh%bcs(n_bc)%jend    =  irange(5,n_121)
                        mesh%bcs(n_bc)%kend    =  irange(6,n_121)
!                        write(161,'("n_bc=",8i8)') n_bc,n_121,mesh%bcs(n_bc)%istart,mesh%bcs(n_bc)%jstart,mesh%bcs(n_bc)%kstart,mesh%bcs(n_bc)%iend,&
!                        mesh%bcs(n_bc)%jend,mesh%bcs(n_bc)%kend
                        !>
                        !>
                        !> the corresponding dimension
                        !>
                        mesh%bcs(n_bc)%istart_cut =  idonor_range(1,n_121)
                        mesh%bcs(n_bc)%jstart_cut =  idonor_range(2,n_121)
                        mesh%bcs(n_bc)%kstart_cut =  idonor_range(3,n_121)
                        mesh%bcs(n_bc)%iend_cut   =  idonor_range(4,n_121)
                        mesh%bcs(n_bc)%jend_cut   =  idonor_range(5,n_121)
                        mesh%bcs(n_bc)%kend_cut   =  idonor_range(6,n_121)
                        !
                        !>
                        !>
                        !> the order of the bcs
                        mesh%bcs(n_bc)%irank  =  abs(itransform(1,n_121))
                        mesh%bcs(n_bc)%jrank  =  abs(itransform(2,n_121))
                        mesh%bcs(n_bc)%krank  =  abs(itransform(3,n_121))
!                        write(161,'("n_bc=",10i4)') n_bc, mesh%bcs(n_bc)%istart, mesh%bcs(n_bc)%iend,&
!                         mesh%bcs(n_bc)%jstart, mesh%bcs(n_bc)%jend, mesh%bcs(n_bc)%kstart, mesh%bcs(n_bc)%kend,&
!                         mesh%bcs(n_bc)%irank,mesh%bcs(n_bc)%jrank,mesh%bcs(n_bc)%krank

#if defined DEBUG
                        open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                        write(61,'(" block_id  bc_num bc_type   is   ie   js   je   ks   ke cut_block   is   ie   js   je   ks   ke ")')
                        write(61,'(i9,i7,1x,a6,6i5,i11,6i5)')   &
                                        mesh%bcs(n_bc)%block_index,&
                                        n_bc,&
                                        mesh%bcs(n_bc)%bc_type,&
                                        mesh%bcs(n_bc)%istart,&
                                        mesh%bcs(n_bc)%iend,&
                                        mesh%bcs(n_bc)%jstart,&
                                        mesh%bcs(n_bc)%jend,&
                                        mesh%bcs(n_bc)%kstart,&
                                        mesh%bcs(n_bc)%kend,&
                                        mesh%bcs(n_bc)%iblock_index,&
                                        mesh%bcs(n_bc)%istart_cut,&
                                        mesh%bcs(n_bc)%iend_cut,&
                                        mesh%bcs(n_bc)%jstart_cut,&
                                        mesh%bcs(n_bc)%jend_cut,&
                                        mesh%bcs(n_bc)%kstart_cut,&
                                        mesh%bcs(n_bc)%kend_cut!,n2p(mesh%bcs(n_bc)%block_index),n2p(mesh%bcs(n_bc)%iblock_index)

#endif
                        !> set the 1-1blocks for coarser meshs in multigrid methods
                        !> set the ijk-dir dimensions
                        !> and the 1-1 blocks of adjacent surfaces
                        !> the ijk-dir dimensions of adjacent surfaces
                        if(mgflags .ne. 0)then
                            do nccg=1,global_level-1
                                n_bc = n_bc + 1
                                !>
                                !>
                                !>
                                mesh%bcs(n_bc)%bc_type        = 'cut'
                                !>
                                !>
                                mesh%bcs(n_bc)%block_index    = mesh%bcs(n_bc-1)%block_index + 1
                                !>
                                !>
                                mesh%bcs(n_bc)%iblock_index   = mesh%bcs(n_bc-1)%iblock_index+ 1
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%istart .eq. 1)then
                                    mesh%bcs(n_bc)%istart   = 1
                                else
                                    mesh%bcs(n_bc)%istart   = mesh%bcs(n_bc-1)%istart/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jstart .eq. 1)then
                                    mesh%bcs(n_bc)%jstart   = 1
                                else
                                    mesh%bcs(n_bc)%jstart   = mesh%bcs(n_bc-1)%jstart/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kstart .eq. 1)then
                                    mesh%bcs(n_bc)%kstart   = 1
                                else
                                    mesh%bcs(n_bc)%kstart   = mesh%bcs(n_bc-1)%kstart/2 + 1
                                end if
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%iend .eq. 1)then
                                    mesh%bcs(n_bc)%iend   = 1
                                else
                                    mesh%bcs(n_bc)%iend   = mesh%bcs(n_bc-1)%iend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jend .eq. 1)then
                                    mesh%bcs(n_bc)%jend   = 1
                                else
                                   mesh%bcs(n_bc)%jend   = mesh%bcs(n_bc-1)%jend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kend .eq. 1)then
                                    mesh%bcs(n_bc)%kend   = 1
                                else
                                    mesh%bcs(n_bc)%kend   = mesh%bcs(n_bc-1)%kend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%istart_cut .eq. 1)then
                                    mesh%bcs(n_bc)%istart_cut   = 1
                                else
                                    mesh%bcs(n_bc)%istart_cut   = mesh%bcs(n_bc-1)%istart_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jstart_cut .eq. 1)then
                                    mesh%bcs(n_bc)%jstart_cut   = 1
                                else
                                    mesh%bcs(n_bc)%jstart_cut   = mesh%bcs(n_bc-1)%jstart_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kstart_cut .eq. 1)then
                                    mesh%bcs(n_bc)%kstart_cut   = 1
                                else
                                    mesh%bcs(n_bc)%kstart_cut   = mesh%bcs(n_bc-1)%kstart_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%iend_cut .eq. 1)then
                                    mesh%bcs(n_bc)%iend_cut   = 1
                                else
                                    mesh%bcs(n_bc)%iend_cut   = mesh%bcs(n_bc-1)%iend_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jend_cut .eq. 1)then
                                    mesh%bcs(n_bc)%jend_cut   = 1
                                else
                                    mesh%bcs(n_bc)%jend_cut   = mesh%bcs(n_bc-1)%jend_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kend_cut .eq. 1)then
                                    mesh%bcs(n_bc)%kend_cut   = 1
                                else
                                    mesh%bcs(n_bc)%kend_cut   = mesh%bcs(n_bc-1)%kend_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                mesh%bcs(n_bc)%irank = mesh%bcs(n_bc-1)%irank
                                mesh%bcs(n_bc)%jrank = mesh%bcs(n_bc-1)%jrank
                                mesh%bcs(n_bc)%krank = mesh%bcs(n_bc-1)%krank
                                !>
                                !>
                                !>
#if defined DEBUG
!                                open(unit=61,file='ccfd.out',form='formatted',status='unknown')
!                                write(61,'(" block bc_num bc_type   is   ie   js   je   ks   ke iblock   is   ie   js   je   ks   ke ")')
!                                write(61,'(i6,i7,1x,a6,6i5,i7,6i5)')   &
!                                                mesh%bcs(n_bc)%block_index,&
!                                                n_bc,&
!                                                mesh%bcs(n_bc)%bc_type,&
!                                                mesh%bcs(n_bc)%istart,&
!                                                mesh%bcs(n_bc)%iend,&
!                                                mesh%bcs(n_bc)%jstart,&
!                                                mesh%bcs(n_bc)%jend,&
!                                                mesh%bcs(n_bc)%kstart,&
!                                                mesh%bcs(n_bc)%kend,&
!                                                mesh%bcs(n_bc)%iblock_index,&
!                                                mesh%bcs(n_bc)%istart_cut,&
!                                                mesh%bcs(n_bc)%iend_cut,&
!                                                mesh%bcs(n_bc)%jstart_cut,&
!                                                mesh%bcs(n_bc)%jend_cut,&
!                                                mesh%bcs(n_bc)%kstart_cut,&
!                                                mesh%bcs(n_bc)%kend_cut
                                open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                                write(61,'(" block_id bc_num bc_type   is   ie   js   je   ks   ke cut_block   is   ie   js   je   ks   ke ")')
                                write(61,'(i9,i7,1x,a6,6i5,i11,6i5)')   &
                                                mesh%bcs(n_bc)%block_index,&
                                                n_bc,&
                                                mesh%bcs(n_bc)%bc_type,&
                                                mesh%bcs(n_bc)%istart,&
                                                mesh%bcs(n_bc)%iend,&
                                                mesh%bcs(n_bc)%jstart,&
                                                mesh%bcs(n_bc)%jend,&
                                                mesh%bcs(n_bc)%kstart,&
                                                mesh%bcs(n_bc)%kend,&
                                                mesh%bcs(n_bc)%iblock_index,&
                                                mesh%bcs(n_bc)%istart_cut,&
                                                mesh%bcs(n_bc)%iend_cut,&
                                                mesh%bcs(n_bc)%jstart_cut,&
                                                mesh%bcs(n_bc)%jend_cut,&
                                                mesh%bcs(n_bc)%kstart_cut,&
                                                mesh%bcs(n_bc)%kend_cut!,n2p(mesh%bcs(n_bc)%block_index),n2p(mesh%bcs(n_bc)%iblock_index)
#endif
                            end do
                            n_bc = n_bc + 1
                        else
                            n_bc = n_bc + 1
                        end if
                    end if
                end do
                !>
                !>
                do nnn = 1,num1to1
                    !> when the block and iblock is to the same block
                    !>
                    if(donorname(nnn) .eq. zonename(n))then
                        !>
                        !>
                        mesh%bcs(n_bc)%block_index    = nb
                        !>
                        do nn=1,nobl
                            if(znname(nnn).eq. zonename(nn) )then
                                if(mgflags .ne. 0)then
                                    mesh%bcs(n_bc)%iblock_index = (nn-1)*(global_level-1) + nn
                                else
                                    mesh%bcs(n_bc)%iblock_index = nn
                                end if
                            end if
                        end do
                        !>
                        !>
                        mesh%bcs(n_bc)%bc_type        = 'cut'
                        !>
                        !>
                        mesh%bcs(n_bc)%istart         = idonor_range(1,nnn) !mesh%bcs(nnn)%istart_cut
                        mesh%bcs(n_bc)%jstart         = idonor_range(2,nnn) !mesh%bcs(nnn)%jstart_cut
                        mesh%bcs(n_bc)%kstart         = idonor_range(3,nnn) !mesh%bcs(nnn)%kstart_cut
                        !>
                        !>
                        mesh%bcs(n_bc)%iend           = idonor_range(4,nnn) !mesh%bcs(nnn)%iend_cut
                        mesh%bcs(n_bc)%jend           = idonor_range(5,nnn) !mesh%bcs(nnn)%jend_cut
                        mesh%bcs(n_bc)%kend           = idonor_range(6,nnn) !mesh%bcs(nnn)%kend_cut
                        !>
                        !>
                        mesh%bcs(n_bc)%istart_cut     = irange(1,nnn) !mesh%bcs(nnn)%istart
                        mesh%bcs(n_bc)%jstart_cut     = irange(2,nnn) !mesh%bcs(nnn)%jstart
                        mesh%bcs(n_bc)%kstart_cut     = irange(3,nnn) !mesh%bcs(nnn)%kstart
                        !>
                        !>
                        mesh%bcs(n_bc)%iend_cut       = irange(4,nnn) !mesh%bcs(nnn)%iend
                        mesh%bcs(n_bc)%jend_cut       = irange(5,nnn) !mesh%bcs(nnn)%jend
                        mesh%bcs(n_bc)%kend_cut       = irange(6,nnn) !mesh%bcs(nnn)%kend
                        !>
                        !>
                        if( abs(itransform(1,nnn)) .eq. 1) mesh%bcs(n_bc)%irank = 1
                        if( abs(itransform(1,nnn)) .eq. 2) mesh%bcs(n_bc)%jrank = 1
                        if( abs(itransform(1,nnn)) .eq. 3) mesh%bcs(n_bc)%krank = 1
                        !>
                        !>
                        if( abs(itransform(2,nnn)) .eq. 1) mesh%bcs(n_bc)%irank = 2
                        if( abs(itransform(2,nnn)) .eq. 2) mesh%bcs(n_bc)%jrank = 2
                        if( abs(itransform(2,nnn)) .eq. 3) mesh%bcs(n_bc)%krank = 2
                        !>
                        !>
                        !>
                        if( abs(itransform(3,nnn)) .eq. 1) mesh%bcs(n_bc)%irank = 3
                        if( abs(itransform(3,nnn)) .eq. 2) mesh%bcs(n_bc)%jrank = 3
                        if( abs(itransform(3,nnn)) .eq. 3) mesh%bcs(n_bc)%krank = 3
!                        write(161,'("n_bc=",10i4)') n_bc, mesh%bcs(n_bc)%istart, mesh%bcs(n_bc)%iend,&
!                         mesh%bcs(n_bc)%jstart, mesh%bcs(n_bc)%jend, mesh%bcs(n_bc)%kstart, mesh%bcs(n_bc)%kend,&
!                         mesh%bcs(n_bc)%irank,mesh%bcs(n_bc)%jrank,mesh%bcs(n_bc)%krank
                        !>
                        !>
#if defined DEBUG
                        open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                        write(61,'(" block_id bc_num bc_type   is   ie   js   je   ks   ke cut_block   is   ie   js   je   ks   ke ")')
                        write(61,'(i9,i7,1x,a6,6i5,i11,6i5)')  &
                                     mesh%bcs(n_bc)%block_index,&
                                     n_bc,&
                                     mesh%bcs(n_bc)%bc_type,&
                                     mesh%bcs(n_bc)%istart,&
                                     mesh%bcs(n_bc)%iend,&
                                     mesh%bcs(n_bc)%jstart,&
                                     mesh%bcs(n_bc)%jend,&
                                     mesh%bcs(n_bc)%kstart,&
                                     mesh%bcs(n_bc)%kend,&
                                     mesh%bcs(n_bc)%iblock_index,&
                                     mesh%bcs(n_bc)%istart_cut,&
                                     mesh%bcs(n_bc)%iend_cut,&
                                     mesh%bcs(n_bc)%jstart_cut,&
                                     mesh%bcs(n_bc)%jend_cut,&
                                     mesh%bcs(n_bc)%kstart_cut,&
                                     mesh%bcs(n_bc)%kend_cut!,n2p(mesh%bcs(n_bc)%block_index),n2p(mesh%bcs(n_bc)%iblock_index)
#endif
                        !> set the multigrid methods bc infotmations
                        !> the blocks = iblock
                        !> set the ijk-dir dimensions
                        if(mgflags .ne. 0)then
                            do nccg=1,global_level-1
                                !>
                                !>
                                n_bc = n_bc  + 1
                                !>
                                !>
                                mesh%bcs(n_bc)%block_index   = mesh%bcs(n_bc-1)%block_index + 1
                                mesh%bcs(n_bc)%iblock_index  = mesh%bcs(n_bc-1)%iblock_index + 1
                                !>
                                !>
                                !>
                                mesh%bcs(n_bc)%bc_type  = mesh%bcs(n_bc-1)%bc_type
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%istart .eq. 1)then
                                   mesh%bcs(n_bc)%istart   = 1
                                else
                                   mesh%bcs(n_bc)%istart   = mesh%bcs(n_bc-1)%istart/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jstart .eq. 1)then
                                   mesh%bcs(n_bc)%jstart   = 1
                                else
                                   mesh%bcs(n_bc)%jstart   = mesh%bcs(n_bc-1)%jstart/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kstart .eq. 1)then
                                   mesh%bcs(n_bc)%kstart   = 1
                                else
                                   mesh%bcs(n_bc)%kstart   = mesh%bcs(n_bc-1)%kstart/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%iend .eq. 1)then
                                    mesh%bcs(n_bc)%iend   = 1
                                else
                                    mesh%bcs(n_bc)%iend   = mesh%bcs(n_bc-1)%iend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jend .eq. 1)then
                                    mesh%bcs(n_bc)%jend  = 1
                                else
                                    mesh%bcs(n_bc)%jend   = mesh%bcs(n_bc-1)%jend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kend .eq. 1)then
                                    mesh%bcs(n_bc)%kend  = 1
                                else
                                    mesh%bcs(n_bc)%kend   = mesh%bcs(n_bc-1)%kend/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%istart_cut .eq. 1)then
                                    mesh%bcs(n_bc)%istart_cut  = 1
                                else
                                    mesh%bcs(n_bc)%istart_cut   = mesh%bcs(n_bc-1)%istart_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jstart_cut .eq. 1)then
                                    mesh%bcs(n_bc)%jstart_cut  = 1
                                else
                                    mesh%bcs(n_bc)%jstart_cut   = mesh%bcs(n_bc-1)%jstart_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kstart_cut .eq. 1)then
                                    mesh%bcs(n_bc)%kstart_cut  = 1
                                else
                                    mesh%bcs(n_bc)%kstart_cut   = mesh%bcs(n_bc-1)%kstart_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%iend_cut .eq. 1)then
                                    mesh%bcs(n_bc)%iend_cut  = 1
                                else
                                    mesh%bcs(n_bc)%iend_cut   = mesh%bcs(n_bc-1)%iend_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%jend_cut .eq. 1)then
                                    mesh%bcs(n_bc)%jend_cut  = 1
                                else
                                    mesh%bcs(n_bc)%jend_cut   = mesh%bcs(n_bc-1)%jend_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                if(mesh%bcs(n_bc-1)%kend_cut .eq. 1)then
                                    mesh%bcs(n_bc)%kend_cut  = 1
                                else
                                    mesh%bcs(n_bc)%kend_cut   = mesh%bcs(n_bc-1)%kend_cut/2 + 1
                                end if
                                !>
                                !>
                                !>
                                mesh%bcs(n_bc)%irank = mesh%bcs(n_bc-1)%irank
                                mesh%bcs(n_bc)%jrank = mesh%bcs(n_bc-1)%jrank
                                mesh%bcs(n_bc)%krank = mesh%bcs(n_bc-1)%krank

#if defined DEBUG
                                open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                                write(61,'("cut_block bc_num bc_type   is   ie   js   je   ks   ke block_id   is   ie   js   je   ks   ke ")')
                                write(61,'(i9,i7,1x,a6,6i5,i10,6i5)') &
                                             mesh%bcs(n_bc)%block_index,&
                                             n_bc,&
                                             mesh%bcs(n_bc)%bc_type,&
                                             mesh%bcs(n_bc)%istart,&
                                             mesh%bcs(n_bc)%iend,&
                                             mesh%bcs(n_bc)%jstart,&
                                             mesh%bcs(n_bc)%jend,&
                                             mesh%bcs(n_bc)%kstart,&
                                             mesh%bcs(n_bc)%kend,&
                                             mesh%bcs(n_bc)%iblock_index,&
                                             mesh%bcs(n_bc)%istart_cut,&
                                             mesh%bcs(n_bc)%iend_cut,&
                                             mesh%bcs(n_bc)%jstart_cut,&
                                             mesh%bcs(n_bc)%jend_cut,&
                                             mesh%bcs(n_bc)%kstart_cut,&
                                             mesh%bcs(n_bc)%kend_cut!,n2p(mesh%bcs(n_bc)%block_index),n2p(mesh%bcs(n_bc)%iblock_index)
#endif
                            end do
                            n_bc = n_bc + 1
                        else
                            n_bc = n_bc  + 1
                        end if
                    end if
                enddo
            enddo
            num_bc = n_bc - 1
            !check boundary ocnditions

            deallocate(irange)
            deallocate(idonor_range)
            deallocate(znname)
            deallocate(donorname)
            deallocate(itransform)
            deallocate(connectname)
            deallocate(zonename)
#if defined PMPI
        endif
		!>
		!>
		!>
        call mpi_bcast(num_bc,1,mpi_integer,myhost,mycomm,ierr)
        if(myid .ne. myhost)then
            allocate( mesh%bcs(num_bc)) !mesh%bcs=>null
        endif
        !>
        !>
        !>
#if defined DEBUG
        write(*,'("bcast myid=",i6," num_bc=",i8)') myid,num_bc
#endif
        call mpi_barrier(mycomm,ierr)
#endif
        return
    end subroutine get_cgns_blocks_boundary
