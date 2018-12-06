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
	subroutine get_grid_num()
        !>
        !>
        !>
        use global_parameter
        use mesh_overlap_module
!        use blocks_module
        use nodes_paras
        use nodes_var
        use nodes_mg
		!>
		!>
		implicit none

        include "cgnslib_f.h"
#if defined PMPI
        include "mpif.h"
       !integer,dimension::istat(mpi_status_size)
#endif
        integer isize(3,3)
        integer:: ibase,nbl,idimc,jdimc,kdimc
        integer:: ncg_flag,stats,nm,grid,iccg
        integer::ier,nbases,n,nn_grid,ierr
        logical:: alive
        character(len=32)::zonename
        integer::zonetype
        !>
        !>
        type(overlap_type),pointer :: mesh
#if defined PMPI
        real(kind=dprec),dimension(:),allocatable:: nobl_send_xyz
#endif
        !>
        !>
        !>
        imesh = 1
        mesh => grids(imesh)
#if defined PMPI
        if(myid == myhost)then
#endif
            !>read cgns file
            inquire(file=grid_files,exist=alive)
            if(.not.alive) then
                write(*,*) 'error: gird file is not exist'
                stop
            endif
            open (11,file =grid_files,status='old')

            call cg_open_f(grid_files,mode_read,iccg,ier)

            if (ier .ne. 0) call cg_error_exit_f
            call cg_nbases_f(iccg,nbases,ier)
            if (ier .ne. 0) call cg_error_exit_f
            if (nbases .ne. 1) then
                write(*,'('' error nbases ='',i5,'' (should be 1)'')') nbases
                write(*,'('' stopping'')')
                stop
            end if
            ibase=1

            call cg_nzones_f(iccg,ibase,nobl,ier)


            !	  write(iunit11,*) 'come to the cgns_reader33'
            if (ier .ne. 0) call cg_error_exit_f
#if defined DEBUG
            open(unit=61,file='ccfd.out',form='formatted',status='unknown')
            write(61,'('' number of blocks from grid files='',i6)') nobl
            write(61,*)'******************************************'
#endif
        !****************************************************
            if(mgflags .ne. 0)then
                nblocks = global_level * nobl
#if defined DEBUG
                open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                write(61,'(" number of blocks about multigrid methods= ",i6)') nblocks
                write(61,*)'******************************************'
#endif
            else
                nblocks = nobl
#if defined DEBUG
                open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                write(61,'(" number of blocks about no-multigrid methods= ",i6)') nblocks
                write(61,*)'******************************************'
#endif
            endif
       !****************************************************
#if defined PMPI
        endif
        call mpi_bcast(nobl,1,mpi_integer,myhost,mycomm,ierr)
        call mpi_bcast(nblocks,1,mpi_integer,myhost,mycomm,ierr)
#endif
        call allocate_nodes_paras(nobl)
#if defined PMPI
        if(myid == myhost)then
#endif
            !>
            !>
            allocate(mesh%blocks(nblocks))
            !>
            !>



            !****************************************************
            nn_grid = 0
#if defined DEBUG
            write(61,*) 'the i,j,k dimension read form grid files'
            write(61,*) '****************************************'
            write(61,*) 'block_id  idim  jdim  kdim    zonename'
#endif
            do n=1,nobl
                !   read general zone information
                call cg_zone_type_f(iccg,ibase,n,zonetype,ier)
                if(zonetype .ne. structured)then
                    write(*,*) 'error! the grid is not a structured mesh!'
                    stop
                end if
                call cg_zone_read_f(iccg, ibase, n, zonename, isize, ier)
                !write(61,*) iccg,ibase,n,ier
                if (ier .ne. 0) call cg_error_exit_f
                nxt(n)=isize(1,1)
                nyt(n)=isize(2,1)
                nzt(n)=isize(3,1)
                nn_grid = nn_grid + nxt(n)*nyt(n)*nzt(n)
#if defined DEBUG
                write(61,'(i8,2x,3i6,6x,a20)') n,nxt(n),nyt(n),nzt(n),zonename
#endif
            enddo
            !>
            !>
#if defined DEBUG
            open(unit=61,file='ccfd.out',form='formatted',status='unknown')
            write(61,*) '*********************************************'
            write(61,*) 'ccfd.3.0 : mesh num of grid = ',nn_grid
            write(61,*) '*********************************************'
#endif
            write(*,*) '*********************************************'
            write(*,*) 'ccfd.3.0 : mesh num of grid = ',nn_grid
            write(*,*) '*********************************************'
       !****************************************************
            num_grid = 0
            ijkmax    = 0
            nbl      = 1
            !create the multi_grid methods blocks,coarse mesh and fine mesh
            if (mgflags .ne. 0) ncg_flag = (global_level-1)*2
            do nm=1,nobl
                idimc = nxt(nm)!-1
                jdimc = nyt(nm)!-1
                kdimc = nzt(nm)!-1
                if(mgflags .ne. 0)then
                    if(mod(idimc-1,ncg_flag) .ne. 0 .or. mod(jdimc-1,ncg_flag) .ne. 0 .or. mod(kdimc-1,ncg_flag) .ne. 0)then
                        write(*,*) '***an error has occurred about level of mg and grid dimension***'
                        write(*,*) '            '
                        write(*,*) 'when the level of multigrid = ',global_level
                        write(*,*) 'such as x-dimension must to be divisible by ',ncg_flag
                        write(*,*) 'and your xyz-dimension = ',idimc,jdimc,kdimc
                        write(*,*) 'please reset the level of mg or input a new mesh '
                        write(*,*) '*****************************************************************'
                        stop
                    end if
                    num_grid  = num_grid + (nxt(nm)-1)*(nyt(nm)-1)*(nzt(nm)-1)
                    ijkmax      = max0(ijkmax,nxt(nm),nyt(nm),nzt(nm))+2
                    call  mem_allocate_block(mesh%blocks(nbl),idimc,jdimc,kdimc)
#if defined DEBUG
                    write(61,'("block_id=",i6," and the dimensions of  i=",i6," j=",i6," k=",i6)') nbl,idimc,jdimc,kdimc
#endif
                    do ncg = 1,global_level-1
                        nbl   = nbl + 1
                        idimc = idimc/2 + 1
                        jdimc = jdimc/2 + 1
                        kdimc = kdimc/2 + 1
                        call  mem_allocate_block(mesh%blocks(nbl),idimc,jdimc,kdimc)
#if defined DEBUG
                        write(61,'("block_id=",i6," and the dimensions of  i=",i6," j=",i6," k=",i6)') nbl,idimc,jdimc,kdimc
#endif
                    end do
                    nbl = nbl + 1
                else
                    num_grid  = num_grid + (nxt(nm)-1)*(nyt(nm)-1)*(nzt(nm)-1)
                    ijkmax      = max0(ijkmax,nxt(nm),nyt(nm),nzt(nm))+2
#if defined DEBUG
                    write(61,'("block_id=",i6," and the dimensions of  i=",i6," j=",i6," k=",i6)') nm,idimc,jdimc,kdimc
                    write(*,'("block_id=",i6," and the dimensions of  i=",i6," j=",i6," k=",i6)') nm,idimc,jdimc,kdimc
#endif
                    call  mem_allocate_block(mesh%blocks(nm),idimc,jdimc,kdimc)
                    !>
                    !>
                    !>
                endif
            end do
#if defined PMPI
        endif
#endif
        call allocate_nodes_mg(nblocks)
        call allocate_nodes_var(nblocks)

#if defined PMPI
        call mpi_bcast(num_grid,1,mpi_integer,myhost,mycomm,ierr)
        call mpi_bcast(ijkmax,1,mpi_integer,myhost,mycomm,ierr)
        call mpi_bcast(nxt,nobl,mpi_integer,myhost,mycomm,ierr)
        call mpi_bcast(nyt,nobl,mpi_integer,myhost,mycomm,ierr)
        call mpi_bcast(nzt,nobl,mpi_integer,myhost,mycomm,ierr)

        if(myid == myhost)then
#endif
            call cg_close_f(iccg,ier)
            if (ier .ne. 0) call cg_error_exit_f
#if defined PMPI
        endif
#endif
!
        return
        !>
    end subroutine get_grid_num
!
!*******************************************************************

