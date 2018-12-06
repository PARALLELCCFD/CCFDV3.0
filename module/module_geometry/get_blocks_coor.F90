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
    subroutine get_blocks_coor(nbl,iccg)
		!>
		!>
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use cell_module
        use coordinate_module
        use nodes_var
        use nodes_paras
		!>
		!>
!		implicit none
        include "cgnslib_f.h"
#if defined PMPI
        include "mpif.h"
        !integer::istat(mpi_status_size)
        !real,dimension(:),allocatable::nobl_send_xyz
#endif
        integer::isize(3,3),irmin(3),irmax(3),nn
        integer::idimc,jdimc,kdimc,xyz_buf,mytag,nd_dest
        integer::iccg,ibase,igrid,idouble,icoord,stats
        integer::nbl1,nbl!,ifi,ifj
        integer::i,j,k,icor,jcor,kcor
        character*32 zonename,coordname,testname(3)
        real(kind = dprec),dimension(:,:,:),allocatable:: xc,yc,zc
!        real(kind = dprec):: xc(icor,jcor,kcor),yc(icor,jcor,kcor),zc(icor,jcor,kcor)
        !>
        !>
        type(blocks_type)    ,pointer :: iblock
        type(overlap_type)   ,pointer :: mesh

        !> allocate the array for xyz dimension
        !>
        mesh   => grids(imesh)
        iblock => mesh%blocks(nbl)
        !>
        !>
        icor =  iblock%idim
        jcor =  iblock%jdim
        kcor =  iblock%kdim
        !>
        !>
        !>
        allocate(xc(icor,jcor,kcor),yc(icor,jcor,kcor),zc(icor,jcor,kcor))
#if defined PMPI
        if(myid == myhost)then

#endif
            !   coordinate names that we are looking for:
            testname(1)='CoordinateX'
            testname(2)='CoordinateY'
            testname(3)='CoordinateZ'
            ibase = 1
#if defined cgns_single
            idouble = 0
#else
            idouble = 1
#endif
            !> the mg blocks and original blocks
            !> and allocate the original blocks at the top level of mg
            !> so must computing the num of original blcoks to
            !> read the xyz coordinate of cgns grid files
            !>
            if(mgflags .ne. 0)then
                nbl1 = nbl - (nbl/global_level)*(global_level-1)
                open(unit=61,file='ccfd.out',form='formatted',status='unknown')
                write(61,'("original blocks=",i4," mg blocks=",i4)') nbl1,nbl
            else
                nbl1 = nbl
            end if

            !>   read general zone information
            call cg_zone_read_f(iccg, ibase, nbl1, zonename, isize, ier)
            if (ier .ne. 0) call cg_error_exit_f

            if(isize(1,1) .ne. iblock%idim .or. isize(2,1) .ne. iblock%jdim .or.&
                isize(3,1) .ne. iblock%kdim) then
                write(*,'('' grid index inconsistencies:  isize='',3i5,&
                    ''idim,jdim,kdim='',3i5)') isize(1,1),isize(2,1),isize(3,1),&
                    iblock%idim,iblock%jdim,iblock%kdim
                write(*,'('' be sure to order the zones alphabetically'',&
                    '' in the input file!'')')
                stop
            end if

            !>   find out how many grid coordinates exist
            call cg_ncoords_f(iccg, ibase, nbl1, ncoords, ier)
            if (ier .ne. 0) call cg_error_exit_f
            if (ncoords .ne. 3) then
                write(*,'('' ncoords='',i5,''.  expecting 3.'')') ncoords
                stop
            end if
            !>   check coordinate names in data base
            do icoord=1,3
                call cg_coord_info_f(iccg, ibase, nbl1, icoord, itype,&
                coordname,ier)
                if (ier .ne. 0) call cg_error_exit_f
                if( coordname .eq. testname(1) .or.&
                    coordname .eq. testname(2) .or.&
                    coordname .eq. testname(3)) then
                    continue
                else
!                    ifi = integer(testname(1))
!                    ifj = integer(coordname)
!                    if(coordname .ne. testname(1)) write(*,*) ifi,ifj
                    write(*,'('' coordname of '',a32,'' unrecognized.'', "the coordname of test",a32)') coordname,testname(icoord)
                    write(*,'('' looking for CoordinateX, CoordinateY, and'',&
                            '' CoordinateZ'')')
                    stop
                end if
            enddo
            !>   set up array bounds:
            call cg_zone_read_f(iccg, 1, nbl1, zonename, isize, ier)
            !write(*,'("nm=",i4," isize=",i4," jsize=",i4," ksize=",i4," im=",i4," jm=",i4," km=",i4)') nbl,isize(1,1),isize(2,1),isize(3,1),dat(nbl)%im +1,dat(nbl)%jm +1,dat(nbl)%km +1
            irmin(1) = 1
            irmin(2) = 1
            irmin(3) = 1

            irmax(1) = icor
            irmax(2) = jcor
            irmax(3) = kcor

            !>   read x-coordinate
            if (idouble .eq. 1) then
                !>ps:the dimension of xc or dat(nbl)%x must mach with the irmin and irmax
                !>if the dimension of xc or dat(nbl)%x more than the irmin and irmax
                !>read the xyz-coordinate will have fault
                !>
                call cg_coord_read_f(iccg, ibase, nbl1, 'CoordinateX',&
                           realdouble, irmin, irmax, xc,ier)!dat(nbl)%x, ier)

            else
                call cg_coord_read_f(iccg, ibase, nbl1, 'CoordinateX',&
                           realsingle, irmin, irmax, xc,ier)!dat(nbl)%x, ier)
            end if
            !>if have any error and exit the subroutine
            if (ier .ne. 0) call cg_error_exit_f
            !>read y-coordinate
            if (idouble .eq. 1) then
                call cg_coord_read_f(iccg, ibase, nbl1, 'CoordinateY',&
                           realdouble, irmin, irmax, yc,ier)!dat(nbl)%y, ier)
            else
                call cg_coord_read_f(iccg, ibase, nbl1, 'CoordinateY',&
                           realsingle, irmin, irmax, yc,ier)!dat(nbl)%y, ier)
            end if
            !>
            if (ier .ne. 0) call cg_error_exit_f
            !> read z-coordinate
            if (idouble .eq. 1) then
                call cg_coord_read_f(iccg, ibase, nbl1, 'CoordinateZ',&
                           realdouble, irmin, irmax, zc,ier)!dat(nbl)%z, ier)
            else
                call cg_coord_read_f(iccg, ibase, nbl1, 'CoordinateZ',&
                           realsingle, irmin, irmax, zc,ier)!dat(nbl)%z, ier)
            end if
            !>
            if (ier .ne. 0) call cg_error_exit_f
            !>
            !deallocate(x,y,z)
!            call cg_close_f(iccg,ier)
!            if (ier .ne. 0) call cg_error_exit_f
            !> restore the xyz-coordinate
            !>

            if(ialph .eq. 0)then
                !> in the x-z plane(with z "up")
                do i=icor+1,2,-1
                    do j=jcor+1,2,-1
                        do k=kcor+1,2,-1
                            iblock%coordinate(i,j,k)%x = xc(i-1,j-1,k-1)
                            iblock%coordinate(i,j,k)%y = yc(i-1,j-1,k-1)
                            iblock%coordinate(i,j,k)%z = zc(i-1,j-1,k-1)
!                            write(118,'(3i4,3e24.16)')i,j,k,iblock%coordinate(i,j,k)%x,iblock%coordinate(i,j,k)%y,iblock%coordinate(i,j,k)%z
                        end do
                    end do
                end do
            else if(ialph .eq. 1)then
                !>in the x-y plane(with y "up")
                do i=icor+1,2,-1
                    do j=jcor+1,2,-1
                        do k=kcor+1,2,-1
                            iblock%coordinate(i,j,k)%x =  xc(i-1,j-1,k-1)
                            iblock%coordinate(i,j,k)%y =  -zc(i-1,j-1,k-1)
                            iblock%coordinate(i,j,k)%z =  yc(i-1,j-1,k-1)
                        end do
                    end do
                end do
            end if
!            call cpu_time(finish_t)
!            write(*,'("in the block=",i4," and the time spent3 =",e20.14)') nbl,finish_t-start_t


#if defined PMPI

            if(nobl < nodes)then
                write(*,*) 'ccfdv3.0 err0r:: the blocks you read from grid files must larger than porcessoers'
                write(*,*) 'plese check you grid or set a small processor'
                stop
            endif

        endif !myhost

#endif
        deallocate(xc)
        deallocate(yc)
        deallocate(zc)
        return
    end subroutine get_blocks_coor

