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
	subroutine  steady_cfd

		use  mesh_overlap_module
		use  global_parameter
		use  nodes_paras
		use  nodes_var
		use  nodes_paras
		use  nodes_var_bc
		use  bc_module
		!>
		!>
		implicit none
#if defined PMPI
		include "mpif.h"
#endif
        type(bc_types),pointer :: ibc
        !>
        !>
        integer :: ierr,min_flags,ii
        !>
        !>
        min_flags  = 1
		!*********************************************************
		!>
		!>
		if(myid == myhost)then
			write(*,*) '**************************************************'
			write(*,*) '      ccfd.3.0 :  read input files ......    '
			write(*,*) '**************************************************'
		end if
		!>
		!>
        !>read in the control parameters and flow parameters
        call input_global
        !>
        !>
		!select case(mesh_type)
		select case('cgns')
		case('cgns')
			call input_grids_cgns()
		case('plot3d')
			!call input_files_plot3d()
		case('overlap')
			!call input_files_overlap()
		end select
		!>
#if defined PMPI
		call mpi_barrier(mycomm,ierr)
		if(myid .eq. myhost)then
#endif
		   write(*,*) '**************************************************'
		   write(*,*) '      ccfd.3.0 :  computing metric ......    '
		   write(*,*) '**************************************************'
#if defined PMPI
		end if
#endif
		call metric
		!>
		!>
		!>
		!> evaluate directed distance from wall face to inner cell
		!> for use in evaluating the baldwin-lomax turbulence model
		if(ivisc_i .eq. 2 .or. ivisc_j .eq. 2 .or. ivisc_k .eq. 2 )then
#if defined PMPI
            if(myid .eq. myhost)then
#endif
                write(*,*) '**************************************************'
                write(*,*) 'ccfd.3.0:computing distance for B-L model.....'
                write(*,*) '**************************************************'
#if defined PMPI
			end if
#endif
            !>
            !>
			call find_distance_bl
			!>
			!>
		elseif(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3 )then
#if defined PMPI
            if(myid .eq. myhost)then
#endif
                write(*,*) '**************************************************'
                write(*,*) 'ccfd.3.0:computing  distance for SA/SST model.....'
                write(*,*) '**************************************************'
#if defined PMPI
			end if
#endif
            !>
            !>
            !>if no walls at all, don't do min distance calculation
            do ii = 1,num_bc
                ibc => grids(imesh)%bcs(ii)
                if(ibc%bc_type .eq. 'wallviscid' )then
                    min_flags  = 1
                    exit
                end if
                min_flags  = 0
            end do
            ibc => null()
            !>
            !>
            !>
			if( min_flags .ne. 0)then
                !>
                !>
                call find_distance_min
            end if
			!>
			!>
        end if
		!>
		!>
		!> initial the turbulent initial conditions
		!>
		if(ivisc_i .ge. 2 .or. ivisc_j .ge. 2 .or. ivisc_k .ge. 2)then
#if defined PMPI
            if(myid .eq. myhost)then
#endif
                write(*,*) '*********************************************'
                write(*,*) 'ccfd.3.0:initial the turbulent conditions ...'
                write(*,*) '*********************************************'
#if defined PMPI
			end if
#endif
            call turbulent_init
        end if
		!>
		!>
#if defined PMPI

		call mpi_barrier(mycomm,ierr)
		!>
		!>
		!>
		if(myid .eq. myhost)then
#endif
			write(*,*) '**************************************************'
			write(*,*) '      ccfd.3.0 :  begin to compute ......    '
			write(*,*) '**************************************************'
#if defined PMPI
		end if
		if(myid .eq. myhost)then
#endif
			open(unit=71,file='ccfd.resid',form='formatted',status='unknown')
			write(71,'("iteration",10x,"res ",21x,"cl",22x,"cd",21x,"time")')
#if defined PMPI
		endif
#endif
		!>
		!> set the mesh sequence methods
		!> the level of multigrid and mesh sequence's methods
		!> must appropriate
		!>
		if(meshseqflags .ne. 0)then
			if(mgflags .eq.  0)then
#if defined PMPI
				if(myid .eq. myhost)then
#endif
					write(*,*) 'have a wrong parameter of multigrid'
					write(*,*) 'when used the mesh sequence must set'
					write(*,*) 'the mgflags no zero '
					stop
#if defined PMPI
				endif
#endif
			end if
			if(meshseque .gt. global_level)then
#if defined PMPI
				if(myid .eq. myhost)then
#endif
				write(*,*)"ccfd3.0:warning!"
                write(*,*)"the mesh sequence's level must "
                write(*,*)"be less or the level of multigrid,so we will "
                write(*,*)"modify the level of meshsequence used.and the level"
                write(*,*)"meshseque = global_level"
                meshseque = global_level
#if defined PMPI
				end if
#endif
			end if
		elseif(meshseqflags .eq. 0)then
            meshseque = 1
		end if
		!>
		!>
        !> call the steady subroutine
        !> computing the steady problem
        do  imeshseque= 1,meshseque
            if(myid .eq. myhost)then
                write(*,*) '***************************************************************'
                write(*,*) 'the ',imeshseque,'-level  of mesh sequence iterative will begin'
                write(*,*) '***************************************************************'
            end if
            !> interpolate from coarse mesh to finer mesh
            !>
            !>
            if(imeshseque .gt.1 .and. imeshseque .le. meshseque )then
                mode  = 0
                if(myid .eq. myhost)then
                    write(*,*) '****************************************************************'
                    write(*,*) 'used the full multigrid methods and prolonged to the finer mesh '
                    write(*,*) '****************************************************************'
                end if
                call interpolate_solution
            end if
            !>
            !>
            call steady()
        end do
        !>
        !>
        !>
	end subroutine  steady_cfd
