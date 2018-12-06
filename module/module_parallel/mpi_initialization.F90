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
  subroutine mpi_initialization()

    !use nodes_var
    use nodes_paras
    !use nodes_var_bc
!
#if defined PMPI
    include "mpif.h"
!***********************************************************************
!
!     function: initializes for parallel (mimd type) computations
!               determines rank of processor
!
!     inputs:    name    type           description
!                pglob   data_set       pointer to all data
!
!     outputs :  name    type           description
!
!
!     called from: amain
!
!     creation by: peter eliasson
!
!     creation date: 1999-12-17
!
!     modifications:
!     date        vers   programmer    description
!
!     references:
!
!----------------------------------------------------------------------
!
!     arguments:
!
!
!     subprograms:
!
!
!     local variables:
!!
!    integer ierror
!------------the type of the actual argument differs from the type of the dummy argument.-----------------------------------------------------------
!     body of mpi_init
!-----------------------------------------------------------------------
!

    call mpi_init(ierror)
    if (ierror/=mpi_success) then
      write(*,9010)' error from mpi_init, ierror=',ierror
      stop
      !call mpi_error_string(ierror,errmsg,len,ier)
    end if
    call mpi_comm_rank(mpi_comm_world,myid,ierror)
    if (ierror/=0) then
      write(*,9010)' error from mpi_comm_rank, ierror=',ierror
      stop
    end if
    call mpi_comm_size(mpi_comm_world,numprocs,ierror)
    if (ierror/=0) then
      write(*,9010)' error from mpi_comm_size, ierror=',ierror
      stop
    end if
    mycomm = mpi_comm_world
    myhost     = 0
    istat_size = mpi_status_size
    if(numprocs .lt. 2)then
        write(*,100)
100     format('stopping the ccfdv3.0 ...... only running a host process;'&
              ' you must run the mpi  versions with more than one process')
        call mpi_abort(mpi_comm_world,myid,mpi_error)
        call mpi_finalize(ierr)

    end if
    nodes      = numprocs -1
!
9010 format(a,i5)
!
#endif
    return
!-------------------------------------end of mpi_init  -------------------
   end subroutine mpi_initialization
!
