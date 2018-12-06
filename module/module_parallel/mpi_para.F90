module nodes_var

!   the dimenson for parallel and  keep the information about mpi's tag!
	integer,dimension(:),allocatable :: n2p
end module nodes_var
!
subroutine allocate_nodes_var(nodes)
      use nodes_var
	  implicit none
	  integer nodes,nerror
	  nerror = 0
	  allocate(n2p(nodes),stat = nerror)
      if(nerror/=0) then
         write(*,*) 'error......'
         stop
      end if
end subroutine allocate_nodes_var
!*********************************************************
module nodes_mg

	integer,dimension(:),allocatable :: level_mg
end module nodes_mg
!
subroutine allocate_nodes_mg(nblocks)
      use nodes_mg
	  implicit none
	  integer nblocks,nerror
	  nerror = 0
	  allocate(level_mg(nblocks),stat = nerror)
      if(nerror/=0) then
         write(*,*) 'error......'
         stop
      end if
end subroutine allocate_nodes_mg
!*********************************************************
module nodes_paras
     implicit none
     integer:: myid,numprocs,myhost,mycomm,nodes,istat_size
	 integer,dimension(:),allocatable :: nxt
	 integer,dimension(:),allocatable :: nyt
	 integer,dimension(:),allocatable :: nzt
end module nodes_paras
subroutine allocate_nodes_paras(nobl)
      use nodes_paras
	  implicit none
	  integer::nobl,nerror
	  nerror  = 0
	  allocate(nxt(nobl),stat=nerror)
	  if(nerror /= 0)then
	        write(*,*) 'ccfd3.0:error.......'
	  end if
	  allocate(nyt(nobl),stat=nerror)
	  if(nerror /= 0)then
	        write(*,*) 'ccfd3.0:error.......'
	  end if
	  allocate(nzt(nobl),stat=nerror)
	  if(nerror /= 0)then
	        write(*,*) 'ccfd3.0:error.......'
	  end if
end subroutine allocate_nodes_paras
!
!*************************************************************
module nodes_var_bc
!    implicit none
!   the dimenson for parallel and  keep the information about mpi's tag!
	real*8,dimension(:),allocatable :: wk
	integer,dimension(:),allocatable :: ireq_r,ireq_s,keep1,keep2,index_r
!    integer:: i_size
end module nodes_var_bc
!
subroutine allocate_nodes_var_bc(num_bc,i_size1,level)
!      use global_bc
      use nodes_var_bc
	  implicit none
	  integer :: num_bc,nerror,i_size1,i,level
      nerror = 0
      i = i_size1*level
      allocate(wk(i_size1*level),stat=nerror)
      if(nerror /= 0)then
         write(*,*) 'ccfd:error1.......'
         stop
      end if
      wk(1:i)=0.0
      !>
      i = num_bc*level
      allocate(ireq_r(num_bc*level),stat=nerror)
      if(nerror/=0) then
         write(*,*) 'ccfd:error2......'
         stop
      end if
      ireq_r(1:i) = 0
      !>
      allocate(ireq_s(num_bc*level),stat=nerror)
      if(nerror/=0) then
         write(*,*) 'ccfd:error3......'
         stop
      end if
      ireq_s(1:i) = 0
      !>
      allocate(keep1(num_bc*level),stat=nerror)
      if(nerror/=0) then
         write(*,*) 'ccfd:error4......'
         stop
      end if
      keep1(1:i) =0
      !>
      allocate(keep2(num_bc*level),stat=nerror)
      if(nerror/=0) then
         write(*,*) 'ccfd:error5......'
         stop
      end if
      keep2(1:i) = 0
      !>
      allocate(index_r(num_bc*level),stat=nerror)
      if(nerror/=0) then
         write(*,*) 'ccfd:error6......'
         stop
      end if
      index_r(1:i) = 0
      !>
end subroutine allocate_nodes_var_bc
subroutine deallocate_nodes_var_bc()
      use nodes_var_bc
	  implicit none
      deallocate(ireq_r)
      deallocate(ireq_s)
      deallocate(keep1)
      deallocate(keep2)
      deallocate(index_r)
      deallocate(wk)
end subroutine deallocate_nodes_var_bc
!*************************************************************
!*************************************************************
