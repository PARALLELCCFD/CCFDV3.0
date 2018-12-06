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
!>  subroutine load_balance
!>  use the k-way domain decomposition methods
!>  to decomposition the weighted graph of blocks
!>  and the weighted can be the grids or the face's area
!>  call the METIS API to domain decomposition the graph.
!>  last edit 2018-05-24
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
    subroutine load_balance
        !>
        !>
        !>
        use global_parameter
        use mesh_overlap_module
        use bc_module
        use blocks_module
        use nodes_var
        use nodes_paras
        !>
        !>
        !>
		implicit none
		!>
		!>
		!>
        integer :: load_flags,nccg,ll,lll
        integer :: objval,nparts,nvtxs
        integer :: nm,nbl
        integer :: num_of_bc,mm,dim1,dim2,dim3
        integer :: idim,jdim,kdim,ncon
        integer,pointer :: vsize
        !>
        !>
        integer,dimension(:),allocatable    :: adjncy
        integer,dimension(:),allocatable    :: options
        integer,dimension(:),allocatable    :: xadj
        integer,dimension(:),allocatable    :: adjwgt
        integer,dimension(:),allocatable    :: vwgt
        integer,dimension(:),allocatable    :: part
        !>
        !>
        !>
        real(kind=dprec),pointer    :: tpwgts,ubvec
        !>
        type(blocks_type),pointer   :: iblock
        type(overlap_type),pointer  :: mesh
        !>
        !>
        !>
        !>
        allocate(xadj(nobl+1))
        allocate(adjncy(num_bc))
        allocate(options(40))
        !>
        !>
        allocate(vwgt(nobl))
        !>
        allocate(adjwgt(num_bc))
        allocate(part(nobl))
        !>
        !>
        !>
        mesh    => grids(imesh)
        vsize   => null()
        tpwgts  => null()
        ubvec   => null()
        !>
        !>
        !>
        load_flags  = 1
        ncon        = 1
        nbl         = 1
        nparts      = nodes
        nvtxs       = nobl
        do nm=1,40
            options(nm) = 0
        end do
        !>
        !>
        if(load_flags .ne. 1)then
            !>
            !>
            !> used the circulation distribution method to mapping the blocks to processors
            !>
            !>
            if(mgflags .ne. 0)then
                if(nodes .ge.1)then
                    do nm=1,nobl
                        n2p(nbl) = mod(nm,nodes)+1
                        do nccg=1,global_level-1
                            nbl      = nbl + 1
                            n2p(nbl) = n2p(nbl-1)
                        end do
                        nbl = nbl + 1
                    end do
                else
                    do nm=1,nobl
                        n2p(nbl) = 0
                        do nccg=1,global_level-1
                            nbl      = nbl + 1
                            n2p(nbl) = n2p(nbl-1)
                        end do
                        nbl = nbl + 1
                    end do
                endif
            else
                if(nodes .ge. 1)then
                    do nm=1,nobl
                        n2p(nm) =  mod(nm,nodes)+1
                    end do
                else
                    do nm=1,nobl
                        n2p(nm) =  0
                    end do
                endif
            end if
        else
            !>call the METIS API to load balance
            !>
            !>
            if(nodes .eq. 0)then
                do nm=1,nblocks
                    n2p(nm) = 0
                end do
                return
            end if
            !>
            if(nodes .eq. 1)then
                do nm=1,nblocks
                    n2p(nm) = 1
                end do
                return
            end if
            !>
            !>
            !>
            !>
            nbl = 1
            if(nodes .eq. nobl)then
                do nm=1,nobl
                    if(mgflags .ne. 0)then
                        n2p(nbl) = nm
                        do nccg = 1,global_level-1
                            nbl = nbl + 1
                            n2p(nbl) = nm
                        end do
                        nbl = nbl + 1
                    else
                        n2p(nbl)    =   nm
                        nbl         =   nbl +1
                    end if
                end do
                !>
                return
            end if
            !>
            !>
            if(nodes .gt. nobl)then
                write(*,*) 'the processor larger than the blocks,check the processor for parallel !!!!'
                stop
            end if
            !> create the graph
            !>
            !>
            !> 1---2---3---4
            !> |   |   |   |
            !> |   |   |   |
            !> 5---6---7---8
            !> |   |   |   |
            !> |   |   |   |
            !> 9--10--11--12
            !>
            !> have 34 edges
            !>
            !> xadj =[1,3,6,9,11,14,18,22,25,27,30,33,35]
            !>
            !>
            !> have 34 elements for the contact points
            !> adjncy=[2,5,1,3,6,2,4,7,3,8,1,6,9,2,5,7,10,3,6,8,11,4,7,12,5,10,6,9,11,7,10,12,8,11]
            !>
            !>
            !>
            !> create the weights of the points
            nbl = 1
            do nm=1,nblocks,global_level
                iblock => mesh%blocks(nm)
                vwgt(nbl)   = (iblock%idim-1)*(iblock%jdim-1)*(iblock%kdim-1)
!                write(113,'("nbl=",i4," size=",i10)')nbl,vwgt(nbl)
                nbl         = nbl + 1
            end do
            !>
            !>
            ll = 0
            xadj(1) = 1
            nbl = 1
            do nm=1,nblocks,global_level
                do num_of_bc=1,num_bc
                    if(mesh%bcs(num_of_bc)%block_index .eq. nm .and. mesh%bcs(num_of_bc)%bc_type == 'cut')then
                        ll = ll + 1
                        adjncy(ll) = (mesh%bcs(num_of_bc)%iblock_index-1)/global_level + 1
!                        write(111,'("bc_num=",i8," ll=",i10," iblock=",i10)') num_of_bc,ll,mesh%bcs(num_of_bc)%iblock_index
                    end if
                end do
                nbl         = nbl + 1
                xadj(nbl)   = ll  + 1
!                write(114,'("nbl=",i4," ll=",i10)') nbl,xadj(nbl)
            end do
            !>
            !> create the weights of the edges
            !>
            lll = 1
            do nm=1,nblocks,global_level
                do num_of_bc = 1,num_bc
                    mm = (adjncy(lll)-1)*global_level + 1
                    if(mesh%bcs(num_of_bc)%block_index .eq. nm .and. mesh%bcs(num_of_bc)%iblock_index .eq. mm)then
                        dim1        = abs(mesh%bcs(num_of_bc)%istart - mesh%bcs(num_of_bc)%iend)+1
                        dim2        = abs(mesh%bcs(num_of_bc)%jstart - mesh%bcs(num_of_bc)%jend)+1
                        dim3        = abs(mesh%bcs(num_of_bc)%kstart - mesh%bcs(num_of_bc)%kend)+1
                        adjwgt(lll) = dim1*dim2*dim3
!                        write(112,'("block=",2i4," iblock=",i4," bc_num=",i6," face size=",3i5)') lll,nm,mm,num_of_bc,dim1,dim2,dim3
                        lll         = lll + 1
                    end if
                    !>
                    !>
                end do
            end do
            !>
            !>
            do nm=1,nobl
                part(nm)= 0
            end do
            !>
            !>
            !>begin call metis api to mapping the blocks to processors
            !>
            call METIS_SetDefaultOptions(options)
            !>
            options(18)=1
            !>
            !> used the METIS API to DD the graph
            !>
            call METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt,vsize,&
                                     adjwgt,nparts,tpwgts,ubvec,options,&
                                     objval,part)
            !>
            !>
            !>
            !> must set the coarse blocks and finer blocks at the same processor
            !> to make the communications minimize
            !>
!            do nm=1,nobl
!                write(1101,'("iblock=",i8," and the processors=",i8)') nm,part(nm)
!            end do
            nbl = 1
            if(mgflags .ne. 0)then
                do nm=1,nobl
                    n2p(nbl)= part(nm)
                    do nccg=1,global_level-1
                        nbl = nbl + 1
                        n2p(nbl) = part(nm)
                    end do
!                    write(1102,'("nbl=",i5)') nbl
                    nbl = nbl + 1
                end do
            else
                do nm=1,nobl
                    n2p(nm) = part(nm)
                end do
            end if
            !>
            !>
        end if
        !>
        !>
        !>
        deallocate(adjncy)
        deallocate(adjwgt)
        deallocate(vwgt)
        deallocate(xadj)
        deallocate(part)
        deallocate(options)
        !>
        !>
        return
        !>
    end subroutine load_balance
