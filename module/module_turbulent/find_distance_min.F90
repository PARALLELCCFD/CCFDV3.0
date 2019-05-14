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
!>  subroutine find_distance_min
!>  last edit 2018-04-18
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
    subroutine find_distance_min()
        !>
        !> minimum distance function for field equation turbulence models
        use global_parameter
        use mesh_overlap_module
        use blocks_module
        use metric_module
        use nodes_var
        use bc_module
        use nodes_paras

		implicit none
#if defined PMPI
        include "mpif.h"
#endif
        type(overlap_type),pointer :: mesh
        type(blocks_type),pointer  :: iblock
        type(blocks_type),pointer  :: iblock_c
        type(bc_types),pointer     :: ibc
        !>
        !>
        !>
        integer :: i,j,k
        integer :: ii,jj,kk,iii,jjj,kkk,inum,inbl,nbl_c
        integer :: iis,iie,jjs,jje,kks,kke
        integer :: idim,jdim,kdim
        integer :: idim_c,jdim_c,kdim_c
        integer :: nbl,num_of_bc,myidnbl
        integer :: is,ie,js,je,ks,ke
        integer :: iend,jend,kend
        integer :: istart,jstart,kstart
        integer :: irange,iindex
        integer :: npts,wallnpts,ipoint,ierr,iipoint,iiipoint,iwall,bb_index,minbox,ilevel,iilevel
        integer :: iright,ileft,ir,il,istack
        integer :: isplit,blockindex(2,1024)
        integer :: chose_point,numtriangles
        !>
        !>
        !>
        character(len=24) :: face
        !>
        !>
        !>
        real(kind =dprec) :: factr,blank1,blank2,blank3
        real(kind =dprec) :: hole,large,sort_temp
        real(kind =dprec) :: px,py,pz,xs,ys,zs
        real(kind =dprec) :: coorminx,coorminy,coorminz,coormaxx,coormaxy,coormaxz
        real(kind =dprec) :: coormin(3),coormax(3)
        real(kind =dprec) :: lcoormin(3),lcoormax(3),rcoormin(3),rcoormax(3)
        real(kind =dprec) :: xrange,yrange,zrange
        real(kind =dprec) :: stack(7,200)
        real(kind =dprec) :: blockrange(6,1024),distance_tmp(1024),distance_min,test_min,smin
        !>
        !>
        integer           ,dimension(:),allocatable :: ipoint_array
        integer           ,dimension(:),allocatable :: ipoint_array_tmp
        integer           ,dimension(:),allocatable :: isplit_array
        integer           ,dimension(:),allocatable :: neighbor_index
        integer           ,dimension(:),allocatable :: iwall_array
        integer           ,dimension(:),allocatable :: bounding_box
        integer           ,dimension(:,:),allocatable :: neighbor_array
        integer           ,dimension(:,:),allocatable :: itemp
        integer           ,dimension(:,:,:),allocatable :: iblock_chose_point
        !>
        !>
        !>
        real(kind = dprec),dimension(:),    allocatable :: send_array
        real(kind = dprec),dimension(:),    allocatable :: surface1
        real(kind = dprec),dimension(:),    allocatable :: surface2
        real(kind = dprec),dimension(:),    allocatable :: surface3
        real(kind = dprec),dimension(:),    allocatable :: surface4
        real(kind = dprec),dimension(:),    allocatable :: test_tmp
        real(kind = dprec),dimension(:,:),  allocatable :: temp

        real(kind = dprec),dimension(:,:,:),allocatable :: iblock_dist
        real(kind = dprec),dimension(:,:,:),allocatable :: iblock_dist_c
        !>
        !>
        hole  = 1.0
        large = 1.e30
        mesh => grids(imesh)
        iwall = 0
        wallnpts = 0
        iiipoint = 0
        coorminx =  1.0e20
        coorminy =  1.0e20
        coorminz =  1.0e20
        coormaxx = -1.0e20
        coormaxy = -1.0e20
        coormaxz = -1.0e20
        !>
        !>
        !>
        if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
            !>
            do num_of_bc=1,num_bc,global_level
                ibc => mesh%bcs(num_of_bc)
                if(ibc%bc_type == 'wallviscid' )then
                    iwall = iwall + 1
                end if
            end do
            ibc => null()
            !>
            !>
            allocate(iwall_array(iwall))
            !>
            !>
            iwall = 0
            do num_of_bc=1,num_bc,global_level
                ibc => mesh%bcs(num_of_bc)
                if(ibc%bc_type == 'wallviscid' )then
                    iwall = iwall + 1
                    iwall_array(iwall) = num_of_bc
                    is = ibc%istart
                    ie = ibc%iend
                    js = ibc%jstart
                    je = ibc%jend
                    ks = ibc%kstart
                    ke = ibc%kend
                    if(is == ie )then
                        wallnpts   = wallnpts + (abs(je-js)+2)*(abs(ke-ks)+2)
                    end if
                    if(js == je )then
                        wallnpts   = wallnpts + (abs(ie-is)+2)*(abs(ke-ks)+2)
                    end if
                    if(ks == ke )then
                        wallnpts   = wallnpts + (abs(je-js)+2)*(abs(ie-is)+2)
                    end if
                end if
            end do
            !>
            !>
            ibc => null()
            !>
            minbox = sqrt(real(wallnpts))
            minbox = max(minbox,50)
            ilevel = 10
            !>
            allocate(surface1(wallnpts))
            allocate(surface2(wallnpts))
            allocate(surface3(wallnpts))
            allocate(surface4(wallnpts))
            allocate(bounding_box(wallnpts))
            allocate(temp(wallnpts,4))
            allocate(itemp(wallnpts,9))
            allocate(neighbor_array(wallnpts,8))
            allocate(neighbor_index(wallnpts))
            allocate(ipoint_array(wallnpts))
            allocate(ipoint_array_tmp(wallnpts))
            allocate(isplit_array(wallnpts))
            allocate(test_tmp(wallnpts))
            !>
            !>
            do i=1,wallnpts
                neighbor_index(i)   = 0
                neighbor_array(i,1) = 0
                neighbor_array(i,2) = 0
                neighbor_array(i,3) = 0
                neighbor_array(i,4) = 0
                neighbor_array(i,5) = 0
                neighbor_array(i,6) = 0
                neighbor_array(i,7) = 0
                neighbor_array(i,8) = 0
            end do
            !>
            !>
            !>
            !>
            iiipoint  = 0
            !>
            do nbl=1 ,nblocks,global_level
                npts     = 0
                ipoint   = 0
                !>
                !>
                iblock => mesh%blocks(nbl)
                !>
                !>
                idim  =   iblock%idim
                jdim  =   iblock%jdim
                kdim  =   iblock%kdim
                !>
                !>
                !> get the solid surface
                !>
                do num_of_bc=1,iwall
                    ibc => mesh%bcs(iwall_array(num_of_bc))
                    if(ibc%block_index == nbl .and. ibc%bc_type == 'wallviscid')then
                        !>
                        !>
                        is = ibc%istart
                        ie = ibc%iend
                        js = ibc%jstart
                        je = ibc%jend
                        ks = ibc%kstart
                        ke = ibc%kend
                        !>
                        !>
                        if(is == ie .and. is .eq. 2)then
                            npts  = npts + 4*(abs(je-js)+2)*(abs(ke-ks)+2)
                        end if
                        if(is == ie .and. ie .gt. 2)then
                            npts  = npts + 4*(abs(je-js)+2)*(abs(ke-ks)+2)
                        end if
                        if(js == je .and. js .eq. 2)then
                            npts  = npts + 4*(abs(ie-is)+2)*(abs(ke-ks)+2)
                        end if
                        if(js == je .and. je .gt. 2)then
                            npts  = npts + 4*(abs(ie-is)+2)*(abs(ke-ks)+2)
                        end if
                        if(ks == ke .and. ks .eq. 2)then
                            npts  = npts + 4*(abs(je-js)+2)*(abs(ie-is)+2)
                        end if
                        if(ks == ke .and. ke .gt. 2)then
                            npts  = npts + 4*(abs(je-js)+2)*(abs(ie-is)+2)
                        end if
                    end if
                end do
                !>
                !>
                if(npts == 0 ) cycle
                allocate(send_array(npts))
!                if(myid .ne. 0 ) write(*,'("the send buffer sizer=",i6,i4)') npts,nbl
                !>
                !>
#if defined PMPI
                myidnbl = n2p(nbl)
                if(myid == n2p(nbl))then

#endif
                    do num_of_bc=1,iwall
                        ibc => mesh%bcs(iwall_array(num_of_bc))
                        if(ibc%block_index == nbl .and. ibc%bc_type == 'wallviscid')then
                            !>
                            !>
                            is = ibc%istart
                            ie = ibc%iend
                            js = ibc%jstart
                            je = ibc%jend
                            ks = ibc%kstart
                            ke = ibc%kend
                            !>
                            !>
                            if(is .gt. ie)then
                                istart = is + 1
                                iend   = ie
                            elseif(is .lt. ie)then
                                istart = is
                                iend   = ie + 1
                            else
                                istart = is
                                iend   = ie
                            end if
                            !>
                            if(js .gt. je)then
                                jstart = js + 1
                                jend   = je
                            elseif(js .lt. je)then
                                jstart = js
                                jend   = je + 1
                            else
                                jstart = js
                                jend   = je
                            end if
                            !>
                            if(ks .gt. ke)then
                                kstart = ks + 1
                                kend   = ke
                            elseif(ks .lt. ke)then
                                kstart = ks
                                kend   = ke + 1
                            else
                                kstart = ks
                                kend   = ke
                            end if
                            !>
                            !>
                            if(is == ie .and. is .eq. 2)then
                                face = 'istart'
                            end if
                            if(is == ie .and. ie .gt. 2)then
                                face = 'iend'
                            end if
                            if(js == je .and. js .eq. 2)then
                                face = 'jstart'
                            end if
                            if(js == je .and. je .gt. 2)then
                                face = 'jend'
                            end if
                            if(ks == ke .and. ks .eq. 2)then
                                face = 'kstart'
                            end if
                            if(ks == ke .and. ke .gt. 2)then
                                face = 'kend'
                            end if
                            !>
                            !>
                            !>
                            if(face == 'istart')then
                                !>
                                !>
                                do kk=kstart,kend
                                    kks = min(kk,kend)
                                    kke = max(kk,kstart)
                                    do jj=jstart,jend
                                        jjs = min(jj,jend)
                                        jje = max(jj,jstart)
                                        blank1    = max(iblock%cell(is,jjs,kks)%blank,iblock%cell(is,jje,kks)%blank)
                                        blank2    = max(iblock%cell(is,jjs,kke)%blank,iblock%cell(is,jje,kke)%blank)
                                        blank3    = max(blank1,blank2)
                                        factr     = (1.0-blank3)*hole*large
                                        !>
                                        !>
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(is,jj,kk)%x
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(is,jj,kk)%y
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(is,jj,kk)%z
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  factr
                                    end do
                                end do
                                !>
                                !>
                            else if(face == 'iend')then
                                !>
                                !>
                                do kk=kstart,kend
                                    kks = min(kk,kend)
                                    kke = max(kk,kstart)
                                    do jj=jstart,jend
                                        jjs = min(jj,jend)
                                        jje = max(jj,jstart)
                                        blank1    = max(iblock%cell(ie,jjs,kks)%blank,iblock%cell(ie,jje,kks)%blank)
                                        blank2    = max(iblock%cell(ie,jjs,kke)%blank,iblock%cell(ie,jje,kke)%blank)
                                        blank3    = max(blank1,blank2)
                                        factr     = (1.0-blank3)*hole*large
                                        !>
                                        !>
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ie+1,jj,kk)%x
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ie+1,jj,kk)%y
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ie+1,jj,kk)%z
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  factr
                                    end do
                                end do
                                !>
                                !>
                            else if(face == 'jstart')then
                                !>
                                !>
                                do kk=kstart,kend
                                    kks = min(kk,kend)
                                    kke = max(kk,kstart)
                                    do ii=istart,iend
                                        iis = min(ii,iend)
                                        iie = max(ii,istart)
                                        blank1    = max(iblock%cell(iis,js,kks)%blank,iblock%cell(iie,js,kks)%blank)
                                        blank2    = max(iblock%cell(iis,js,kke)%blank,iblock%cell(iie,js,kke)%blank)
                                        blank3    = max(blank1,blank2)
                                        factr     = (1.0-blank3)*hole*large
                                        !>
                                        !>
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,js,kk)%x
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,js,kk)%y
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,js,kk)%z
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  factr
                                    end do
                                end do
                                !>
                                !>
                            else if(face == 'jend')then
                                !>
                                !>
                                do kk=kstart,kend
                                    kks = min(kk,kend)
                                    kke = max(kk,kstart)
                                    do ii=istart,iend
                                        iis = min(ii,iend)
                                        iie = max(ii,istart)
                                        blank1    = max(iblock%cell(iis,je,kks)%blank,iblock%cell(iie,je,kks)%blank)
                                        blank2    = max(iblock%cell(iis,je,kke)%blank,iblock%cell(iie,je,kke)%blank)
                                        blank3    = max(blank1,blank2)
                                        factr     = (1.0-blank3)*hole*large
                                        !>
                                        !>
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,je+1,kk)%x
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,je+1,kk)%y
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,je+1,kk)%z
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  factr
!                                        if(nbl==18) write(92802,'(i8,4e24.16)') ipoint/4,iblock%coordinate(ii,je+1,kk)%x,iblock%coordinate(ii,je+1,kk)%y,iblock%coordinate(ii,je+1,kk)%z,factr
                                    end do
                                end do
                                !>
                                !>
                            else if(face == 'kstart')then
                                !>
                                !>
                                do jj=jstart,jend
                                    jjs = min(jj,jend)
                                    jje = max(jj,jstart)
                                    do ii=istart,iend
                                        iis = min(ii,iend)
                                        iie = max(ii,istart)
                                        blank1    = max(iblock%cell(iis,jjs,ks)%blank,iblock%cell(iie,jjs,ks)%blank)
                                        blank2    = max(iblock%cell(iis,jjs,ks)%blank,iblock%cell(iie,jje,ks)%blank)
                                        blank3    = max(blank1,blank2)
                                        factr     = (1.0-blank3)*hole*large
                                        !>
                                        !>
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,jj,ks)%x
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,jj,ks)%y
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,jj,ks)%z
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  factr
                                    end do
                                end do
                                !>
                                !>
                            else if(face == 'kend')then
                                !>
                                !>
                                do jj=jstart,jend
                                    jjs = min(jj,je)
                                    jje = max(jj,js)
                                    do ii=istart,iend
                                        iis = min(ii,iend)
                                        iie = max(ii,istart)
                                        blank1    = max(iblock%cell(iis,jjs,ke)%blank,iblock%cell(iie,jjs,ke)%blank)
                                        blank2    = max(iblock%cell(iis,jjs,ke)%blank,iblock%cell(iie,jje,ke)%blank)
                                        blank3    = max(blank1,blank2)
                                        factr     = (1.0-blank3)*hole*large
                                        !>
                                        !>
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,jj,ke+1)%x
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,jj,ke+1)%y
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  iblock%coordinate(ii,jj,ke+1)%z
                                        ipoint             =  ipoint + 1
                                        send_array(ipoint) =  factr
                                    end do
                                end do
                                !>
                                !>
                            else
                                write(*,*) 'have some error at the boundary conditions '
                                write(*,*) 'at the subroutine find_distance_min.F90'
                                stop
                            end if
                            !>
                            !>
                        end if
                    end do  !num_of_bc=1,iwall
!                    write(*,'("send the sizer of data=",i6,2i4)') ipoint,myid,nbl
                    !>
                    !>
#if defined PMPI
                end if
                !> send the buffer to the others processors
                !> receive the data on non-local nodes
                !>
                !>
                call mpi_bcast(send_array,npts,mpi_double_precision,myidnbl,mycomm,ierr)
#endif
                !> collect the surface solid data on non-local nodes
                !>
                !>
                iipoint  = 1
                do num_of_bc = 1,iwall
                    ibc => mesh%bcs(iwall_array(num_of_bc))
                    if(ibc%block_index == nbl .and. ibc%bc_type == 'wallviscid')then
                        !>
                        !>
                        is = ibc%istart
                        ie = ibc%iend
                        js = ibc%jstart
                        je = ibc%jend
                        ks = ibc%kstart
                        ke = ibc%kend
                        !>
                        !>
                        if(is .gt. ie)then
                            istart = is + 1
                            iend   = ie
                        elseif(is .lt. ie)then
                            istart = is
                            iend   = ie + 1
                        else
                            istart = is
                            iend   = ie
                        end if
                        !>
                        if(js .gt. je)then
                            jstart = js + 1
                            jend   = je
                        elseif(js .lt. je)then
                            jstart = js
                            jend   = je + 1
                        else
                            jstart = js
                            jend   = je
                        end if
                        !>
                        if(ks .gt. ke)then
                            kstart = ks + 1
                            kend   = ke
                        elseif(ks .lt. ke)then
                            kstart = ks
                            kend   = ke + 1
                        else
                            kstart = ks
                            kend   = ke
                        end if
                        !>
                        !>
                        if(is == ie )then
                            face = 'i'
                        end if
                        if(js == je)then
                            face = 'j'
                        end if
                        if(ks == ke )then
                            face = 'k'
                        end if
                        !>
                        !>
                        !>each interior point has 8 neighbors,defined as
                        !>
                        !> p1 --- p2 --- p3
                        !> |    /  |  \   |
                        !> |   /   |   \  |
                        !> |  /    |    \ |
                        !> p4 --- ip --- p5
                        !> |  \    |    / |
                        !> |   \   |   /  |
                        !> |    \  |  /   |
                        !> p6 --- p7 --- p8
                        !>
                        !>
                        if(face == 'i')then
                            irange   = abs(jstart-jend)+1
                            do kk = kstart,kend
                                do jj= jstart,jend
                                    iindex                  = 0
                                    iiipoint                = iiipoint + 1
                                    ipoint_array(iiipoint)  = iiipoint
                                    surface1(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    surface2(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    surface3(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    surface4(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    do kkk = max(kk-1,kstart),min(kk+1,kend)
                                        do jjj =max(jj-1,jstart),min(jj+1,jend)
                                            if(kkk .ne. kk .and. jjj .ne. jj)then
                                                iindex  = iindex + 1
                                                neighbor_array(iiipoint,2*iindex-1) =  iiipoint+(jjj-jj)+1
                                                neighbor_array(iiipoint,2*iindex)   =  iiipoint+(kkk-kk)*irange+1
                                            end if
                                        end do
                                    end do
                                    neighbor_index(iiipoint) = iindex
                                    !>
                                    !>
                                end do
                            end do
                        end if
                        if(face == 'j')then
                            irange   = abs(istart-iend)+1
                            do kk = kstart,kend
                                do ii= istart,iend
                                    iindex                  = 0
                                    iiipoint                = iiipoint + 1
                                    ipoint_array(iiipoint)  = iiipoint
                                    surface1(iiipoint)      = send_array(iipoint)
!                                    if(myid .ne. 0) write(180709,'("surface1=",e24.16,i8)') surface1(iiipoint),ipoint_array(iiipoint)
                                    iipoint                 = iipoint  + 1
                                    surface2(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    surface3(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    surface4(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    do kkk = max(kk-1,kstart),min(kk+1,kend)
                                        do iii =max(ii-1,istart),min(ii+1,iend)
                                            if(kkk .ne. kk .and. iii .ne. ii)then
                                                iindex  = iindex + 1
                                                neighbor_array(iiipoint,2*iindex-1) =  iiipoint+(iii-ii) + 1
                                                neighbor_array(iiipoint,2*iindex)   =  iiipoint+(kkk-kk)*irange + 1
!                                                if(myid .ne. 0) write(*,'("neighbors points=",5i4)') iiipoint,2*iindex-1,2*iindex,neighbor_array(iiipoint,2*iindex-1),neighbor_array(iiipoint,2*iindex)
                                            end if
                                        end do
                                    end do
                                    neighbor_index(iiipoint) = iindex
                                    !>
                                    !>
                                end do
                            end do
                        end if
                        if(face == 'k')then
                            irange  = abs(istart-iend)+1
                            do jj = jstart,jend
                                do ii= istart,iend
                                    iindex                  = 0
                                    iiipoint                = iiipoint + 1
                                    ipoint_array(iiipoint)  = iiipoint
                                    surface1(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    surface2(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    surface3(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    surface4(iiipoint)      = send_array(iipoint)
                                    iipoint                 = iipoint  + 1
                                    do jjj = max(jj-1,jstart),min(jj+1,jend)
                                        do iii =max(ii-1,istart),min(ii+1,iend)
                                            if(jjj .ne. jj .and. iii .ne. ii)then
                                                iindex  = iindex + 1
                                                neighbor_array(iiipoint,2*iindex-1) =  iiipoint+(iii-ii) + 1
                                                neighbor_array(iiipoint,2*iindex)   =  iiipoint+(jjj-jj)*irange + 1
                                            end if
                                        end do
                                    end do
                                    neighbor_index(iiipoint) = iindex
                                    !>
                                    !>
                                end do
                            end do
                        end if
                        !>
                        !>
                        !>
                        !>
                    end if
                end do !num_of_bc=1,iwall
                !>
                !>check the data and buffer size
!                if(iipoint .gt. 1 .and. (iipoint-1) .ne. npts)then
                if((iipoint-1) .ne. npts)then
                    write(*,*) '******************************************************************************'
                    write(*,*) 'the buffer of surface point error occurred  at find_distance_min subroutine!!!'
                    write(*,'(" receive the data is ",i10," but the buffer size=",i10," at block=",i8)') iipoint-1,npts,nbl
                    write(*,*) '******************************************************************************'
                    stop
                end if
                !>
                !>
                deallocate(send_array)
                !>
                !>
            end do !nbl=1,nblocks
            !>
            !>
            !>Heap sort methods for surface points with respect to x coordinate
            !>
            !>
!            do i=1,wallnpts
!               write(92801,'("heap=",i8,4e24.16)')i,surface1(i),surface2(i),surface3(i),surface4(i)
!            end do
            if(wallnpts .gt. 1)then
                call heap_sort(wallnpts,surface1,ipoint_array)
            end if
!            do i=1,wallnpts
!               write(9280,'("heap=",i8,e24.16)') ipoint_array(i),surface1(ipoint_array(i))
!            end do
            !>
            !>
            !>
            do i=1,wallnpts
!                write(*,'("heap_sort=",e24.16,i8)') surface1(ipoint_array(i)),ipoint_array(i)
                temp(i,1)  = surface1(ipoint_array(i))
                temp(i,2)  = surface2(ipoint_array(i))
                temp(i,3)  = surface3(ipoint_array(i))
                temp(i,4)  = surface4(ipoint_array(i))
                itemp(i,1) = neighbor_index(ipoint_array(i))
                itemp(i,2) = neighbor_array(ipoint_array(i),1)
                itemp(i,3) = neighbor_array(ipoint_array(i),2)
                itemp(i,4) = neighbor_array(ipoint_array(i),3)
                itemp(i,5) = neighbor_array(ipoint_array(i),4)
                itemp(i,6) = neighbor_array(ipoint_array(i),5)
                itemp(i,7) = neighbor_array(ipoint_array(i),6)
                itemp(i,8) = neighbor_array(ipoint_array(i),7)
                itemp(i,9) = neighbor_array(ipoint_array(i),8)
            end do
            !>
            !>
            do i=1,wallnpts
                surface1(i) = temp(i,1)
                surface2(i) = temp(i,2)
                surface3(i) = temp(i,3)
                surface4(i) = temp(i,4)
                neighbor_index(i)   =  itemp(i,1)
                neighbor_array(i,1) =  itemp(i,2)
                neighbor_array(i,2) =  itemp(i,3)
                neighbor_array(i,3) =  itemp(i,4)
                neighbor_array(i,4) =  itemp(i,5)
                neighbor_array(i,5) =  itemp(i,6)
                neighbor_array(i,6) =  itemp(i,7)
                neighbor_array(i,7) =  itemp(i,8)
                neighbor_array(i,8) =  itemp(i,9)
            end do
            !>
            !>
            !>
            do i=1,wallnpts
                ii                   = ipoint_array(i)
                ipoint_array_tmp(ii) = i
            end do
            !>
            !>
            do i=1,wallnpts
                numtriangles = neighbor_index(i)
                do jj=1,2*numtriangles
                    ii                   = neighbor_array(i,jj)
                    neighbor_array(i,jj) = ipoint_array_tmp(ii-1)
                end do
            end do
!            do i=1,wallnpts
!                write(1807201,'("neighbor=",10i6)')neighbor_index(i),i,neighbor_array(i,1),neighbor_array(i,2),neighbor_array(i,3),neighbor_array(i,4),&
!                neighbor_array(i,5),neighbor_array(i,6),neighbor_array(i,7),neighbor_array(i,8)
!            end do
            !>
            !>
            !> create a bounding box contain all the surface points
            !>
            !>
            !>
            !>
            coormin(1) = coorminx
            coormin(2) = coorminy
            coormin(3) = coorminz
            coormax(1) = coormaxx
            coormax(2) = coormaxy
            coormax(3) = coormaxz
            do i=1,wallnpts
                !>
                if(surface1(i) .lt. coormin(1)) coormin(1) = surface1(i)
                if(surface2(i) .lt. coormin(2)) coormin(2) = surface2(i)
                if(surface3(i) .lt. coormin(3)) coormin(3) = surface3(i)
                !>
                !>
                if(surface1(i) .gt. coormax(1)) coormax(1) = surface1(i)
                if(surface2(i) .gt. coormax(2)) coormax(2) = surface2(i)
                if(surface3(i) .gt. coormax(3)) coormax(3) = surface3(i)
            end do
            !>
            !>
            iilevel  = 1
            istack   = 1
            isplit   = 0
            inum     = 0

            !>
            !>search the coordinate to find all points that fall within a bounding box
            !>
            bb_index = 0
            do i=1,wallnpts
                if( surface1(i) .ge. coormin(1) .and. &
                    surface1(i) .le. coormax(1) .and. &
                    surface2(i) .ge. coormin(2) .and. &
                    surface2(i) .le. coormax(2) .and. &
                    surface3(i) .ge. coormin(3) .and. &
                    surface3(i) .le. coormax(3)   )then
                    bb_index               = bb_index + 1
                    bounding_box(bb_index) = i
                end if
            end do
            stack(1,istack) = coormin(1)
            stack(2,istack) = coormin(2)
            stack(3,istack) = coormin(3)
            stack(4,istack) = coormax(1)
            stack(5,istack) = coormax(2)
            stack(6,istack) = coormax(3)
            stack(7,istack) = iilevel
            !>
            !>
            !>
            !> subdivide the bounding box
            !>
            do while(istack .ge. 1)
                !>
                bb_index = 0
                coormin(1) =  stack(1,istack)
                coormin(2) =  stack(2,istack)
                coormin(3) =  stack(3,istack)
                coormax(1) =  stack(4,istack)
                coormax(2) =  stack(5,istack)
                coormax(3) =  stack(6,istack)
                !>
                !>
                iilevel    =  stack(7,istack)
                istack     =  istack - 1
                !>search the coordinate to find all points that fall within a bounding box
                !>
                do i=1,wallnpts
                    if( surface1(i) .ge. coormin(1) .and. &
                        surface1(i) .le. coormax(1) .and. &
                        surface2(i) .ge. coormin(2) .and. &
                        surface2(i) .le. coormax(2) .and. &
                        surface3(i) .ge. coormin(3) .and. &
                        surface3(i) .le. coormax(3)   )then
                        bb_index               = bb_index + 1
                        bounding_box(bb_index) = i
!                        if(myid == 0) write(1807201,'(2i6)') bb_index,bounding_box(bb_index)
                    end if
                end do
!                if(myid == 0) write(*,'("in stack=",i8,2i4)') bb_index,iilevel,istack
                !>
                !>
                !>
                if(bb_index  .ge. minbox .and. iilevel .lt. ilevel)then
                    !>
                    !>
                    !>
                    iilevel = iilevel + 1
                    !>
                    !>
                    xrange = coormax(1) - coormin(1)
                    yrange = coormax(2) - coormin(2)
                    zrange = coormax(3) - coormin(3)
                    !>
                    if(xrange .gt. yrange .and. xrange .gt. zrange)then
                        !>
                        !>
                        call shell_sort(bb_index,surface1,bounding_box)
                        !>
                        iright = bb_index/2
                        do i = 1,bb_index/2
                            il     = bounding_box(iright)
                            ir     = bounding_box(iright+1)
                            if(surface1(il) .ne. surface1(ir) .or. iright .ge. bb_index) exit
                            iright = iright + 1
                            if(iright .ge. bb_index) exit
                        end do
!                        if(myid == 0) write(*,'("iright=",i6)') iright
                        !>
                        !>
                    else if(yrange .gt. xrange .and.  yrange .gt. zrange)then
                        !>
                        !>
                        call shell_sort(bb_index,surface2,bounding_box)
                        !>
                        iright = bb_index/2
                        do i = 1,bb_index/2
                            il     = bounding_box(iright)
                            ir     = bounding_box(iright+1)
                            if(surface2(il) .ne. surface2(ir) .or. iright .ge. bb_index) exit
                            iright = iright + 1
                            if(iright .ge. bb_index) exit
                        end do
!                        if(myid == 0) write(*,'("jright=",i6)') iright
                        !>
                        !>
                    else
                        !>
                        !>
                        call shell_sort(bb_index,surface3,bounding_box)
                        !>
                        iright = bb_index/2
                        do i = 1,bb_index/2
                            il     = bounding_box(iright)
                            ir     = bounding_box(iright+1)
                            if(surface3(il) .ne. surface3(ir) .or. iright .ge. bb_index) exit
                            iright = iright + 1
                            if(iright .ge. bb_index) exit
                        end do
!                        if(myid == 0) write(*,'("kright=",i6)') iright
                        !>
                        !>
                    end if
                    !>
                    !>
                    !> subdivide the bounding box
                    lcoormin(1) = coorminx
                    lcoormin(2) = coorminy
                    lcoormin(3) = coorminz
                    lcoormax(1) = coormaxx
                    lcoormax(2) = coormaxy
                    lcoormax(3) = coormaxz
                    !>
                    do i=1,iright
                        !>
                        if(surface1(bounding_box(i)) .lt. lcoormin(1)) lcoormin(1) = surface1(bounding_box(i))
                        if(surface2(bounding_box(i)) .lt. lcoormin(2)) lcoormin(2) = surface2(bounding_box(i))
                        if(surface3(bounding_box(i)) .lt. lcoormin(3)) lcoormin(3) = surface3(bounding_box(i))
                        !>
                        !>
                        if(surface1(bounding_box(i)) .gt. lcoormax(1)) lcoormax(1) = surface1(bounding_box(i))
                        if(surface2(bounding_box(i)) .gt. lcoormax(2)) lcoormax(2) = surface2(bounding_box(i))
                        if(surface3(bounding_box(i)) .gt. lcoormax(3)) lcoormax(3) = surface3(bounding_box(i))
                    end do

                    !>
                    !>
                    rcoormin(1) = coorminx
                    rcoormin(2) = coorminy
                    rcoormin(3) = coorminz
                    rcoormax(1) = coormaxx
                    rcoormax(2) = coormaxy
                    rcoormax(3) = coormaxz
                    do i=iright+1,bb_index
                        !>
                        if(surface1(bounding_box(i)) .lt. rcoormin(1)) rcoormin(1) = surface1(bounding_box(i))
                        if(surface2(bounding_box(i)) .lt. rcoormin(2)) rcoormin(2) = surface2(bounding_box(i))
                        if(surface3(bounding_box(i)) .lt. rcoormin(3)) rcoormin(3) = surface3(bounding_box(i))
                        !>
                        !>
                        if(surface1(bounding_box(i)) .gt. rcoormax(1)) rcoormax(1) = surface1(bounding_box(i))
                        if(surface2(bounding_box(i)) .gt. rcoormax(2)) rcoormax(2) = surface2(bounding_box(i))
                        if(surface3(bounding_box(i)) .gt. rcoormax(3)) rcoormax(3) = surface3(bounding_box(i))
                    end do
                    !>
                    !> push the coordinate max and min in the stack
                    !>
                    istack          = istack + 1
                    stack(1,istack) = rcoormin(1)
                    stack(2,istack) = rcoormin(2)
                    stack(3,istack) = rcoormin(3)
                    stack(4,istack) = rcoormax(1)
                    stack(5,istack) = rcoormax(2)
                    stack(6,istack) = rcoormax(3)
                    stack(7,istack) = iilevel
                   ! if(myid == 0) write(*,'("istack=",i4," left_indx",i6," left_bbmaxmin=",6e24.16)') istack,iright,lcoormin(1),lcoormin(2),lcoormin(3),lcoormax(1),lcoormax(2),lcoormax(3)
                    !>
                    !>
                    istack          = istack + 1
                    stack(1,istack) = lcoormin(1)
                    stack(2,istack) = lcoormin(2)
                    stack(3,istack) = lcoormin(3)
                    stack(4,istack) = lcoormax(1)
                    stack(5,istack) = lcoormax(2)
                    stack(6,istack) = lcoormax(3)
                    stack(7,istack) = iilevel
                    !if(myid == 0) write(*,'("istack=",i4," right_indx",i6," right_bbmaxmin=",6e24.16)') istack,bb_index-iright,rcoormin(1),rcoormin(2),rcoormin(3),rcoormax(1),rcoormax(2),rcoormax(3)
                    !>
                else
                    !>
                    !>
!                    if(myid == 0) write(*,'("  box stack=",i8,2i4)') bb_index,iilevel,istack
                    if(bb_index .gt. 0 )then
                        isplit   = isplit  + 1
                        blockrange(1,isplit) =  coormin(1)
                        blockrange(2,isplit) =  coormin(2)
                        blockrange(3,isplit) =  coormin(3)
                        blockrange(4,isplit) =  coormax(1)
                        blockrange(5,isplit) =  coormax(2)
                        blockrange(6,isplit) =  coormax(3)
                        !>
                        blockindex(1,isplit) = bb_index
                        blockindex(2,isplit) = inum + 1
                       ! if(myid == 0) write(*,'("  blockindex=",4i10)') minbox,wallnpts,bb_index,blockindex(2,isplit)
                        !>
                        do jjj=1,bb_index
                            inum = inum + 1
                            isplit_array(inum) = bounding_box(jjj)
                        end do
                        !>
                        !>
                    end if
                    !>
                    !>
                end if
                !>
                !>
            end do !> do while(istack .ge. 1 )
            !>
            !>
            !> calculate distance from each point in the field to points on
            !> the viscous surfaces and then look at neighboring triangles
            !> and calculate the finer meshes
            !>
            iblock => null()
            do nbl = 1,nblocks,global_level
#if defined PMPI
                if(myid .eq. n2p(nbl))then
#endif
!                    write(*,'("calculate the distance @myid=",i4," and #iblock=",i5)')myid,nbl
                    iblock => mesh%blocks(nbl)
                    idim  = iblock%idim
                    jdim  = iblock%jdim
                    kdim  = iblock%kdim
                    !>
                    !>
                    allocate(iblock_dist(2:idim+1,2:jdim+1,2:kdim+1))
                    allocate(iblock_chose_point(2:idim+1,2:jdim+1,2:kdim+1))
                    !>
                    !>
!                    write(*,'("axis=",4i6)') idim,jdim,kdim,global_level
                    !>
                    do k=2,kdim+1
                        do j=2,jdim+1
                            do i=2,idim+1
                                !>
                                do inum = 1,isplit
                                    px = iblock%coordinate(i,j,k)%x
                                    if(real(px) .le. blockrange(1,inum)) px = blockrange(1,inum)
                                    if(real(px) .ge. blockrange(4,inum)) px = blockrange(4,inum)
                                    !>
                                    py = iblock%coordinate(i,j,k)%y
                                    if(real(py) .le. blockrange(2,inum)) py = blockrange(2,inum)
                                    if(real(py) .ge. blockrange(5,inum)) py = blockrange(5,inum)
                                    !>
                                    pz = iblock%coordinate(i,j,k)%z
                                    if(real(pz) .le. blockrange(3,inum)) pz = blockrange(3,inum)
                                    if(real(pz) .ge. blockrange(6,inum)) pz = blockrange(6,inum)
                                    !>
!                                    write(*,'("inum=",4i4,3e24.16)') inum,i,j,k,px,py,pz
                                    distance_tmp(inum) = (iblock%coordinate(i,j,k)%x - px)**2 + &
                                                         (iblock%coordinate(i,j,k)%y - py)**2 + &
                                                         (iblock%coordinate(i,j,k)%z - pz)**2
!                                    write(*,'("bb=",i4,7e24.16)') inum,iblock%coordinate(i,j,k)%x,iblock%coordinate(i,j,k)%y,&
!                                    iblock%coordinate(i,j,k)%z,px,py,pz,distance_tmp(inum)
                                end do
                                !>
                                !> find nearest bounding box
                                !>
                                smin        = 1.0e34
                                chose_point = 0
                                do inum=1,isplit
                                    distance_min = distance_tmp(1)
                                    do jjj=2,isplit
                                        distance_min = min(distance_min,distance_tmp(jjj))
                                    end do
                                    do kk=1,isplit
                                        if(distance_min .eq. distance_tmp(kk)) exit
                                    end do
                                    !>
                                    distance_tmp(kk)  = 2.0e34
                                    !>
                                    if(distance_min .le. smin)then
                                        !>
                                        !>
                                        jjj = blockindex(2,kk)
                                        test_min = 1.0e34
                                        do kkk=1,blockindex(1,kk)
                                            xs = (surface1(isplit_array(kkk+jjj-1)))
                                            ys = (surface2(isplit_array(kkk+jjj-1)))
                                            zs = (surface3(isplit_array(kkk+jjj-1)))
                                            test_tmp(kkk) = (iblock%coordinate(i,j,k)%x-xs)**2 + &
                                                            (iblock%coordinate(i,j,k)%y-ys)**2 + &
                                                            (iblock%coordinate(i,j,k)%z-zs)**2 + &
                                                             surface4(isplit_array(kkk+jjj-1))
                                            !>
                                            !>
                                            test_min = min(test_tmp(kkk),test_min)
                                            !write(*,'("kkk=",3i6)') kk,kkk,kkk+jjj-1!,xs,ys,zs!,iblock%coordinate(i,j,k)%x
!                                            iblock%coordinate(i,j,k)%y,iblock%coordinate(i,j,k)%z,surface4(isplit_array(kkk+jjj-1))
                                        end do
!                                        write(*,'("in 100=",e24.16)') test_min
                                        !>
                                        !>
                                        if(test_min .lt. smin)then
                                            !>
                                            do kkk=1,blockindex(1,kk)
                                                if(test_min .eq. test_tmp(kkk)) exit
                                            end do
                                            smin = test_min
                                            chose_point = isplit_array(kkk+jjj-1)
                                        end if
                                        !>
                                        !>
                                    else
                                        exit
                                    end if
!                                    write(*,'("in 200=",i4,e24.16,i6)') inum,smin,chose_point
                                    !>
                                    !>
                                    !>
                                end do
                                !>
                                !>
                                iblock_dist(i,j,k)        = sqrt(smin)
                                iblock_chose_point(i,j,k) = chose_point
!                                write(1807212,'("iblock_dist=",e24.16,i6)') iblock_dist(i,j,k),iblock_chose_point(i,j,k)
                            end do
                        end do
                    end do
                    !>
                    !>
                    !> calculation to triangles
                    !>
                    do k=2,kdim+1
                        do j=2,jdim+1
                            do i=2,idim+1
                                chose_point  = iblock_chose_point(i,j,k)
                                numtriangles = neighbor_index(chose_point)
!                                write(1807213,'("chose=",2i6)') chose_point,numtriangles
!                                write(1807208,'("coordinate=",3e24.16)') iblock%coordinate(i,j,k)%x,iblock%coordinate(i,j,k)%y,iblock%coordinate(i,j,k)%z
                                if(numtriangles .eq. 4)then
                                    do inum =1,numtriangles
                                        il = neighbor_array(chose_point,2*inum-1)
                                        ir = neighbor_array(chose_point,2*inum)
!                                        write(1807213,'("chose=",3i6)') inum,il,ir
!                                        write(1807208,'("triangle1=",12e24.16)') surface1(chose_point),surface2(chose_point),surface3(chose_point),surface4(chose_point),&
!                                                       surface1(il),surface2(il),surface3(il),surface4(il),&
!                                                       surface1(ir),surface2(ir),surface3(ir),surface4(ir)
                                        call triangles(iblock_dist(i,j,k),iblock%coordinate(i,j,k)%x,iblock%coordinate(i,j,k)%y,iblock%coordinate(i,j,k)%z,&
                                                       surface1(chose_point),surface2(chose_point),surface3(chose_point),surface4(chose_point),&
                                                       surface1(il),surface2(il),surface3(il),surface4(il),&
                                                       surface1(ir),surface2(ir),surface3(ir),surface4(ir))
!                                        write(18072018,'("triangle_sol1=",3i4,e24.16)') i-1,j-1,k-1,iblock_dist(i,j,k)
                                    end do
                                else
                                    !>
                                    do inum = chose_point,1,-1
                                        if(surface1(inum) .lt. surface1(chose_point) ) exit
                                    end do
                                    if(inum .eq. 0) inum = 0
                                    !>
                                    !>
                                    inum = inum + 1
                                    !>
                                    do jjj= chose_point,wallnpts
                                        if(surface1(jjj) .gt. surface1(chose_point)) exit
                                    end do
                                    if(jjj .eq. wallnpts+1) jjj = wallnpts + 1
                                    jjj = jjj  - 1
                                    !>
                                    !>
!                                    write(1807214,'("chose1=",3i6)') inum,jjj,chose_point
                                    do kkk= inum,jjj
                                        if(surface2(kkk) .eq. surface2(chose_point) .and.  &
                                           surface3(kkk) .eq. surface3(chose_point) .and.  &
                                           neighbor_index(kkk) .ne. 4)then
                                           !>
                                           !>
                                           numtriangles = neighbor_index(kkk)
                                           do iii = 1,numtriangles
                                                il = neighbor_array(kkk,2*iii-1)
                                                ir = neighbor_array(kkk,2*iii)
!                                                write(1807214,'("chose2=",3i6)') iii,il,ir
                                                !>
                                                !>
!                                                write(1807208,'("triangle2=",12e24.16)') surface1(kkk),surface2(kkk),surface3(kkk),surface4(kkk),&
!                                                               surface1(il) ,surface2(il) ,surface3(il) ,surface4(il),&
!                                                               surface1(ir) ,surface2(ir) ,surface3(ir) ,surface4(ir)
                                                call triangles(iblock_dist(i,j,k),iblock%coordinate(i,j,k)%x,iblock%coordinate(i,j,k)%y,iblock%coordinate(i,j,k)%z,&
                                                               surface1(kkk),surface2(kkk),surface3(kkk),surface4(kkk),&
                                                               surface1(il) ,surface2(il) ,surface3(il) ,surface4(il),&
                                                               surface1(ir) ,surface2(ir) ,surface3(ir) ,surface4(ir))
!                                                write(18072018,'("triangle_sol2=",3i4,e24.16)') i-1,j-1,k-1,iblock_dist(i,j,k)
                                           end do
                                           !>
                                           !>
                                        end if
                                    end do
                                end if
                                !>
                                !>
                            end do
                        end do
                    end do
                    !>
                    !>
                    kkk = 1
                    do k=2,kdim
                        kkk = kkk +1
                        jjj = 1
                        do j=2,jdim
                            jjj = jjj +1
                            iii = 1
                            do i=2,idim
                                iii = iii + 1
                                iblock%turbulent(iii,jjj,kkk)%distance  = 0.125*(iblock_dist(i  ,j  ,k  ) + iblock_dist(i+1,j  ,k  ) + &
                                                                                 iblock_dist(i  ,j+1,k  ) + iblock_dist(i  ,j  ,k+1) + &
                                                                                 iblock_dist(i+1,j+1,k  ) + iblock_dist(i+1,j  ,k+1) + &
                                                                                 iblock_dist(i  ,j+1,k+1) + iblock_dist(i+1,j+1,k+1) )
!                              write(nbl+9000,'("axis=",3i4,e24.16)')iii,jjj,kkk,iblock%turbulent(iii,jjj,kkk)%distance
                            end do
                        end do
                    end do
!                    iii = 1
!                    do i=2,idim
!                        iii = iii +1
!                        kkk = 1
!                        do k=2,kdim
!                            kkk = kkk +1
!                            jjj = 1
!                            do j=2,jdim
!                                jjj = jjj + 1
!                                iblock%turbulent(iii,jjj,kkk)%distance  = 0.125*(iblock_dist(i  ,j  ,k  ) + iblock_dist(i+1,j  ,k  ) + &
!                                                                                 iblock_dist(i  ,j+1,k  ) + iblock_dist(i  ,j  ,k+1) + &
!                                                                                 iblock_dist(i+1,j+1,k  ) + iblock_dist(i+1,j  ,k+1) + &
!                                                                                 iblock_dist(i  ,j+1,k+1) + iblock_dist(i+1,j+1,k+1) )
!                                write(1807231,'("distance=",3i4,e24.16)') iii-1,jjj-1,kkk-1,iblock%turbulent(iii,jjj,kkk)%distance
!                            end do
!                        end do
!                    end do

                    !>
                    !> coarser level
                    !> to collocate the finer meshes mini-distance to a coarser meshes
                    !> and the coarser mesh and finer mesh at the same processor so
                    !> can used the finer mesh mini distance not used send/recv the data
                    if(mgflags .ne. 0 .and. global_level .gt. 1)then
                        nbl_c = nbl
                        do inbl=nbl_c+1,global_level+nbl-1
                            !>
                            !>
                            iblock   => mesh%blocks(inbl-1)
                            iblock_c => mesh%blocks(inbl  )
                            !>
                            !>
                            idim     = iblock%idim
                            jdim     = iblock%jdim
                            kdim     = iblock%kdim
                            idim_c   = iblock_c%idim
                            jdim_c   = iblock_c%jdim
                            kdim_c   = iblock_c%kdim
                            !>
                            !>
                            allocate(iblock_dist_c(2:idim_c+1,2:jdim_c+1,2:kdim_c+1))
                            !>
                            !>
                            kkk = 1
                            !>
                            if(idim .gt. 2)then
                                !>
                                !> for three dimension coordinate
                                !>
                                do k=2,kdim+1,2
                                    kkk = kkk + 1
                                    jjj = 1
                                    do j=2,jdim+1,2
                                        jjj = jjj + 1
                                        iii = 1
                                        do i=2,idim+1,2
                                            iii = iii + 1
                                            iblock_dist_c(iii,jjj,kkk) = iblock_dist(i,j,k)
                                        end do
                                    end do
                                end do
                                !>
                                !>
                            else
                                !>
                                !> for two dimension coordinate
                                !>
                                do k=2,kdim+1,2
                                    kkk = kkk + 1
                                    jjj = 1
                                    do j=2,jdim+1,2
                                        jjj = jjj + 1
                                        iii = 1
                                        do i=2,idim+1,1
                                            iii = iii + 1
                                            iblock_dist_c(iii,jjj,kkk) = iblock_dist(i,j,k)
                                        end do
                                    end do
                                end do
                                !>
                                !>
                            end if
                            !>
                            !>
                            !>
                            !>
                            kkk = 1
                            do k=2,kdim_c
                                kkk = kkk +1
                                jjj = 1
                                do j=2,jdim_c
                                    jjj = jjj +1
                                    iii = 1
                                    do i=2,idim_c
                                        iii = iii + 1
                                        iblock_c%turbulent(iii,jjj,kkk)%distance  = 0.125*(iblock_dist_c(i  ,j  ,k  ) + iblock_dist_c(i+1,j  ,k  ) + &
                                                                                           iblock_dist_c(i  ,j+1,k  ) + iblock_dist_c(i  ,j  ,k+1) + &
                                                                                           iblock_dist_c(i+1,j+1,k  ) + iblock_dist_c(i+1,j  ,k+1) + &
                                                                                           iblock_dist_c(i  ,j+1,k+1) + iblock_dist_c(i+1,j+1,k+1) )
                                    end do
                                end do
                            end do
                            !>
                            !>
                            do k=2,kdim_c+1
                                do j=2,jdim_c+1
                                    do i=2,idim_c+1
                                        iblock_dist(i,j,k)= iblock_dist_c(i,j,k)
                                    end do
                                end do
                            end do
                            !>
                            deallocate(iblock_dist_c)
                            !>
                            !>
                        end do ! inbl=nbl_c+1,global_level+nbl-1
                        !>
                        !>
                    end if !if(mgflags .ne. 0 .and. global_level .gt. 1)then
                    !>
                    !>
                    deallocate(iblock_dist)
                    deallocate(iblock_chose_point)
#if defined PMPI
                end if ! if(myid .eq. n2p(nbl))then
#endif
            end do !do nbl = 1,nblocks,global_level
            !>
            !>
            !>
            !>
        end if !if(ivisc_i .ge. 3 .or. ivisc_j .ge. 3 .or. ivisc_k .ge. 3)then
        !>
        !>
        deallocate(ipoint_array)
        deallocate(ipoint_array_tmp)
        deallocate(isplit_array)
        deallocate(neighbor_index)
        deallocate(surface1)
        deallocate(surface2)
        deallocate(surface3)
        deallocate(surface4)
        deallocate(test_tmp)
        deallocate(bounding_box)
        deallocate(temp)
        deallocate(itemp)
        deallocate(neighbor_array)
        !>
        !>
        return

	end subroutine
