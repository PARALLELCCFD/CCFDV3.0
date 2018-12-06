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
    subroutine input_global
        !>
        use global_parameter
        use nodes_paras
        implicit none
        !>
        !>
#if defined PMPI
        include "mpif.h"
#endif

        logical :: alive
        integer :: i,ierr
        integer :: unit
        real(kind=dprec) :: tmp
#if defined PMPI
        !>
        integer ::send_buff_i
        integer ::send_buff_r
        parameter(send_buff_i=34)
        parameter(send_buff_r=21)
        !>
        integer,dimension(:),allocatable::ii_send_buffer
        real(kind=dprec),dimension(:),allocatable::ii_send_buffer_real
        !>
        !>
        !>
        allocate(ii_send_buffer(send_buff_i))
        allocate(ii_send_buffer_real(send_buff_r))
#endif
        !>
        !>
#if defined PMPI
        if(myid ==myhost)then
#endif
            !>   open the input file
            inquire(file='input.dat',exist=alive)
            if(.not.alive) then
                write(*,*)'error:input file not exist'
                stop
            endif
            open(31,file='input.dat',form='formatted',status='old')

            !>  read the name of grid file
            !>  read the head of the input file
            do i=1,9
                read(31,*)
            enddo
            !>
            !>*Grid File*(file name)
            read(31,*) grid_files
            !> open the output files
            !>
            open(unit=61,file='ccfd.out',form='formatted',status='unknown')
            open(unit=71,file='ccfd.resid',form='formatted',status='unknown')
            open(unit=63,file='ccfd.restart',form='formatted',status='unknown')


            !>  read the physical parameter
            !>  the coordinate  and axis directions
            read(31,*)
            !> mesh type
            !> grids type
            read(31,*) mesh_type
            read(31,*)
            !>
            !> =3:three dimension
            !> =2:two dimension
            read(31,*) grid_dim
            read(31,*)
            !>
            !> =0:steady;
            !> =1:unsteady;
            read(31,*) unstd
            read(31,*)
            !>
            !> =0: x-z plane with z-down;
            !> =1: x-y plane with y-up;
            read(31,*) ialph
#if defined DEBUG
            write(61,*) grid_files
            write(61,*) mesh_type
            write(61,'('' griddim='',i6)') grid_dim
            write(61,'('' unsteady='',i6)') unstd
            write(61,'('' ialph='',i6)') ialph
#endif
            !>
            !>physical conditions
            !>
            read(31,*)
            read(31,*)
            !>
            !>Mach number,non-dimensional, based on the free-stream values
            read(31,*) mach !mach number
            read(31,*)
            !>
            !> Angle of attack,degrees, only for compressible flows
            read(31,*) alpha
            !>
            alpha = alpha/57.295779513082
            !>
            !>
            read(31,*)
            !>
            !>Side slip angle,degrees, only for compressible flows
            read(31,*) beta
            !>
            !>
            beta = beta/57.295779513082
            !>
            !>
            read(31,*)
            !>
            !> Reynolds number,free-stream Reynolds number per unit grid length (millions) re/L)
            read(31,*) reue
            read(31,*)
            !>
            !>
            read(31,*) height
            read(31,*)
            !>
            !> Free-stream temperature (288.15 K by default)
            read(31,*) tinf
            !>
            !>
            if(tinf  .lt. 0.0)then
                write(*,'("warning::the Free-stream temperature  is negative")')
                write(*,'("warning::and default temperature is 288.15")')
                tinf = 288.15
            end if
            !>
            !>
#if defined DEBUG
            write(61,'('' mach='',f15.5)') mach
            write(61,'('' alpha='',f15.5)') alpha
            write(61,'('' beta='',f15.5)') beta
            write(61,'('' height='',f15.5)') height
            write(61,'('' reue='',f15.5)') reue
            write(61,'('' tinf='',f15.5)') tinf
#endif

            !>model of  geometry parameter
            read(31,*)
            read(31,*)
            !>
            !>Reference area for force coefficients
            read(31,*) a_ref !
            read(31,*)
            !>
            !>Reference length
            read(31,*) c_ref
            read(31,*)
            !>
            !> reference span used to compute non-dimensional moments
            read(31,*) b_ref !
            read(31,*)
            !>
            !>Moment center , Reference origin for moment computation
            read(31,*) x_ref,y_ref,z_ref !
            if(ialph ==1 )then
                tmp   = y_ref
                y_ref = -z_ref
                z_ref = tmp
            end if

#if defined DEBUG
            write(61,'('' sref='',f15.5)') a_ref
            write(61,'('' cref='',f15.5)') c_ref
            write(61,'('' bref='',f15.5)') b_ref
            write(61,'('' x_ref='',f15.5)') x_ref
            write(61,'('' y_ref='',f15.5)') y_ref
            write(61,'('' z_ref='',f15.5)') z_ref
#endif
            !>
            !>
            if (height .ne. -1) then
                call sta_aero(height,reue,mach,c_ref)
            end if
            !>
            !>
            if(reue .lt. 0.0)then
                write(*,'("warning::the reue parameter  is negative")')
                write(*,'("warning::and default reue is 1.0")')
                reue = 1.0
            end if
            reue = reue * 1.0e+06
            !>
            !>
            !> Time
            !>
            read(31,*)
            read(31,*)
            !>
            !> <0:steady computing；>0，dual-time schemes unsteady methods,
            !> this version just can used steady computing ,the values keep <0.0
            !>
            read(31,*) dtt
            read(31,*)
            !>
            !>CFL number
            read(31,*) cfl
            delta_dt= 0.0
            taucfl  = 7.500
            ntstep  = 1
            ita     = -1

#if defined DEBUG
            write(61,'('' dt='',f15.5)') dtt
            write(61,'('' cfl_num='',f15.5)') cfl
#endif

            !> Numerical parameter
            !>
            !>
            read(31,*)
            read(31,*)
            !> limiter
            !> =0:no limiter;=1: van albada;=2:min-mod limiter;=3: SPEKREIJSE VANKAT
            read(31,*) limiter
            read(31,*)
            !>
            !> value<-1.0:no limiter;
            !> when xkap≠1/3,used the min-mod,
            !> xkap=1/3 can used van albada or SPEKREIJSE VANKAT
            read(31,*) xkap
            read(31,*)
            !>
            !> correction the eigenvalue,0.01<epsa_r<0.4,when <0.0 no corrections
            !> when mach > 0.7 used the schemes
            !>
            read(31,*) epsa_r
            read(31,*)
            !>
            !>turbulent
            !> =0: euler;
            !> =1: n-s;
            !> =2: B-L model;
            !> =3: SA model;
            read(31,*) nvisc
            !>
            !>
            if(nvisc .eq. 2)then
                ivisc_i = 0
                ivisc_j = nvisc
                ivisc_k = 0
            else
                ivisc_i = nvisc
                ivisc_j = nvisc
                ivisc_k = nvisc
            end if
            !>
            !>
            ifds = 1
            spadif = 0.33333
            !>
            !>
#if defined DEBUG
            write(61,'('' limiter='',i6)') limiter
            write(61,'('' xkap='',e20.12)') xkap
            write(61,'('' epsa_r='',e20.12)') epsa_r
            write(61,'('' turb='',i6)') nvisc

#endif

            !>Restart
            read(31,*)
            read(31,*)
            !>
            !> restart flags
            !> =0:no restart;=1: output the restart files for restart computing
            read(31,*) restart
            read(31,*)
            !>
            !>iterations interval for each output the restart files
            !>
            read(31,*) ntres
            read(31,*)

#if defined DEBUG
            write(61,'('' irestart='',i6)') restart
            write(61,'('' ntrestart='',i6)') ntres
#endif

            !>Multi-grid methods
            read(31,*)
            !>
            !> Multi-grid Flag
            !> =0:not used the multi-gird methods; ≠0 used the multi-grid method
            read(31,*) mgflags
            read(31,*)
            !>
            !> multi-grid level
            !> level of multi-grids method ,when the value equal 1,mean not used the MG methods
            read(31,*) global_level
            read(31,*)
            !>
            !> VCYCLE
            !> =0:used the w-cycle schemes; =1:used the v-cycle
            read(31,*) vcycle
            read(31,*)
            !>
            !> WCYCLE
            !> =0:used the v-cycle schemes;=1;used the w-cycle
            read(31,*) wcycle
            read(31,*)
            !>
            !> implicit residual smoother
            !> =0:not used the IRS methods;=1:used the IRS methods
            read(31,*) irs
            read(31,*)
            !>
            !> SMOOTH_R_COE
            !> IRS schemes coefficient
            read(31,*) smooth_r_coe
            read(31,*)
            !>
            !>implicit correction smoother
            !>=0:not used the ICS methods;=1:used the ICS methods
            read(31,*) ics
            read(31,*)
            !>
            !>SMOOTH_C_COE
            !>ICS schemes coefficient
            read(31,*) smooth_c_coe
            read(31,*)
            !>
            !>
            additer = 0
            !>
            !> Mesh Sequencing methods
            !>
            read(31,*)
            !>
            !>Mesh Sequencing Flag
            !>=0:not used the mesh Sequencing methods;
            !>≠0:used the mesh Sequencing and must set the iterations for each level
            read(31,*) meshseqflags
            read(31,*)
            !>
            !>
            read(31,*) meshseque
            read(31,*)
            !>
            !> Mesh Sequencing Level
            !> =1:just used the multi-grid methods;
            !> >1:the level of mesh Sequencing methods
            if(meshseqflags .le. 0 )then
                meshseque  = 1
            end if
            !>
            !>
            if(meshseqflags .ne. 0)then
                do i =1,meshseque
                    read(31,*) ncycle(i)
                end do
                !>
                total_ncyc  = 0
                do i=1,meshseque
                    total_ncyc = total_ncyc + ncycle(i)
                end do
            else
                read(31,*) ncycle(1)  !
                total_ncyc = ncycle(1)
            endif

#if defined DEBUG
            write(61,'('' seque='',i6)') meshseque
            write(61,'('' seqflags='',i6)') meshseqflags
            write(61,'('' mgflags='',i6)') mgflags
            write(61,'('' irs='',i6)') irs
            write(61,'('' ics='',i6)') ics
            if(meshseqflags .ne. 0)then
                write(61,'('' ncycle='',i6)') (ncycle(i),i=1,meshseque)
            else
                write(61,'('' ncycle='',i6)') ncycle(1)
            end if
            write(61,'('' total iterations='',i6)') total_ncyc
            write(61,'('' global_level='',i6)') global_level
            write(61,'('' vcycle='',i6)') vcycle
            write(61,'('' wcycle='',i6)') wcycle
#endif

            !>
            !>
            dyngrid = 0
            irigb   = 0
            icat    = 7
            overlap = 0
            if(overlap .eq. 0 )then
                meshes =  1
            else
                read(31,*)
                !>
                !>overlap meshes number
                read(31,*) meshes
            end if
#if defined DEBUG
            write(*,*) 'the input files read end-files'
#endif
            !> when not used the mg method guarantee the global_level=1
            !> the global_levle will be used at the solver
            if(.not. mgflags) global_level = 1
            if(.not. meshseqflags) meshseque=1
            !>
            read(31,*)
            close(31)
#if defined PMPI
            !>
            !>
            !**********************************************
            !bcast the global parameter to other nodes from host node
            !**********************************************

            ii_send_buffer(1)  = grid_dim
            ii_send_buffer(2)  = unstd
            ii_send_buffer(3)  = ialph
            ii_send_buffer(4)  = ntstep
            ii_send_buffer(5)  = ita
            ii_send_buffer(6)  = ifds
            ii_send_buffer(7)  = limiter
            ii_send_buffer(8)  = nvisc
            ii_send_buffer(9)  = restart
            ii_send_buffer(10) = ntres
            ii_send_buffer(11) = meshseque
            ii_send_buffer(12) = meshseqflags
            ii_send_buffer(13) = global_level
            ii_send_buffer(14) = additer
            !ii_send_buffer(15) = ncycle
            ii_send_buffer(16) = mgflags
            ii_send_buffer(17) = dyngrid
            ii_send_buffer(18) = irigb
            ii_send_buffer(19) = icat
            ii_send_buffer(20) = overlap
            ii_send_buffer(21) = irs
            ii_send_buffer(22) = ics
            ii_send_buffer(23) = vcycle
            ii_send_buffer(24) = wcycle
            if(meshseqflags .ne. 0)then
                ii_send_buffer(25)= ncycle(1)
                ii_send_buffer(26)= ncycle(2)
                ii_send_buffer(27)= ncycle(3)
                ii_send_buffer(28)= ncycle(4)
                ii_send_buffer(29)= ncycle(5)
            else
                ii_send_buffer(25)= ncycle(1)
            end if
            ii_send_buffer(30) =  ivisc_i
            ii_send_buffer(31) =  ivisc_j
            ii_send_buffer(32) =  ivisc_k
            ii_send_buffer(33) =  total_ncyc
            ii_send_buffer(34) =  meshes
            !>
            !>
            !>
            ii_send_buffer_real(1)  = mach
            ii_send_buffer_real(2)  = alpha
            ii_send_buffer_real(3)  = beta
            ii_send_buffer_real(4)  = height
            ii_send_buffer_real(5)  = reue
            ii_send_buffer_real(6)  = tinf
            ii_send_buffer_real(7)  = a_ref
            ii_send_buffer_real(8)  = c_ref
            ii_send_buffer_real(9)  = b_ref
            ii_send_buffer_real(10) = x_ref
            ii_send_buffer_real(11) = y_ref
            ii_send_buffer_real(12) = z_ref
            ii_send_buffer_real(13) = dtt
            ii_send_buffer_real(14) = delta_dt
            ii_send_buffer_real(15) = cfl
            ii_send_buffer_real(16) = taucfl
            ii_send_buffer_real(17) = spadif
            ii_send_buffer_real(18) = smooth_r_coe
            ii_send_buffer_real(19) = smooth_c_coe
            ii_send_buffer_real(20) = xkap
            ii_send_buffer_real(21) = epsa_r
            !>
            !>
            !>
            call mpi_bcast(ii_send_buffer,34,mpi_integer,myhost,mycomm,ierr)
            call mpi_bcast(ii_send_buffer_real,21,mpi_double_precision,myhost,mycomm,ierr)
            !>
            !>
            !>
            !>
        !> the myid == myhost
        else
            !**********************************************
            !no-myhost recv the global parameter
            !**********************************************
            !>
            !>
            !>
            call mpi_bcast(ii_send_buffer,34,mpi_integer,myhost,mycomm,ierr)
            call mpi_bcast(ii_send_buffer_real,21,mpi_double_precision,myhost,mycomm,ierr)
            grid_dim     = ii_send_buffer(1)
            unstd        = ii_send_buffer(2)
            ialph        = ii_send_buffer(3)
            ntstep       = ii_send_buffer(4)
            ita          = ii_send_buffer(5)
            ifds         = ii_send_buffer(6)
            limiter      = ii_send_buffer(7)
            nvisc        = ii_send_buffer(8)
            restart      = ii_send_buffer(9)
            ntres        = ii_send_buffer(10)
            meshseque    = ii_send_buffer(11)
            meshseqflags = ii_send_buffer(12)
            global_level = ii_send_buffer(13)
            additer      = ii_send_buffer(14)
!           ncycle       = ii_send_buffer(15)
            mgflags      = ii_send_buffer(16)
            dyngrid      = ii_send_buffer(17)
            irigb        = ii_send_buffer(18)
            icat         = ii_send_buffer(19)
            overlap      = ii_send_buffer(20)
            irs          = ii_send_buffer(21)
            ics          = ii_send_buffer(22)
            vcycle       = ii_send_buffer(23)
            wcycle       = ii_send_buffer(24)
            if(meshseqflags .ne. 0)then
                ncycle(1)       = ii_send_buffer(25)
                ncycle(2)       = ii_send_buffer(26)
                ncycle(3)       = ii_send_buffer(27)
                ncycle(4)       = ii_send_buffer(28)
                ncycle(5)       = ii_send_buffer(29)
            else
                ncycle(1)       = ii_send_buffer(25)
            end if
            ivisc_i          = ii_send_buffer(30)
            ivisc_j          = ii_send_buffer(31)
            ivisc_k          = ii_send_buffer(32)
            total_ncyc       = ii_send_buffer(33)
            meshes           = ii_send_buffer(34)
            !>
            !>
            !>
            mach         = ii_send_buffer_real(1)
            alpha        = ii_send_buffer_real(2)
            beta         = ii_send_buffer_real(3)
            height        = ii_send_buffer_real(4)
            reue         = ii_send_buffer_real(5)
            tinf         = ii_send_buffer_real(6)
            a_ref        = ii_send_buffer_real(7)
            c_ref        = ii_send_buffer_real(8)
            b_ref         = ii_send_buffer_real(9)
            x_ref        = ii_send_buffer_real(10)
            y_ref        = ii_send_buffer_real(11)
            z_ref        = ii_send_buffer_real(12)
            dtt          = ii_send_buffer_real(13)
            delta_dt     = ii_send_buffer_real(14)
            cfl          = ii_send_buffer_real(15)
            taucfl       = ii_send_buffer_real(16)
            spadif       = ii_send_buffer_real(17)
            smooth_r_coe = ii_send_buffer_real(18)
            smooth_c_coe = ii_send_buffer_real(19)
            xkap         = ii_send_buffer_real(20)
            epsa_r       = ii_send_buffer_real(21)

        endif
#endif
        !
        gamma = 1.4
        gm1   = gamma -1.
        gm2   = 2.- gamma
        gm    = 1./(gamma -1.)
        g0gm1 = gamma / gm1
        !> prandtl for laminar flow
        prandtl  = 0.7200
        !> prandtl for turbulent
        prandtl_tur   = 0.900
        !>
        !>
        atm_t = 288.15
        atm_r = 1.225
        atm_p = 101325
        !>
        !>
        nvidi = 3
        diim  = 1.05
        !>

        atm_mu = 1.7161e-5*(atm_t/273.16)**1.5*(273.16+124.0)/(atm_t+124.0)
        !>
        !>
        !>
#if defined PMPI
        !>
        deallocate(ii_send_buffer)
        deallocate(ii_send_buffer_real)
#endif
        !>
    end subroutine  input_global
