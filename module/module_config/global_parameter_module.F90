!   _______________________________________________________________________________
!    _______/\\\\\\\\\_______/\\\\\\\\\______/\\\\\\\\\\\\\\__/\\\\\\\\\\\\\\\___________
!     _____/\/////////______/\\\///////______\/\\\//////////__\/\/////////////___________
!      ___/\\\/____________/\\\/______________\/\\\____________\/\\\________/\\\\________
!       __/\\\____________/\\\_________________\/\\\\\\\\\\\\\__\/\\\_______\////________
!        _\/\\\___________\/\\\_________________\/\\\/////////___\/\\\________/\\\\______
!         _\//\\\__________\//\\\________________\/\\\____________\/\\\_______\////______
!          __\///\\\_________\///\\\______________\/\\\____________\/\\\______/\\\\______
!           ____\///\\\\\\\\\___\////\\\\\\\\______\/\\\____________\/\\\\\\\\\////______
!            ______\/////////_______\////////_______\///_____________\/////////__________
!            __________________________________________________________________________________
!----------------------------------------------------------------------------------------------
!>  subroutine input_files_cgns
!>  last edit 2016-05-24
!>  last edit by liuxz
!----------------------------------------------------------------------------------------------
!********************************config类类型成员变量的声明和定义*****************************************
!**************************全局变量设置module，包括了来流参数，多重网格参数，限制器等重要参数*************
!*********************开发者可以自由的再次添加全局变量，当然也可以自定义module类型来添加需要的全局变量****

    module global_parameter
        !>!来流参数
        !>
        !>
        implicit none
        !the global parameter from the input files
        !and check by the programming
        !>
        !>
        !>mach number and the
        !>modify mach number if necessary
        real*8:: mach,xmach
        !>
        !> angle of attack
        real*8:: alpha
        !>
        !> angle of sideslip
        real*8:: beta
        !>
        !> the initial variable for input flow
        !>
        real*8:: r0,p0,h0,u0,v0,w0
        !>
        !> height and if height not equit -1
        !> must re-calculate
        real*8:: height
        !>
        !> reynolds number
        real*8:: reue
        !>
        !> temperature
        real*8:: tinf
        !>
        !>prandtl number for laminar flow
        real*8:: prandtl
        !>
        !>prandtl number for turbulent
        real*8:: prandtl_tur
        !>
        !> adiabatic wall twall=0
        !> have heat transfer wall twall = 1
        real*8:: twall
        !>
        !> rmre = xmach/reue
        real*8:: rmre
        !>
        !>reference value  of temperature
        real*8:: atm_t
        !>
        !>reference value  of pressures
        real*8:: atm_p
        !>
        !>reference value  of density
        real*8:: atm_r
        !>
        real*8 :: atm_mu
        !>reference value  of area
        real*8:: a_ref
        !>
        !>reference value  of length
        real*8:: c_ref
        !>
        !>
        real*8 :: b_ref
        !>
        !> moment point at x-directions
        real*8:: x_ref
        !>
        !> moment point at y-directions
        real*8:: y_ref
        !>
        !> moment point at z-directions
        real*8:: z_ref
        !> gamma constant number
        real*8:: gamma
        !>
        !>gm1 = gamma -1.
        real*8:: gm1
        !>
        !> 2.- gamma
        real*8:: gm2
        !>
        !> gm    = 1./(gamma -1.)
        real*8:: gm
        !>
        !>g0gm1 =  gamma / gm1
        real*8:: g0gm1
        !>
        !> cfl number
        real*8:: cfl
        !>
        !>constant
        real*8:: diim
        !>
        !> values for limiter
        !>一般>-1为二阶迎风格式,<-1一阶空间精度
        real*8:: xkap
        !>
        !>
        !>!熵修正参数，钝体参数范围
        !> 0.01  < epsa_r < 0.4
        real*8:: epsa_r
        !>
        !>icyc=1,the total error at blocks
        real*8:: uniterr
        !>
        !> dt < 0  steady computing
        !> dt > 0 unsteady computing
        real*8:: dtt
        !>
        !>
        real*8 :: delta_dt
        !>
        !>双时间步tau-ts迭代cfl数
        real*8:: taucfl
        !>
        !>!格式精度
        real*8:: spadif
        !>
        !>!implicit residual smoothing coefficient for residual
        real*8:: smooth_r_coe
        !>
        !>!implicit residual smoothing coefficient for solutions
        real*8:: smooth_c_coe
        !>
        !>
        real*8 :: cl
        !>
        !>
        real*8 :: cd
        !***********************************************************
        !***********************************************************
        !***********************************************************
        !***********************************************************
        !***********************************************************
        !***********************************************************
        !>
        integer dprec
#if defined single
        !> used  single precision,as the real*4
        parameter(dprec = 4)
#else
        !> used the double precision,as the real*8
        parameter(dprec = 8)
#endif
        !>
        !>the max coordinate at the dimension
        !>
        integer :: ijkmax
        !>
        !> 2d grid_dim = 2
        !> 3d grid_dim = 3
        integer :: grid_dim
        !>
        !> limiter=0,no limiter
        !> limiter=1,van ablaba limiter
        !> limiter=2,min-mod limiter
        !> limiter=3
        integer :: limiter
        !>
        !>!空间格式选择
        integer :: ifds
        !***********************************************************
        !***********************************************************
        !> the parameter for the full multi-grid methods
        !>
        !> the parameter for the full  multi-grid methods
        !> mesh-seq-flags = 1 used the FMG methods
        integer :: meshseqflags
        !>
        !> the level of the finer mesh to the coarser mesh
        !> when not used the multi-grid methods then mesh-sequence ==1
        !> parameter for the multi-grid methods
        integer :: meshseque
        !>
        !> a parameter for the level iterations of FMG methods
        integer :: imeshseque
        !***********************************************************
        !***********************************************************
        !> the parameter for the multi-grid methods
        !>
        !> the parameter flags for the multi-grid methods
        !>
        integer :: mgflags
        !>
        !> the parameter for the full multi-grid methods
        !> when not used the multi-grid methods the global_level =1
        !>
        integer :: global_level
        !>
        !> a parameter fot the mg methods
        integer :: ncg
        !>
        !>vcycle =1 the mg used the v-cycle else used w-cycle
        integer :: vcycle
        !>
        !> wcycle = 1 for one-w-cycle =2 for double-w-cycle
        integer :: wcycle
        !>
        !>for multi-grid and full multi-grid about v-cycle and w-cycle
        integer :: kode
        integer :: mode
        !***********************************************************
        !***********************************************************
        !>
        !> the parameter for the 1-1 cut boundary conditions
        !> when the nodes want to exchange the data
        !> the array size of the buffer for the exchange data
        integer :: i_size1to1
        !>
        !>  附加迭代标记
        integer :: additer
        !>
        !> the number of the blocks of the grid
        !>
        integer :: nblocks
        !>
        !>
        integer :: num_bc
        !>
        !>it used for implicit residual methods for residual
        integer :: irs
        !>
        !> it used for correction the solutions
        integer :: ics

        integer :: imp
        !>
        !>
        !> the parameter for the lu-sgs
        integer :: nvidi
        !>
        !>
        !> the flags for the viscous and inviscid
        !> nvisc == 0 : euler computing
        !> nvisc == 1 : n-s  computing laminar flow
        !> nvisc == 2 : baldwin_lomax turbulent model
        !> nvisc == 3 : sa  turbulent model
        !>
        integer :: nvisc
        !>
        !> the computing  model of the three directions
        !> like when the boundary conditions is no-slip wall
        !> so the ivisc_i >= 1
        !>
        integer ::ivisc_i,ivisc_j,ivisc_k
        !>
        !> the parameter for steady computing
        !> steady == 0 then the steady computing
        !> steady == 1 then the unsteady computing
        integer :: unstd
        !>
        !> 0:the symmetry plane is x-z plane with z up/
        !> 1:the symmetry plane is x-y plane with y up
        integer :: ialph
        !>
        !> !时间精度控制参数
        integer :: ita
        !>
        !> !续算开关
        integer :: restart
        !>
        !>!更新续算文件的间隔
        integer :: ntres
        !>
        !> the icyc of the restart computing
        !>续算的迭代开始步数
        integer :: icyc_restart
        !>
        !>!非定常迭代步数 for the unsteady
        integer :: ntstep
        !>
        !>
        !>
        integer :: icyc
        !>
        !> !多重网格迭代次数
        integer :: ncycle(5)
        !>
        !>
        !>
        integer :: total_ncyc
        !>
        !>
        integer :: ides
        !>
        !>
        integer :: nobl
        !>
        !> the total number of grid mesh
        !>
        integer :: num_grid
        !>
        !!动网格开关
        integer :: dyngrid
        !>
        !>!重叠网格标记
        integer :: overlap
        !>
        !> number of overlap meshes
        integer :: meshes
        !> number of overlap meshes
        integer :: imesh
        !>
        !> !刚体运动标记
        integer :: irigb
        !>
        !> !运动网格分组标记
        integer :: icat
        !>
        !> the grid file input
        !>
        character*80 grid_files
        !>
        !> the mesh type of the grid file
        !> right now 3.0 just can used the CGNS type grid to computing
        !> and we will edit the parallel CGNS type code for a large-scale case
        !>
        character*25 mesh_type
        !>
        !>
        !>
        contains

    end module global_parameter



