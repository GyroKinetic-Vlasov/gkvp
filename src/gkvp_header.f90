MODULE GKV_header
!-------------------------------------------------------------------------------
!
!    Header for general use in the fluxtube code
!
!    Notes
!    -----
!      There are some restrictions on numerical parameters:
!         mod( global_nz, nprocz ) = 0, due to z parallelization.
!         mod( global_nv, nprocv ) = 0, due to v parallelization.
!         mod( global_nm+1, nprocm ) = 0, due to m parallelization.
!         nm>=3, due to fft and colli.
!         nzb>=2, due to 4th-oreder z derivative.
!      There are some recommendations on numerical parameters:
!         mod( nxw, nprocw ) = 0, due to w parallelization.
!         mod( global_ny+1, nprocw ) = 0, due to w parallelization.
!         nzb<=nz due to data copy in zfilter.
!
!    Update history
!    --------------
!      gkvp_f0.61 (S. Maeyama, Mar 2021)
!        - equib_type is extended from len=8 to len=15.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

  implicit none

  public

  integer, parameter :: DP = selected_real_kind(14)

!--------------------------------------
!  Dimension size (grid numbers)
!--------------------------------------
!  Global simulation domain 
!  in  x, y,z,v,m (0:2*nxw-1,  0:2*nyw-1,-global_nz:global_nz-1,1:2*global_nv,0:global_nm)
!  in kx,ky,z,v,m (   -nx:nx,0:global_ny,-global_nz:global_nz-1,1:2*global_nv,0:global_nm)

  integer, parameter :: nxw = 196, nyw = 20
  integer, parameter :: nx = 128, global_ny = 1 ! 2/3 de-aliasing rule
  integer, parameter :: global_nz = 12, global_nv = 24, global_nm = 7

  integer, parameter :: nzb = 2, &  ! the number of ghost grids in z
                        nvb = 2     ! the number of ghost grids in v and m

!--------------------------------------
!  Data distribution for MPI
!--------------------------------------

  integer, parameter :: nprocw = 1, nprocz = 2, nprocv = 4, nprocm = 2, nprocs = 1

!--------------------------------------
!  Parameters for variable sizes
!--------------------------------------
!  Local simulation domain 
!  in kx,ky,z,v,m (divided in ky,z,v,m) (   -nx:nx,      0:ny,-nz:nz-1,1:2*nv,0:nm)
!  in  x,ky,z,v,m (divided in ky,z,v,m) (0:2*nxw-1,      0:ny,-nz:nz-1,1:2*nv,0:nm)
!  in  y, x,z,v,m (divided in  x,z,v,m) (    0:nyw,0:nxw_size,-nz:nz-1,1:2*nv,0:nm)

  integer, parameter :: nxw_size = (2*nxw-1)/nprocw     ! local allocation size (0:nxw_size)
  integer, parameter :: ny       = global_ny / nprocw   ! local allocation size (0:ny)

  integer, parameter :: nz = global_nz / nprocz,          &
                        nv = global_nv / nprocv,          &
                        nm = (global_nm + 1) / nprocm - 1,&
                        ns = nprocs

  integer, parameter :: nxyz = (2*nx+1)*(ny+1)*(2*nz), &
                        nxy  = (2*nx+1)*(ny+1)

  integer, parameter :: nnx = nxw*2, nny = nyw*2


!--------------------------------------
!  Constants
!--------------------------------------

  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, &
                                 twopi = pi * 2._DP,         &
                                 eps = 0.0000000001_DP
  complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )


!--------------------------------------
!  Index Range
!--------------------------------------

! ---- y dimension -------
  integer :: ist_y                           ! local start index of y
  integer :: iend_y                          ! local end   index of y
  integer :: nsize_y                         ! local size of y
  integer :: ist1_y                          ! local start index of y for global start index 1 

  integer :: ist_y_g                         ! global start index of y
  integer :: iend_y_g                        ! global end   index of y

! ---- xw dimension  -------
  integer :: ist_xw                           ! local start index of xw
  integer :: iend_xw                          ! local end   index of xw
  integer :: nsize_xw                         ! local size of xw

  integer :: ist_xw_g                         ! global start index of xw
  integer :: iend_xw_g                        ! global end   index of xw


!--------------------------------------
!  Parameters for time
!--------------------------------------

  real(kind=DP) :: e_limit                           ! elapsed time limit of a job
  real(kind=DP) :: tend                              ! end time
  real(kind=DP) :: dtout_fxv, dtout_ptn, dtout_eng   ! time-spacing for output
  real(kind=DP) :: dtout_dtc                         ! time-spacing for dt control


!--------------------------------------
!  Configuration parameters to be 
!    initialized in init subroutine
!--------------------------------------

!  real(kind=DP), &
!    dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)  :: kvd, kvs
  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm)  :: vdx, vdy, vsy
  real(kind=DP), &
    dimension(-nx:nx,0:ny,-nz:nz-1,0:nm)         :: j0, j1, j2
  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: g0, ksq
  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: fct_poisson, fct_e_energy
  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: fct_ampere, fct_m_energy
  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm) :: fmx
  real(kind=DP), dimension(-nx:nx)               :: fctgt

!!!%%% Parameters for colli_full %%%
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm) :: xxa
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm,0:ns-1,0:ns-1) :: nu_h, nu_g, nu_d, nu_p
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm)               :: nu_hs, nu_gs, nu_ds, nu_ps
!!!  real(kind=DP), dimension(1:6,-nz:nz-1,1:2*nv,0:nm,0:ns-1,0:ns-1) :: x_tst, y_fld 
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm,0:ns-1,0:ns-1,1:2) :: c_t0
!!!!  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm,0:ns-1,1:6) :: vfunc
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm,0:ns-1,1:6) :: vfunc
!!!  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:6) :: jfunc
!!!  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: adbtc
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real(kind=DP), dimension(0:ns-1,0:ns-1) :: ctauiv, calpha, ctheta, cgamma, ceta, cxi

  real(kind=DP), dimension(-nx:nx)          :: kx
  real(kind=DP), dimension(0:ny)            :: ky
  real(kind=DP), dimension(-nz:nz-1)        :: zz, omg
  real(kind=DP), dimension(1:2*nv)          :: vl
  real(kind=DP), dimension(0:nm)            :: mu
  real(kind=DP), dimension(-nz:nz-1,0:nm)   :: vp, mir
  real(kind=DP), dimension(-nz:nz-1)        :: dvp
  real(kind=DP), dimension(-nz:nz-1)        :: dpara, rootg
  real(kind=DP), dimension(1:3,1:3,-nz:nz-1) :: gg_g !%%% For shearflow_lagrange
  real(kind=DP) :: s_hat_g, eps_r_g                  !%%% For shearflow_lagrange->running

  complex(kind=DP), dimension(0:ny)         :: ck
  integer, dimension(0:ny)                  :: dj

  real(kind=DP) :: dt_max, dt
  logical :: adapt_dt
                                                     !!! Parameters for the R0 units
  real(kind=DP), dimension(0:ns-1) ::   R0_Ln,  &    ! R0/Lns
                                        R0_Lt,  &    ! R0/Lts
                                           nu,  &    ! collision freq.   
                                         Anum,  &    ! mass number
                                         Znum,  &    ! charge number     
                                          fcs,  &    ! charge-density fraction 
                                          sgn,  &    ! signs of charge   
                                          tau,  &    ! T-ratio
                                         dns1        ! initial perturbation amp.
  real(kind=DP) :: dv, cfsrf, lambda_i, q_0, q_bar, beta, tau_ad, vmax
  real(kind=DP) :: mach, uprime, gamma_e, kxmin_g, kymin_g, tlim_exb
  real(kind=DP) :: Nref, Lref, Tref, Zeff
  integer       :: iFLR, icheck, ibprime, nx0
  real(kind=DP) :: baxfactor

  real(kind=DP) :: courant_num 

!--------------------------------------
!  Type of calculation
!--------------------------------------

  character(9)  :: calc_type, &  ! "linear", "lin_freq", "nonlinear"
                   z_bound,   &  ! "zerofixed", "outflow", "mixed"
                   z_filt,    &  ! "on", "off"
                   z_calc,    &  ! "cf4", "up5"
                   col_type,  &  ! "LB", "full", "lorentz"
                   time_advnc    ! "rkg4", "imp_colli", "auto_init"
  real(kind=DP) :: art_diff

  integer :: num_triad_diag


!--------------------------------------
!  Parameters for numerical settings
!--------------------------------------

  integer :: inum
  logical :: ch_res, init_random
  character(512) :: f_log, f_hst, f_phi, f_fxv, f_cnt
  character(15)  :: equib_type  ! "analytic", "s-alpha", "s-alpha-shift",
                                ! "circ-MHD", "vmec", "eqdsk", "slab"

! --- unit numbers for I/O
  integer, parameter :: inml = 5,  & 
                        olog = 10, &
                        icnt = 20, &
                        ophi = 30, &
                        oAl  = 31, &
                        omom = 32, &
                        otrn = 33, &
                        otri = 34, &
                        ofxv = 40, &
                        ocnt = 50, &
                        odtc = 59, &
                        oeng = 60, &
                        omen = 61, &
                        owes = 62, &
                        owem = 63, &
                        oges = 64, &
                        ogem = 65, &
                        oqes = 66, &
                        oqem = 67, &
                        obln = 68, &
                        ofrq = 69, &
                        odsp = 70, &
                        ocst = 71, &
                        inbz = 14, &
                        ivmc = 15, &
                        omtr = 16, &
                        ovmc = olog


END MODULE GKV_header
