MODULE GKV_set
!-------------------------------------------------------------------------------
!
!    Set file I/O, and read parameters from namelist
!
!    Update history of gkvp_set.f90
!    --------------
!      gkvp_f0.61 (S. Maeyama, Mar 2021)
!        - equib_type = "s-alpha-shift" is added.
!        - Initial random phase rr is set by global (mx,gmy) indices.
!      gkvp_f0.60 (S. Maeyama, Jan 2021)
!        - Use fileio module to switch Fortran/NetCDF binary output.
!      gkvp_f0.58 (S. Maeyama, Oct 2020)
!        - init_random is added to switch random number for initialization.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - equib_type="slab" is added for shearless slab geometry.
!        - Set ky=0, ksq=0 for padding iend_y<my.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_math,   only: math_j0, math_j1, math_j2, math_g0, math_random
  use GKV_intgrl, only: intgrl_fsrf, intgrl_v0_moment_ms
  use GKV_fld,    only: fld_esfield, fld_emfield_ff, fld_ff2hh
  use GKV_bndry,  only: bndry_zvm_bound_f
  use GKV_advnc,  only: caldlt_rev
  use GKV_dtc,    only: dtc_init
! for vmec equilibrium
! use GKV_vmecin, only: vmecin_fileopen, vmecin_coeff, vmecin_read
! for vmec equilibrium w/ Booz_xform by M. Nakata & M. Nunami  (Aug. 2016)
  use GKV_vmecbzx, only: vmecbzx_boozx_read, vmecbzx_boozx_coeff
! for tokamak(eqdsk) equilibrium
  use GKV_igs,    only: igs_read, igs_coeff
  use GKV_colli,  only: colli_set_param
  use GKV_colliimp,  only: colliimp_set_param
  use GKV_tips,   only: tips_reality
  !fj start 202010
  use GKV_fileio
  !fj end 202010
  use GKV_shearflow, only: shearflow_lagrange_init, & !call after reading "time"
                           shearflow_running_init !call after reading "time"

  implicit none

  private

  public   set_init, set_close


CONTAINS

!--------------------------------------
  SUBROUTINE set_init( ff, phi, Al, hh, time )
!--------------------------------------

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    real(kind=DP), intent(out) :: time


      call set_start
      call set_param


      if ( trim(equib_type) == "slab"     .OR. &
           trim(equib_type) == "analytic" .OR. &
           trim(equib_type) == "s-alpha"  .OR. &
           trim(equib_type) == "s-alpha-shift"  .OR. &
           trim(equib_type) == "circ-MHD" .OR. &
           trim(equib_type) == "vmec"     .OR. &
           trim(equib_type) == "eqdsk" ) then

        call set_cnfig

      else

        if ( rankg == 0 ) then
          write(*,*) "set_cnfig_error!! on namelist: equib"
        end if
        call MPI_Finalize (ierr_mpi)
        stop

      end if

      call set_value( ff, phi, Al, hh, time )

    return


  END SUBROUTINE set_init


!--------------------------------------
  SUBROUTINE set_start
!--------------------------------------

    character(128) :: memo

    character(6)   :: crank
    character(1)   :: srank
    character(3)   :: cold, cnew

    character(10)   :: cdate, ctime

    namelist /cmemo/ memo
    namelist /calct/ calc_type, z_bound, z_filt, z_calc, art_diff, &
                     init_random, num_triad_diag
    namelist /equib/ equib_type
    namelist /run_n/ inum, ch_res
    namelist /files/ f_log, f_hst, f_phi, f_fxv, f_cnt

    character(256)   :: env_string       !fj


      call getenv ( 'fu05',env_string )  !fj
      open(inml, file=env_string )       !fj


      call date_and_time( cdate, ctime )

      read(inml,nml=cmemo)


      read(inml,nml=calct)
      read(inml,nml=equib)
      if (trim(z_calc) == "up5") art_diff = 0._DP


      inum = 1
      ch_res = .false.
      read(inml,nml=run_n)


      read(inml,nml=files)

      write( crank, fmt="(i6.6)" ) rankg
      write( srank, fmt="(i1.1)" ) ranks
      write( cold,  fmt="(i3.3)" ) inum-1
      write( cnew,  fmt="(i3.3)" ) inum


      open( olog, file=trim(f_log)//crank//"."//srank//".log."//cnew )

      if ( inum > 1 ) then
        !fj start 202010
        !open( icnt, file=trim(f_cnt)//crank//".cnt."//cold, &
        !      form="unformatted", status="old", action="read" )
        call fileio_open_icnt( trim(f_cnt) )
        !fj end 202010
      end if

      !fj start 202010
      !open( ofxv, file=trim(f_fxv)//crank//"."//srank//".fxv."//cnew, form="unformatted" )
      !open( ocnt, file=trim(f_cnt)//crank//".cnt."//cnew, form="unformatted" )
      call fileio_open_fxv( trim(f_fxv) )
      call fileio_open_cnt( trim(f_cnt) )
      !fj end 202010

      !fj start 202011
      !if ( vel_rank == 0 ) then
      !  open( omom, file=trim(f_phi)//crank//"."//srank//".mom."//cnew, form="unformatted" )
      !end if
      call fileio_open_mom( trim(f_phi) )
      !fj end 202011

      !fj start 202011
      !if ( ranks == 0 .AND. vel_rank == 0 ) then
      !  open( ophi, file=trim(f_phi)//crank//"."//srank//".phi."//cnew, form="unformatted" )
      !  open(  oAl, file=trim(f_phi)//crank//"."//srank//".Al."//cnew, form="unformatted" )
      !end if
      call fileio_open_phi( trim(f_phi) )
      call fileio_open_Al( trim(f_phi) )
      !fj end 202011

      !fj start 202011
      !if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
      !  open( otrn, file=trim(f_phi)//crank//"."//srank//".trn."//cnew, form="unformatted" )
      !end if
      call fileio_open_trn( trim(f_phi) )
      !fj end 202011

      if( rankg == 0 ) then
        open( omtr, file=trim(f_hst)//"mtr."//cnew )
        open( odtc, file=trim(f_hst)//"dtc."//cnew )
        open( oeng, file=trim(f_hst)//"eng."//cnew )
        open( omen, file=trim(f_hst)//"men."//cnew )
        open( owes, file=trim(f_hst)//"wes."//cnew )
        open( owem, file=trim(f_hst)//"wem."//cnew )
        if ( trim(calc_type) == "lin_freq" ) then
          open( ofrq, file=trim(f_hst)//"frq."//cnew )
          open( odsp, file=trim(f_hst)//"dsp."//cnew )
        end if
      end if
      if( rankg == nprocz/2 ) then
        open( ocst, file=trim(f_hst)//"cst."//cnew )
      end if

      if( rank == 0 ) then
        open( obln, file=trim(f_hst)//"bln."//srank//"."//cnew )
        open( oges, file=trim(f_hst)//"ges."//srank//"."//cnew )
        open( ogem, file=trim(f_hst)//"gem."//srank//"."//cnew )
        open( oqes, file=trim(f_hst)//"qes."//srank//"."//cnew )
        open( oqem, file=trim(f_hst)//"qem."//srank//"."//cnew )
      end if

      write( olog, * ) "##### ", trim(memo), " #####"
      write( olog, * ) ""
      write( olog, * ) "# Date : ", cdate
      write( olog, * ) "# Time : ", ctime
      write( olog, * ) ""
      write( olog, * ) "# Type of calc. : ", trim(calc_type)
      write( olog, * ) "# Boundary condition in zz       : ", trim(z_bound)
      write( olog, * ) "# 4th-order filter in zz         : ", trim(z_filt)
      write( olog, * ) "# Finite difference scheme in zz : ", trim(z_calc)
      write( olog, * ) "# Artificial diffusion in zz     : ", art_diff
      write( olog, * ) "# Number of triad transfer diag. : ", num_triad_diag
      write( olog, * ) "# Type of equib. : ", trim(equib_type)
      write( olog, * ) ""
      write( olog, * ) "# Run number = ", inum
      write( olog, * ) "# Resolution change = ", ch_res
      write( olog, * ) ""


    return
  

  END SUBROUTINE set_start


!--------------------------------------
  SUBROUTINE set_close
!--------------------------------------

     close( olog )

     !fj start 202010
     !close( icnt )
     !close( ofxv )
     !close( ocnt )
     call fileio_close_icnt
     call fileio_close_fxv
     call fileio_close_cnt
     !fj end 202010

     !fj start 202011
     !if ( vel_rank == 0 ) then
     !  close( omom )
     !end if
     call fileio_close_mom
     !fj end 202011

     !fj start 202011
     !if ( ranks == 0 .AND. vel_rank == 0 ) then
     !  close( ophi )
     !  close( oAl  )
     !end if
     call fileio_close_phi
     call fileio_close_Al
     !fj end 202011

     !fj start 202011
     !if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
       !close( otrn )
     !end if
     call fileio_close_trn
     !fj end 202011

     if( rankg == 0 ) then
       close( omtr )
       close( odtc )
       close( oeng )
       close( omen )
       close( owes )
       close( owem )
       if ( trim(calc_type) == "lin_freq" ) then
         close( ofrq )
         close( odsp )
       end if
     end if
      if( rankg == nprocz/2 ) then
       close( ocst )
     end if

     if( rank == 0 ) then
       close( obln )
       close( oges )
       close( ogem )
       close( oqes )
       close( oqem )
     end if

  END SUBROUTINE set_close


!--------------------------------------
  SUBROUTINE set_param
!--------------------------------------

    namelist /runlm/ e_limit
    namelist /times/ tend, dtout_fxv, dtout_ptn, dtout_eng, dtout_dtc
    namelist /deltt/ dt_max, adapt_dt, courant_num, time_advnc


      e_limit   = 5._DP * 3600._DP - 300._DP

      read(inml,nml=runlm)


      tend      = 10.00_DP
      dtout_fxv = 5._DP
      dtout_ptn = 0.5_DP
      dtout_eng = 0.005_DP
      dtout_dtc = 0.005_DP

      read(inml,nml=times)


      dt_max    = 0.005_DP
      adapt_dt  = .false.
      courant_num  = 0.6_DP

      read(inml,nml=deltt)


        write( olog, * ) " # Numerical parameters "
        write( olog, * ) ""
        write( olog, * ) " # nxw, nyw  = ", nxw, nyw
        write( olog, * ) " # global_ny = ", global_ny
        write( olog, * ) " # global_nz = ", global_nz
        write( olog, * ) " # global_nv, global_nm = ", global_nv, global_nm
        write( olog, * ) ""
        write( olog, * ) " # nx, ny, nz   = ", nx, ny, nz
        write( olog, * ) " # nv, nm       = ", nv, nm
        write( olog, * ) " # nzb, nvb     = ", nzb, nvb
        write( olog, * ) " # nxw_size (local xw allocation size) = ", nxw_size
        write( olog, * ) " # number of species  = ", nprocs
        write( olog, * ) " # ranks=0: Electron"
        write( olog, * ) " # ranks=1: main ion"
        write( olog, * ) " # ranks>1: other ions"
        write( olog, * ) " # Note that proton mass mp and main ion tmep. Ti is used for normalizations"
        write( olog, * ) " # kx, ky are normalized with rho_tp = mp*vtp/e/Baxi, where vtp = sqrt(Ti/mp)"
        write( olog, * ) " # time t is normalized with Raxi/vtp"
        write( olog, * ) ""
        write( olog, * ) " # e_limit      = ", e_limit
        write( olog, * ) " # tend         = ", tend
        write( olog, * ) " # dtout_fxv, dtout_ptn = ", dtout_fxv, dtout_ptn
        write( olog, * ) " # dtout_eng, dtout_dtc = ", dtout_eng, dtout_dtc
        write( olog, * ) ""
        write( olog, * ) " # maximum time step dt_max = ", dt_max
        write( olog, * ) " # adaptive time step = ", adapt_dt
        write( olog, * ) ""

        write( olog, * ) " # MPI parallelization parameters "
        write( olog, * ) ""
        write( olog, * ) " # nproc , rankg = ", nproc , rankg
        write( olog, * ) " # nprocw, rankw = ", nprocw, rankw
        write( olog, * ) " # nprocz, rankz = ", nprocz, rankz
        write( olog, * ) " # nprocv, rankv = ", nprocv, rankv
        write( olog, * ) " # nprocm, rankm = ", nprocm, rankm
        write( olog, * ) " # nprocs, rank  = ", nprocs , rank
        write( olog, * ) " # izup, izdn    = ", izup, izdn
        write( olog, * ) " # ivup, ivdn    = ", ivup, ivdn
        write( olog, * ) " # imup, imdn    = ", imup, imdn
        write( olog, * ) ""
        write( olog, * ) " # fft_nproc , fft_rank  = ", fft_nproc , fft_rank
        write( olog, * ) " # zsp_nproc , zsp_rank  = ", zsp_nproc , zsp_rank
        write( olog, * ) " # vel_nproc , vel_rank  = ", vel_nproc , vel_rank
        write( olog, * ) " # ranks                 = ", ranks
        write( olog, * ) ""

        write( olog, * ) " # ist_y     = ", ist_y
        write( olog, * ) " # iend_y    = ", iend_y
        write( olog, * ) " # nsize_y   = ", nsize_y
        write( olog, * ) " # ist1_y    = ", ist1_y
        write( olog, * ) " # ist_y_g   = ", ist_y_g
        write( olog, * ) " # iend_y_g  = ", iend_y_g
        write( olog, * ) ""
        write( olog, * ) " # ist_xw    = ", ist_xw
        write( olog, * ) " # iend_xw   = ", iend_xw
        write( olog, * ) " # nsize_xw  = ", nsize_xw
        write( olog, * ) " # ist_xw_g  = ", ist_xw_g
        write( olog, * ) " # iend_xw_g = ", iend_xw_g
        write( olog, * ) ""


  END SUBROUTINE set_param


!--------------------------------------
  SUBROUTINE set_cnfig
!--------------------------------------

    real(kind=DP) :: s_hat

    real(kind=DP) :: eps_r

    real(kind=DP) :: rdeps00, eps_hor, lprd, mprd, lmmq, malpha
    real(kind=DP) :: eps_mor, eps_por, lprdm1, lprdp1, lmmqm1, lmmqp1
    real(kind=DP) :: eps_rnew, rdeps1_0, rdeps1_10, rdeps2_10, rdeps3_10

! for s-alpha model with Shafranov shift
    real(kind=DP) ::  p_total, dp_totaldx, beta_total, alpha_MHD

! for circular MHD
    real(kind=DP), dimension(1:3,1:3) :: gg
    real(kind=DP) :: kkx, kky, domgdz, domgdx, domgdy


!! for vmec equilibrium
!    real(kind=DP) :: rho2R_0, q_input, theta
!    real(kind=DP) :: r_0
!    real(kind=DP) :: gdwss, gdwtt, gdwzz, gdwst, gdwsz, gdwtz, &
!                     gupss, guptt, gupzz, gupst, gupsz, guptz, &
!                     babs, Bs  , Bth , Bzt , dBds, dBdt, dBdz, &
!                     dBdt_mir, vmec_rootg, rootgft, rootgbz
    real(kind=DP) :: theta


    real(kind=DP) :: lx, ly, lz, kxmin, kymin, dz, mmax, dm, del_c
    real(kind=DP) :: lz_l, z0, z0_l
    integer       :: n_tht, m_j

    real(kind=DP) :: gg0

    real(kind=DP) :: bb, kmo
    real(kind=DP) :: cfsrf_l

    integer       :: global_iv, global_im
    integer       :: mx, my, iz, iv, im, is, is1, is2, ierr_mpi


    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: nw
    real(kind=DP), dimension(:,:,:), allocatable :: ww

!    real(kind=DP) :: rad_a, r_minor, eps_b, rho_unit, r_a
!    real(kind=DP) :: R0_unit, r_edge, b0b00, alpha_fix

    real(kind=DP), dimension(0:ns-1) :: eta
    real(kind=DP), dimension(0:ns-1,0:ns-1) :: nust
    real(kind=DP) :: r_major
   
    real(kind=DP) :: s_input, s_0      ! radial label of fluxtube center 
    integer       :: mc_type           ! 0:Axisym., 1:Boozer, 2:Hamada
    integer       :: q_type            ! 0:use q and s_hat value in confp, 1:calclated by IGS
    integer       :: isw, nss, ntheta, nzeta
    real(kind=DP) :: phi_ax            ! axisymetric toroidal angle 

    integer, parameter :: num_omtr = 13
    real(kind=DP) :: metric_l(1:num_omtr,-nz:nz-1), metric_g(1:num_omtr,-global_nz:global_nz-1)




    namelist /physp/ R0_Ln,  &    ! R0/Lns
                     R0_Lt,  &    ! R0/Lts
                        nu,  &    ! factor for collision freq. in LB model    
                      Anum,  &    ! mass number
                      Znum,  &    ! charge number 
                       fcs,  &    ! charge-density fraction 
                       sgn,  &    ! signs of charge 
                       tau,  &    ! T-ratio Ts/T0, T0=reference ion temp. of ranks=1
                      dns1,  &    ! initial perturbation amplitude
                    tau_ad,  &    ! Ti/Te for ITG-ae, Te/Ti for ETG-ai
                  lambda_i,  &    ! (Debye/rho_tp)^2 
                      beta,  &    ! mu0*ni*Ti/B^2
                   ibprime,  &    ! flag for finite beta-prime effect on kvd
                      vmax,  &    ! maximum v_para in unit of v_ts
                       nx0        ! mode number for the initial perturbation

    namelist /rotat/ mach, uprime, gamma_e

    namelist /nperi/ n_tht, kymin, m_j, del_c
    namelist /confp/ eps_r, eps_rnew,                       &
                     q_0, s_hat,                            &
                     lprd, mprd, eps_hor, eps_mor, eps_por, &
                     rdeps00, rdeps1_0, rdeps1_10,          & 
                     rdeps2_10, rdeps3_10, malpha
!    namelist /vmecp/ q_0, rad_a,                            &
!                     R0_unit, r_edge,                       &
!                     b0b00, alpha_fix
    namelist /vmecp/ s_input, nss, ntheta, nzeta

    namelist /igsp/ s_input, mc_type, q_type, nss, ntheta

    namelist /nu_ref/ Nref,     & ! reference (electron) density in m^(-3)
                      Lref,     & ! reference length (=R_axis) in m
                      Tref,     & ! reference main-ion (ranks=1) temperature in keV 
                      col_type, & ! flag for collision type: LB or full
                      iFLR,     & ! flag for GK- or DK-limit in collision
                      icheck      ! flag for Maxwellain anihilation test (w/ iFLR=0)

      tau(:)   = 1.0_DP
      nu(:)    = 0.002_DP
      R0_Ln(:) = 2.5_DP
      R0_Lt(:) = 7.5_DP


      read(inml,nml=physp)


        do is = 0, ns-1
          if( R0_Ln(is) /= 0._DP ) then
            eta(is) = R0_Lt(is) / R0_Ln(is)
          else
            eta(is) = 1.d+20
          end if
        end do


        write( olog, * ) " # Physical parameters"
        write( olog, * ) ""
        write( olog, * ) " # r_major/L_ns = ", R0_Ln(:)
        write( olog, * ) " # r_major/L_ts = ", R0_Lt(:)
        write( olog, * ) " # eta          = ", eta(:)
        write( olog, * ) " # nu           = ", nu(:)
        write( olog, * ) " # A-number     = ", Anum(:)
        write( olog, * ) " # Z-number     = ", Znum(:)
        write( olog, * ) " # fcs          = ", fcs(:)
        write( olog, * ) " # sgn          = ", sgn(:)
        write( olog, * ) " # tau          = ", tau(:)
        write( olog, * ) " # dns1         = ", dns1(:)
        write( olog, * ) " # tau_ad       = ", tau_ad
        write( olog, * ) " # lambda_i^2   = ", lambda_i
        write( olog, * ) " # beta_i       = ", beta
        write( olog, * ) " # ibprime      = ", ibprime
        write( olog, * ) " # nx0          = ", nx0
        write( olog, * ) ""


      mach = 0._DP
      uprime = 0._DP
      gamma_e = 0._DP

      read(inml,nml=rotat)

        write( olog, * ) " # Mean rotation parameters"
        write( olog, * ) ""
        write( olog, * ) " # Mach number  = ", mach
        write( olog, * ) " # uptime       = ", uprime
        write( olog, * ) " # gamma_ExB    = ", gamma_e
        write( olog, * ) ""


      n_tht = 1

      read(inml,nml=nperi)


      if( trim(equib_type) == "slab") then

        read(inml,nml=confp)

        lprdm1   = 0._DP
        lprdp1   = 0._DP

        lmmq     = 0._DP
        lmmqm1   = 0._DP
        lmmqp1   = 0._DP

        q_0      = 1._DP ! For now, fixed q_0=1. Changing q_0 can extend parallel z-box size.
        s_hat    = 0._DP ! only shear less slab
        eps_r    = 1._DP

        eps_hor  = 0._DP
        lprd     = 0._DP
        mprd     = 0._DP
        malpha   = 0._DP

        rdeps00  = 0._DP
        eps_mor  = 0._DP
        eps_por  = 0._DP

        write( olog, * ) " # Configuration parameters"
        write( olog, * ) ""
        write( olog, * ) " # q_0          = ", q_0
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) ""

        write( olog, * ) " # eps_hor      = ", eps_hor
        write( olog, * ) " # lprd         = ", lprd
        write( olog, * ) " # mprd         = ", mprd
        write( olog, * ) " # malpha       = ", malpha
        write( olog, * ) " # rdeps00      = ", rdeps00

        write( olog, * ) " # eps_mor      = ", eps_mor
        write( olog, * ) " # lprdm1       = ", lprdm1
        write( olog, * ) " # eps_por      = ", eps_por
        write( olog, * ) " # lprdp1       = ", lprdp1
        write( olog, * ) ""

      else if( trim(equib_type) == "analytic"  .OR.  &
               trim(equib_type) == "s-alpha"   .OR.  &
               trim(equib_type) == "s-alpha-shift"   .OR.  &
               trim(equib_type) == "circ-MHD" ) then


        read(inml,nml=confp)


        lprdm1   = lprd - 1.0_DP
        lprdp1   = lprd + 1.0_DP

        lmmq     = lprd   - mprd * q_0
        lmmqm1   = lprdm1 - mprd * q_0
        lmmqp1   = lprdp1 - mprd * q_0


        write( olog, * ) " # Configuration parameters"
        write( olog, * ) ""
        write( olog, * ) " # q_0          = ", q_0
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) ""

        write( olog, * ) " # eps_hor      = ", eps_hor
        write( olog, * ) " # lprd         = ", lprd
        write( olog, * ) " # mprd         = ", mprd
        write( olog, * ) " # malpha       = ", malpha
        write( olog, * ) " # rdeps00      = ", rdeps00

        write( olog, * ) " # eps_mor      = ", eps_mor
        write( olog, * ) " # lprdm1       = ", lprdm1
        write( olog, * ) " # eps_por      = ", eps_por
        write( olog, * ) " # lprdp1       = ", lprdp1
        write( olog, * ) ""



!      else if( trim(equib_type) == "vmec" ) then
!
!
!! --- Paramters at rho=0.65 (shot#088343 at t = 1.833 [s])
!!
!!        Ln_unit =-4.230701_DP                      ! Ln [m]
!!        Lt_unit = 0.3135611_DP                     ! Lt [m]
!!        R0_unit = 3.599858_DP                      ! R0 [m]
!!        r_edge  = 0.6362872D0                      ! r_edge [m]
!!        b0b00   = 2.88846853946973647d0/2.940307D0 ! b00mode/B0
!!        alpha_fix = 0.314159253589793d0            ! pi/10
!!
!
!        read(inml,nml=vmecp)
!
!        eps_b      = r_edge / R0_unit   ! --- a / R                  ! by nunami (2010.04.21)
!
!        rho2R_0  = eps_b / rad_a               ! --- rho / R_0
!        rho_unit = rho2R_0 * R0_unit           ! --- rho
!        r_a      = rad_a * rho2R_0             ! --- rad_a / R_0
!
!        eps_r = 0.1115200537d0
!
!        call vmecin_fileopen
!
!
!        call vmecin_read
!
!
!        q_input = q_0
!        theta   = 0._DP
!
!        call vmecin_coeff( rad_a, R0_unit, rho2R_0, q_input, theta,  &
!                           alpha_fix, r_0, r_minor, s_hat,           &
!                           gdwss, gdwtt, gdwzz, gdwst, gdwsz, gdwtz, &
!                           gupss, guptt, gupzz, gupst, gupsz, guptz, &
!                           babs, Bs, Bth, Bzt, dBds, dBdt, dBdz,     &
!                           dBdt_mir, vmec_rootg, rootgft, rootgbz )
!
!
!        write( olog, * ) " # Configuration parameters"
!        write( olog, * ) ""
!        write( olog, * ) " # q_0          = ", q_0
!        write( olog, * ) " # s_hat        = ", s_hat
!        write( olog, * ) ""
!        write( olog, * ) " # eps_r        = ", eps_r
!        write( olog, * ) ""


     else if( trim(equib_type) == "vmec" ) then


        read(inml,nml=confp)

        read(inml,nml=vmecp)

        call vmecbzx_boozx_read( nss, ntheta, nzeta )

         isw = 0
         iz = 0
         call vmecbzx_boozx_coeff( isw,  nss,  ntheta,  nzeta,  s_input, iz, 0._DP,  lz_l,   &  ! input 
                           s_0,           q_0,     s_hat,    eps_r,  phi_ax,             &  ! output
                           omg(iz), rootg(iz),    domgdx,   domgdz,  domgdy,             &
                           gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),                      &
                           gg(2,3),   gg(3,3)  )



        write( olog, * ) " # Configuration parameters"
        write( olog, * ) ""
        write( olog, * ) " # r_major/L_ns = ", R0_Ln(:)
        write( olog, * ) " # r_major/L_ts = ", R0_Lt(:)
        write( olog, * ) " # eta          = ", eta(:)
        write( olog, * ) " # q_0          = ", q_0
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) " # s_input, s_0 = ", s_input, s_0
        write( olog, * ) " # nss, ntheta, nzeta  = ", nss, ntheta, nzeta


     else if( trim(equib_type) == "eqdsk" ) then


        read(inml,nml=confp)

        read(inml,nml=igsp)

        call igs_read( mc_type, nss, ntheta )

        if ( q_type == 1 ) then
         isw = 0
         iz = 0
         call igs_coeff( isw,  mc_type,   nss,    ntheta,  s_input,  0._DP, lz_l,   &  ! input 
                         s_0,       q_0,     s_hat,    eps_r,  theta,               &  ! output
                           omg(iz), rootg(iz),    domgdx,   domgdz, domgdy,         &
                           gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),                 &
                           gg(2,3),   gg(3,3)  )
        end if



        write( olog, * ) " # Configuration parameters"
        write( olog, * ) ""
        write( olog, * ) " # r_major/L_ns = ", R0_Ln(:)
        write( olog, * ) " # r_major/L_ts = ", R0_Lt(:)
        write( olog, * ) " # eta          = ", eta(:)
        write( olog, * ) " # q_0          = ", q_0
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) " # s_input, s_0 = ", s_input, s_0
        write( olog, * ) " # nss, ntheta  = ", nss, ntheta

      else

        write( olog, * ) " # wrong choice of the equilibrium "
        call flush(olog)
        call MPI_Finalize(ierr_mpi)
        stop

      end if


! --- coordinate settings ---

        if (abs(s_hat) < 1.d-10) then ! When s_hat == ZERO
          m_j = 0
          kxmin = kymin
        else if (m_j == 0) then
          kxmin = kymin
        else
          kxmin    = abs(2._DP * pi * s_hat * kymin / real(m_j, kind=DP))
        end if
        lx       = pi / kxmin
        ly       = pi / kymin
                      ! kymin=pi/ly=pi/[r_minor*pi/(q0*n_alp)]=q0*n_alp/r_minor

        lz       = real( n_tht, kind=DP ) * pi        ! total z-length
        lz_l     = lz / real( nprocz, kind=DP )       ! local z-length

        do mx = -nx, nx
          kx(mx)   = kxmin * real( mx, kind=DP )
        end do

        ky(:) = 0._DP
        do my = ist_y_g, iend_y_g
          ky(my-ist_y_g)   = kymin * real( my, kind=DP )
        end do

        kxmin_g = kxmin
        kymin_g = kymin 

        z0       = - lz                     ! global lower boundary
        z0_l     = 2._DP * lz_l * real( rankz, kind=DP ) + z0
                                            ! local lower boundary

        dz       = lz_l / real( nz, kind=DP )

        do iz = -nz, nz-1
          zz(iz)   = dz * real( iz + nz, kind=DP ) + z0_l
        end do


        dv       = 2._DP * vmax / real( 2 * nv * nprocv -1, kind=DP )

        do iv = 1, 2*nv
          global_iv = 2 * nv * rankv + iv
          vl(iv)    = dv * ( real( global_iv - nv * nprocv - 1, kind=DP ) + 0.5_DP )
        end do
                                            ! --- debug
                                            !   write( olog, * ) " *** iv, vl "
                                            !   do iv = 1, 2*nv
                                            !     global_iv = 2 * nv * rankv + iv
                                            !     write( olog, * ) iv, global_iv, vl(iv)
                                            !   end do
                                            !   write( olog, * ) ""

        mmax     = vmax
        dm       = mmax / real( nprocm * ( nm+1 ) - 1, kind=DP )
                                            ! --- equal spacing in vperp

        do im = 0, nm
          global_im = ( nm+1 ) * rankm + im
          mu(im)    = 0.5_DP * ( dm * real( global_im, kind=DP ) )**2
        end do


        write( olog, * ) " # Numerical parameters"
        write( olog, * ) ""
        write( olog, * ) " # n_tht = ", n_tht
        write( olog, * ) " # lx, ly, lz   = ", lx, ly, lz
        write( olog, * ) " # lz,   z0     = ", lz, z0
        write( olog, * ) " # lz_l, z0_l   = ", lz_l, z0_l
        write( olog, * ) " # kxmin, kymin = ", kxmin, kymin
        write( olog, * ) " # kxmax, kymax = ", kxmin*nx, kymin*global_ny
        write( olog, * ) " # kperp_max    = ", sqrt((kxmin*nx)**2+(kymin*global_ny)**2)
        write( olog, * ) " # m_j, del_c   = ", m_j, del_c
        write( olog, * ) " # dz           = ", dz
        write( olog, * ) " # dv, vmax     = ", dv, vmax
        write( olog, * ) " # dm, mmax     = ", dm, mmax
        write( olog, * ) ""

        if (gamma_e == 0._DP) then
          tlim_exb = 999999.d0
        else
          tlim_exb = (kxmin*(nx-nx0))/(kymin*global_ny*abs(gamma_e))
        end if
        write( olog, * ) " # ExB limit time tlim_exb  = ", tlim_exb
        write( olog, * ) " # for (mx=nx0,my=global_ny) initial perturbation: "
        write( olog, * ) " # tlim_exb = kxmin*(nx-nx0)/(kymax*|gamma_e|)"
        write( olog, * ) ""


! --- coordinate settings ---

 
! --- operator settings ---


        do iz = -nz, nz-1

!!! for slab model
          if ( trim(equib_type) == "slab") then

            q_bar = q_0 
            r_major = 1._DP ! in the R0 unit
            theta = zz(iz)

            omg(iz)   = 1._DP
            rootg(iz) = q_0*r_major
            dpara(iz) = dz * q_0 * r_major

            do im = 0, nm
              vp(iz,im)  = sqrt( 2._DP * mu(im) )!* omg(iz) )
              mir(iz,im) = 0._DP
              do iv = 1, 2*nv
                vdx(iz,iv,im) = 0._DP
                vdy(iz,iv,im) = 0._DP
                vsy(iz,iv,im) =                                           &
                  - sgn(ranks) * tau(ranks) / Znum(ranks)                 & 
                  * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                           + omg(iz)*mu(im) - 1.5_DP ) )
              end do
            end do   ! im loop ends

            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = kx(mx)**2 + ky(my)**2
              end do
            end do

            baxfactor    = 1._DP

           !- for OUTPUT hst/*.mtr.* -
            domgdz = 0._DP
            domgdy = 0._DP
            domgdx = 0._DP
            gg(1,1) = 1._DP
            gg(1,2) = 0._DP
            gg(1,3) = 0._DP
            gg(2,1) = gg(1,2)
            gg(2,2) = 1._DP
            gg(2,3) = 0._DP
            gg(3,1) = gg(1,3)
            gg(3,2) = gg(2,3)
            gg(3,3) = 1._DP
            metric_l( 1,iz) = zz(iz)   ! [ 1]
            metric_l( 2,iz) = theta    ! [ 2]
            metric_l( 3,iz) = omg(iz)  ! [ 3]
            metric_l( 4,iz) = domgdx   ! [ 4]
            metric_l( 5,iz) = domgdy   ! [ 5]
            metric_l( 6,iz) = domgdz   ! [ 6]
            metric_l( 7,iz) = gg(1,1)  ! [ 7]
            metric_l( 8,iz) = gg(1,2)  ! [ 8]
            metric_l( 9,iz) = gg(1,3)  ! [ 9]
            metric_l(10,iz) = gg(2,2)  ! [10]
            metric_l(11,iz) = gg(2,3)  ! [11]
            metric_l(12,iz) = gg(3,3)  ! [12]
            metric_l(13,iz) = rootg(iz)! [13]
           !-------------------------


!!! for the concentric and large-aspect-ratio model !!!
          else if( trim(equib_type) == "analytic" ) then

            q_bar = q_0 
            r_major = 1._DP ! in the R0 unit

            theta = zz(iz)

            omg(iz)   = 1._DP                                             &
                      - eps_r * ( cos( zz(iz) )                           &
                              + eps_hor * cos( lmmq   * zz(iz) - malpha ) &
                              + eps_mor * cos( lmmqm1 * zz(iz) - malpha ) &
                              + eps_por * cos( lmmqp1 * zz(iz) - malpha ) )

            rootg(iz) = q_0*r_major/omg(iz)
            dpara(iz) = dz * q_0 * r_major

                                            ! --- debug
                                            !    write( olog, * ) " *** z, omg "
                                            !    do iz = -nz, nz-1
                                            !      write( olog, * ) zz(iz), omg(iz)
                                            !    end do
                                            !    write( olog, * )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )
              mir(iz,im) = mu(im) * eps_r / ( q_0 * r_major )  &
                       * ( sin(zz(iz))                                          &
                         + eps_hor * lmmq   * sin( lmmq   * zz(iz) - malpha )   &
                         + eps_mor * lmmqm1 * sin( lmmqm1 * zz(iz) - malpha )   &
                         + eps_por * lmmqp1 * sin( lmmqp1 * zz(iz) - malpha ) )

              do iv = 1, 2*nv
                !do my = ist_y, iend_y
                !  do mx = -nx, nx
                                            ! kvd and kvs are revised November 2011
                                            !             into general species forms. 
                    !!!!kvd(mx,my,iz,iv,im)=                                            &
                    !!!!  - ( vl(iv)**2 + omg(iz)*mu(im) ) * eps_rnew / r_major         &
                    !!!!  * ( ky(my) * ( rdeps00 + rdeps1_0 * cos( zz(iz) )             &
                    !!!!         + rdeps2_10 * cos( lmmq   * zz(iz) - malpha )          &
                    !!!!         + rdeps1_10 * cos( lmmqm1 * zz(iz) - malpha )          &
                    !!!!         + rdeps3_10 * cos( lmmqp1 * zz(iz) - malpha ) )        &
                    !!!!         + ( kx(mx) + s_hat * zz(iz) * ky(my) )                 &
                    !!!!         * ( sin( zz(iz) )                                      &
                    !!!!         + eps_hor * lprd   * sin( lmmq   * zz(iz) - malpha )   &
                    !!!!         + eps_mor * lprdm1 * sin( lmmqm1 * zz(iz) - malpha )   &
                    !!!!         + eps_por * lprdp1 * sin( lmmqp1 * zz(iz) - malpha ) ) &
                    !!!!     ) * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    !!!!kvs(mx,my,iz,iv,im) =                                     &
                    !!!!  - sgn(ranks) * tau(ranks) / Znum(ranks) * ky(my)        & 
                    !!!!  * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                    !!!!                           + omg(iz)*mu(im) - 1.5_DP ) )

                    !%%% kvd = kx(mx) * vdx(iz,iv,im) + ky(my) * vdy(iz,iv,im) %%%
                    !%%% kvs = ky(my) * vsy(iz,iv,im)                          %%%
                    vdx(iz,iv,im)=                                            &
                      - ( vl(iv)**2 + omg(iz)*mu(im) ) * eps_rnew / r_major         &
                      * ( 0._DP * ( rdeps00 + rdeps1_0 * cos( zz(iz) )             &
                             + rdeps2_10 * cos( lmmq   * zz(iz) - malpha )          &
                             + rdeps1_10 * cos( lmmqm1 * zz(iz) - malpha )          &
                             + rdeps3_10 * cos( lmmqp1 * zz(iz) - malpha ) )        &
                             + ( 1._DP + s_hat * zz(iz) * 0._DP )                 &
                             * ( sin( zz(iz) )                                      &
                             + eps_hor * lprd   * sin( lmmq   * zz(iz) - malpha )   &
                             + eps_mor * lprdm1 * sin( lmmqm1 * zz(iz) - malpha )   &
                             + eps_por * lprdp1 * sin( lmmqp1 * zz(iz) - malpha ) ) &
                         ) * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vdy(iz,iv,im)=                                            &
                      - ( vl(iv)**2 + omg(iz)*mu(im) ) * eps_rnew / r_major         &
                      * ( 1._DP * ( rdeps00 + rdeps1_0 * cos( zz(iz) )             &
                             + rdeps2_10 * cos( lmmq   * zz(iz) - malpha )          &
                             + rdeps1_10 * cos( lmmqm1 * zz(iz) - malpha )          &
                             + rdeps3_10 * cos( lmmqp1 * zz(iz) - malpha ) )        &
                             + ( 0._DP + s_hat * zz(iz) * 1._DP )                 &
                             * ( sin( zz(iz) )                                      &
                             + eps_hor * lprd   * sin( lmmq   * zz(iz) - malpha )   &
                             + eps_mor * lprdm1 * sin( lmmqm1 * zz(iz) - malpha )   &
                             + eps_por * lprdp1 * sin( lmmqp1 * zz(iz) - malpha ) ) &
                         ) * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 & 
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )

                !  end do
                !end do
              end do

            end do   ! im loop ends


            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = ( kx(mx) + s_hat * zz(iz) * ky(my) )**2 + ky(my)**2
              end do
            end do

            baxfactor    = 1._DP

           !- for OUTPUT hst/*.mtr.* - !%%% under benchmark %%%
            domgdz = eps_r * ( sin(zz(iz))                                          &
                             + eps_hor * lmmq   * sin( lmmq   * zz(iz) - malpha )   &
                             + eps_mor * lmmqm1 * sin( lmmqm1 * zz(iz) - malpha )   &
                             + eps_por * lmmqp1 * sin( lmmqp1 * zz(iz) - malpha ) )
            domgdy = - eps_rnew / r_major * (                                  &
                     - ( sin( zz(iz) )                                         &
                        + eps_hor * lprd   * sin( lmmq   * zz(iz) - malpha )   &
                        + eps_mor * lprdm1 * sin( lmmqm1 * zz(iz) - malpha )   &
                        + eps_por * lprdp1 * sin( lmmqp1 * zz(iz) - malpha )   &
                     ) - (-1._DP/eps_r) * domgdz )
            domgdx = eps_rnew / r_major * (                                    &
                     - (                                                       &
                          rdeps00                                              &
                        + rdeps1_0 * cos( zz(iz) )                             &
                        + rdeps2_10 * cos( lmmq   * zz(iz) - malpha )          &
                        + rdeps1_10 * cos( lmmqm1 * zz(iz) - malpha )          &
                        + rdeps3_10 * cos( lmmqp1 * zz(iz) - malpha )          &
                        + s_hat * zz(iz) * ( sin( zz(iz) )                     &
                        + eps_hor * lprd   * sin( lmmq   * zz(iz) - malpha )   &
                        + eps_mor * lprdm1 * sin( lmmqm1 * zz(iz) - malpha )   &
                        + eps_por * lprdp1 * sin( lmmqp1 * zz(iz) - malpha ) ) &
                     ) - (-s_hat*zz(iz)/eps_r) * domgdz )
            gg(1,1) = 1._DP
            gg(1,2) = s_hat*zz(iz)
            gg(1,3) = 0._DP
            gg(2,1) = gg(1,2)
            gg(2,2) = 1._DP + (s_hat*zz(iz))**2
            gg(2,3) = 1._DP/(r_major*eps_r)
            gg(3,1) = gg(1,3)
            gg(3,2) = gg(2,3)
            gg(3,3) = 1._DP/((r_major*eps_r)**2)
            metric_l( 1,iz) = zz(iz)   ! [ 1]
            metric_l( 2,iz) = theta    ! [ 2]
            metric_l( 3,iz) = omg(iz)  ! [ 3]
            metric_l( 4,iz) = domgdx   ! [ 4]
            metric_l( 5,iz) = domgdy   ! [ 5]
            metric_l( 6,iz) = domgdz   ! [ 6]
            metric_l( 7,iz) = gg(1,1)  ! [ 7]
            metric_l( 8,iz) = gg(1,2)  ! [ 8]
            metric_l( 9,iz) = gg(1,3)  ! [ 9]
            metric_l(10,iz) = gg(2,2)  ! [10]
            metric_l(11,iz) = gg(2,3)  ! [11]
            metric_l(12,iz) = gg(3,3)  ! [12]
            metric_l(13,iz) = rootg(iz)! [13]
           !-------------------------

!!! for s-alpha !!! <--- the current version is the same as "analytic"
          else if( trim(equib_type) == "s-alpha" .or. trim(equib_type) == "s-alpha-shift" ) then

            q_bar = q_0
            r_major = 1._DP ! in the R0 unit

            if (trim(equib_type) == "s-alpha") then
              !--- s-alpha model without Shafranov shift -
              alpha_MHD = 0._DP
            else if (trim(equib_type) == "s-alpha-shift") then
              !--- s-alpha model with Shafranov shift ----
              p_total = 0._DP
              dp_totaldx = 0._DP
              beta_total = 0._DP
              do is = 0, ns-1
                p_total = p_total + fcs(is) * tau(is) / Znum(is)
                dp_totaldx = dp_totaldx - fcs(is) * tau(is) / Znum(is) * (R0_Ln(is) + R0_Lt(is))
                beta_total = beta_total + 2._DP * beta * fcs(is) * tau(is) / Znum(is)
              end do
              alpha_MHD = - q_0**2 * r_major * beta_total * dp_totaldx / p_total
            end if

            theta = zz(iz)

            omg(iz)   = 1._DP - eps_r * cos( theta )        ! s-alpha with eps-expansion
            !omg(iz)   = 1._DP / (1._DP + eps_r * cos( theta )) ! for benchmark

            rootg(iz) = q_0*r_major/omg(iz)
            dpara(iz) = dz* q_0 * r_major

            domgdz  = eps_r * sin( theta )
            !domgdz  = eps_r * sin( theta ) * omg(iz)**2 ! for benchmark
            domgdx  = -cos( theta ) / r_major
            domgdy  = 0._DP


            gg(1,1) = 1._DP
            gg(1,2) = s_hat*zz(iz) - alpha_MHD*sin(zz(iz)) ! with Shafranov shift
            gg(1,3) = 0._DP
            gg(2,1) = gg(1,2)
            gg(2,2) = 1._DP + (s_hat*zz(iz) - alpha_MHD*sin(zz(iz)))**2 ! with Shafranov shift
            gg(2,3) = 1._DP/(r_major*eps_r)
            gg(3,1) = gg(1,3)
            gg(3,2) = gg(2,3)
            gg(3,3) = 1._DP/((r_major*eps_r)**2)

            kkx = -r_major * (q_0/q_bar) &
                           * ( gg(1,1)*gg(2,3) - gg(1,2)*gg(1,3) )*domgdz
            kky =  r_major * (q_bar/q_0) &
                           * ( domgdx - ( gg(1,2)*gg(2,3) - gg(2,2)*gg(1,3) )*domgdz )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              mir(iz,im) = mu(im) * (q_0/q_bar) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                !do my = ist_y, iend_y
                !  do mx = -nx, nx

                    !!!!kvd(mx,my,iz,iv,im) =                               &
                    !!!!   ( vl(iv)**2 + omg(iz)*mu(im) ) / r_major         &   
                    !!!!  * ( kkx*kx(mx) + kky*ky(my) )                     &
                    !!!!  * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    !!!!kvs(mx,my,iz,iv,im) =                                     &
                    !!!!  - sgn(ranks) * tau(ranks) / Znum(ranks) * ky(my)        &
                    !!!!  * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                    !!!!                           + omg(iz)*mu(im) - 1.5_DP ) )  &
                    !!!!  * (q_bar/q_0)

                    !%%% kvd = kx(mx) * vdx(iz,iv,im) + ky(my) * vdy(iz,iv,im) %%%
                    !%%% kvs = ky(my) * vsy(iz,iv,im)                          %%%
                    vdx(iz,iv,im) =                                     &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / r_major         &   
                      * kkx                                             &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vdy(iz,iv,im) =                                     &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / r_major         &   
                      * kky                                             &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )  &
                      * (q_bar/q_0)

                !  end do
                !end do
              end do

            end do   ! im loop ends


            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = ( kx(mx) + ( s_hat * zz(iz) - alpha_MHD*sin(zz(iz)) ) &
                                * ky(my) )**2 + ky(my)**2 ! with Shafranov shift
              end do
            end do

            baxfactor    = 1._DP

           !- for OUTPUT hst/*.mtr.* -
            metric_l( 1,iz) = zz(iz)   ! [ 1]
            metric_l( 2,iz) = theta    ! [ 2]
            metric_l( 3,iz) = omg(iz)  ! [ 3]
            metric_l( 4,iz) = domgdx   ! [ 4]
            metric_l( 5,iz) = domgdy   ! [ 5]
            metric_l( 6,iz) = domgdz   ! [ 6]
            metric_l( 7,iz) = gg(1,1)  ! [ 7]
            metric_l( 8,iz) = gg(1,2)  ! [ 8]
            metric_l( 9,iz) = gg(1,3)  ! [ 9]
            metric_l(10,iz) = gg(2,2)  ! [10]
            metric_l(11,iz) = gg(2,3)  ! [11]
            metric_l(12,iz) = gg(3,3)  ! [12]
            metric_l(13,iz) = rootg(iz)! [13]
           !-------------------------


!!! for circular MHD equilibrium !!!
          else if( trim(equib_type) == "circ-MHD" ) then

            q_bar = dsqrt( 1._DP - eps_r**2 )*q_0
            r_major = 1._DP ! in the R0 unit

            theta = 2._DP*atan( sqrt( (1._DP+eps_r)/(1._DP-eps_r) ) &
                              * tan(zz(iz)/2._DP) )

            omg(iz)   = sqrt( q_bar**2 + eps_r**2 ) &
                      / ( 1._DP + eps_r*cos( theta ) ) / q_bar

            rootg(iz) = q_0*r_major*( 1._DP+eps_r*cos(theta) )**2
            dpara(iz) = dz * omg(iz) * rootg(iz)


            domgdz  = eps_r * sin(theta) * sqrt( q_bar**2 + eps_r**2 ) &
                       / ( 1._DP + eps_r * cos( theta ) )**2           &
                       / ( 1._DP - eps_r * cos( zz(iz)) ) / q_0

            domgdx  = -( cos(theta)                                            &
                         - eps_r*(1._DP-s_hat+eps_r**2*q_0**2/q_bar**2)        &
                                *(1._DP+eps_r*cos(theta))/(q_bar**2+eps_r**2)  &
                         - eps_r*sin(theta)**2/(1._DP-eps_r**2)                &
                         ) / ((1._DP + eps_r*cos(theta))**2)                   &
                           * sqrt(q_bar**2+eps_r**2) / q_bar / r_major

            domgdy  = 0._DP

            gg(1,1) = (q_0/q_bar)**2
            gg(1,2) = ( s_hat*zz(iz)*q_0/q_bar - eps_r*sin(zz(iz))/(1._DP-eps_r**2) )*q_0/q_bar
            gg(1,3) = - sin(zz(iz))/(1._DP-eps_r**2)/r_major*q_0/q_bar
            gg(2,1) = gg(1,2)
            gg(2,2) = (s_hat*zz(iz)*q_0/q_bar)**2 - 2._DP*q_0/q_bar*s_hat*zz(iz)*eps_r*sin(zz(iz))/(1._DP-eps_r**2) &
                        + (q_bar**2+eps_r**2)/((1._DP+eps_r*cos(theta))**2)/(q_0**2)   &
                        + (eps_r*sin(zz(iz)))**2/(1._DP-eps_r**2)**2
            gg(2,3) = ( -s_hat*zz(iz)*q_0/q_bar*sin(zz(iz))/(1._DP-eps_r**2)           &
                           + ((q_bar/q_0)**2)/((1._DP+eps_r*cos(theta))**2)/eps_r      &
                           + eps_r*(sin(zz(iz))**2)/((1._DP-eps_r**2)**2)              &
                         ) / r_major
            gg(3,1) = gg(1,3)
            gg(3,2) = gg(2,3)
            gg(3,3) = ( ((q_bar/q_0)**2)/((1._DP+eps_r*cos(theta))**2)/(eps_r**2)  &
                           + (sin(zz(iz))**2)/((1._DP-eps_r**2)**2)                &
                         ) / (r_major**2)

            kkx =  r_major*( -domgdy + (gg(1,3)*gg(1,2) - gg(1,1)*gg(2,3))*domgdz/omg(iz)**2 )
            kky =  r_major*(  domgdx + (gg(1,3)*gg(2,2) - gg(1,2)*gg(2,3))*domgdz/omg(iz)**2 )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              mir(iz,im) = mu(im) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                !do my = ist_y, iend_y
                !  do mx = -nx, nx

                    !!!!kvd(mx,my,iz,iv,im)=                                        &
                    !!!!   ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )     &
                    !!!!  * ( kkx*kx(mx) + kky*ky(my) )                             &
                    !!!!  * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    !!!!kvs(mx,my,iz,iv,im) =                                     &
                    !!!!  - sgn(ranks) * tau(ranks) / Znum(ranks) * ky(my)        &
                    !!!!  * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                    !!!!                           + omg(iz)*mu(im) - 1.5_DP ) )  

                    !%%% kvd = kx(mx) * vdx(iz,iv,im) + ky(my) * vdy(iz,iv,im) %%%
                    !%%% kvs = ky(my) * vsy(iz,iv,im)                          %%%
                    vdx(iz,iv,im)=                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )     &
                      * kkx                                                     &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vdy(iz,iv,im)=                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )     &
                      * kky                                                     &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )  


                !  end do
                !end do
              end do

            end do   ! im loop ends


            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)         &
                              + 2._DP*kx(mx)*ky(my)*gg(1,2) &
                              + (ky(my)**2)*gg(2,2)
              end do
            end do

            baxfactor    = 1._DP

           !- for OUTPUT hst/*.mtr.* -
            metric_l( 1,iz) = zz(iz)   ! [ 1]
            metric_l( 2,iz) = theta    ! [ 2]
            metric_l( 3,iz) = omg(iz)  ! [ 3]
            metric_l( 4,iz) = domgdx   ! [ 4]
            metric_l( 5,iz) = domgdy   ! [ 5]
            metric_l( 6,iz) = domgdz   ! [ 6]
            metric_l( 7,iz) = gg(1,1)  ! [ 7]
            metric_l( 8,iz) = gg(1,2)  ! [ 8]
            metric_l( 9,iz) = gg(1,3)  ! [ 9]
            metric_l(10,iz) = gg(2,2)  ! [10]
            metric_l(11,iz) = gg(2,3)  ! [11]
            metric_l(12,iz) = gg(3,3)  ! [12]
            metric_l(13,iz) = rootg(iz)! [13]
           !-------------------------

!!!! for VMEC equilibrium !!!
!          else if( trim(equib_type) == "vmec" ) then
!
!            q_bar = q_0 
!            theta = zz(iz)
!
!            call vmecin_coeff( rad_a, R0_unit, rho2R_0, q_input, theta,  &
!                               alpha_fix, r_0, r_minor, s_hat,           &
!                               gdwss, gdwtt, gdwzz, gdwst, gdwsz, gdwtz, &
!                               gupss, guptt, gupzz, gupst, gupsz, guptz, &
!                               babs, Bs, Bth, Bzt, dBds, dBdt, dBdz,     &
!                               dBdt_mir, vmec_rootg, rootgft, rootgbz )
!
!            omg(iz)   = babs
!
!            rootg(iz) = vmec_rootg * R0_unit * R0_unit * R0_unit
!            dpara(iz) = dz * babs * rootgft * b0b00
!
!
!
!            do im = 0, nm
!
!              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )
!              mir(iz,im) = mu(im) * dBdt_mir / babs / rootgft / b0b00
!
!              do iv = 1, 2*nv
!                do my = ist_y, iend_y
!                  do mx = -nx, nx
!
!                    kvd(mx,my,iz,iv,im) =                                              &
!                      - (( vl(iv)**2 + omg(iz)*mu(im) ) / rootgbz /babs/babs/babs )    &
!                      * ((r_0/q_0) * ky(my)                                            &
!                      * ( ( (Bs/r_a) + Bzt * (q_0/r_0) * s_hat * zz(iz) ) * dBdt       &
!                         +( (Bs/r_a) * q_0 - Bth * (q_0/r_0) * s_hat * zz(iz) ) * dBdz &
!                         -( Bth + Bzt * q_0 ) * dBds / r_a )                           &
!                      + kx(mx) * ( Bzt * dBdt - Bth * dBdz ))                          &
!                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )
!                                              ! --- k*v_d term
!
!                    kvs(mx,my,iz,iv,im) =                                               &
!                      - sgn(ranks) * ky(my)                                             &
!                      * ((r_0/q_0) * (Bth + Bzt * q_0) / rootgbz / babs / babs)         &
!                      * ( R0_Ln(ranks)                                                  &
!                        + R0_Lt(ranks) * (0.5_DP*vl(iv)**2 + omg(iz)*mu(im) - 1.5_DP) ) &
!                      * tau(ranks) / Znum(ranks)
!                                              ! --- k*v_* term
!                  end do
!                end do
!              end do
!
!            end do   ! im loop ends
!
!
!            do my = ist_y, iend_y
!              do mx = -nx, nx
!                ksq(mx,my,iz)  = (r_a * kx(mx))**2  * gupss                             &
!                   + ky(my)**2  * ( (r_0/q_0)**2                                        &
!                   * ( gupzz + guptt * q_0 **2 - guptz * 2._DP * q_0 )                  &
!                   + 2._DP * s_hat * (r_0/q_0) * zz(iz) * r_a * ( gupst * q_0 - gupsz ) &
!                   + r_a * r_a * gupss * (s_hat**2) * (zz(iz)**2) )                     &
!                   + (r_a * kx(mx)) * ky(my) * 2._DP * ( (r_0/q_0)                      &
!                   * ( gupst * q_0 - gupsz ) + r_a * gupss * s_hat * zz(iz) ) 
!                                              ! --- squere of k_perp 
!              end do
!            end do
!
!            baxfactor    = b0b00                ! --- For the use in caldlt
!

!  this is new vmec-BoozXform interface  by M. Nakata & M. Nunami  (Aug. 2016)
          else if( trim(equib_type) == "vmec" ) then

            q_bar = q_0
            isw = 1
            r_major = 1._DP ! in the R0 unit

            call vmecbzx_boozx_coeff( isw,  nss,  ntheta,  nzeta,  s_input,     iz, zz(iz),  lz_l,   &  ! input 
                                       s_0,       q_0,     s_hat,    eps_r, phi_ax,          &  ! output
                                   omg(iz), rootg(iz),    domgdx,   domgdz, domgdy,          &
                                   gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),                  &
                                   gg(2,3),   gg(3,3)  )

            dpara(iz) = dz * omg(iz) * rootg(iz)

            kkx =  r_major*( -domgdy + (gg(1,3)*gg(1,2) - gg(1,1)*gg(2,3))*domgdz/omg(iz)**2 )
            kky =  r_major*(  domgdx + (gg(1,3)*gg(2,2) - gg(1,2)*gg(2,3))*domgdz/omg(iz)**2 )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              mir(iz,im) = mu(im) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                !do my = ist_y, iend_y
                !  do mx = -nx, nx

                    !!!!kvd(mx,my,iz,iv,im) =                                        &
                    !!!!   ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                    !!!!  * ( kkx*kx(mx) + kky*ky(my) )                              &
                    !!!!  * ( sgn(ranks) * tau(ranks) / Znum(ranks) )                &
                    !!!! - real(ibprime,kind=DP) * vl(iv)**2 / r_major / omg(iz)**2  & ! grad-p (beta-prime) term 
                    !!!!  * ( beta*(R0_Ln(ranks) + R0_Lt(ranks))*ky(my) )            &
                    !!!!  * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    !!!!kvs(mx,my,iz,iv,im) =                                     &
                    !!!!  - sgn(ranks) * tau(ranks) / Znum(ranks) * ky(my)        &
                    !!!!  * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                    !!!!                           + omg(iz)*mu(im) - 1.5_DP ) )  

                    !%%% kvd = kx(mx) * vdx(iz,iv,im) + ky(my) * vdy(iz,iv,im) %%%
                    !%%% kvs = ky(my) * vsy(iz,iv,im)                          %%%
                    vdx(iz,iv,im) =                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                      * kkx                                                      &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )


                    vdy(iz,iv,im) =                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                      * kky                                                      &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )                &
                     - real(ibprime,kind=DP) * vl(iv)**2 / r_major / omg(iz)**2  & ! grad-p (beta-prime) term 
                      * ( beta*(R0_Ln(ranks) + R0_Lt(ranks)) )                   &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )  


                !  end do
                !end do
              end do

            end do   ! im loop ends


            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)         &
                              + 2._DP*kx(mx)*ky(my)*gg(1,2) &
                              + (ky(my)**2)*gg(2,2)
              end do
            end do

            baxfactor    = 1._DP

           !- for OUTPUT hst/*.mtr.* -
            metric_l( 1,iz) = zz(iz)   ! [ 1]
            metric_l( 2,iz) = phi_ax   ! [ 2] Axisymetric toroidal angle
            metric_l( 3,iz) = omg(iz)  ! [ 3]
            metric_l( 4,iz) = domgdx   ! [ 4]
            metric_l( 5,iz) = domgdy   ! [ 5]
            metric_l( 6,iz) = domgdz   ! [ 6]
            metric_l( 7,iz) = gg(1,1)  ! [ 7]
            metric_l( 8,iz) = gg(1,2)  ! [ 8]
            metric_l( 9,iz) = gg(1,3)  ! [ 9]
            metric_l(10,iz) = gg(2,2)  ! [10]
            metric_l(11,iz) = gg(2,3)  ! [11]
            metric_l(12,iz) = gg(3,3)  ! [12]
            metric_l(13,iz) = rootg(iz)! [13]
           !-------------------------


          else if( trim(equib_type) == "eqdsk" ) then

            q_bar = q_0
            isw = 1
            r_major = 1._DP ! in the R0 unit

            call igs_coeff( isw,  mc_type,   nss,    ntheta,  s_input, zz(iz), lz_l, &  ! input 
                                       s_0,       q_0,     s_hat,    eps_r,  theta,       &  ! output
                                   omg(iz), rootg(iz),    domgdx,   domgdz, domgdy,       &
                                   gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),               &
                                   gg(2,3),   gg(3,3)  )

            dpara(iz) = dz * omg(iz) * rootg(iz)

            kkx =  r_major*( -domgdy + (gg(1,3)*gg(1,2) - gg(1,1)*gg(2,3))*domgdz/omg(iz)**2 )
            kky =  r_major*(  domgdx + (gg(1,3)*gg(2,2) - gg(1,2)*gg(2,3))*domgdz/omg(iz)**2 )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              mir(iz,im) = mu(im) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                !do my = ist_y, iend_y
                !  do mx = -nx, nx

                    !!!!kvd(mx,my,iz,iv,im) =                                        &
                    !!!!   ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                    !!!!  * ( kkx*kx(mx) + kky*ky(my) )                              &
                    !!!!  * ( sgn(ranks) * tau(ranks) / Znum(ranks) )                &
                    !!!! - real(ibprime,kind=DP) * vl(iv)**2 / r_major / omg(iz)**2  & ! grad-p (beta-prime) term 
                    !!!!  * ( beta*(R0_Ln(ranks) + R0_Lt(ranks))*ky(my) )            &
                    !!!!  * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    !!!!kvs(mx,my,iz,iv,im) =                                     &
                    !!!!  - sgn(ranks) * tau(ranks) / Znum(ranks) * ky(my)        &
                    !!!!  * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                    !!!!                           + omg(iz)*mu(im) - 1.5_DP ) )  

                    !%%% kvd = kx(mx) * vdx(iz,iv,im) + ky(my) * vdy(iz,iv,im) %%%
                    !%%% kvs = ky(my) * vsy(iz,iv,im)                          %%%
                    vdx(iz,iv,im) =                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                      * kkx                                                      &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vdy(iz,iv,im) =                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                      * kky                                                      &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )                &
                     - real(ibprime,kind=DP) * vl(iv)**2 / r_major / omg(iz)**2  & ! grad-p (beta-prime) term 
                      * ( beta*(R0_Ln(ranks) + R0_Lt(ranks)) )                   &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )  


                !  end do
                !end do
              end do

            end do   ! im loop ends


            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)         &
                              + 2._DP*kx(mx)*ky(my)*gg(1,2) &
                              + (ky(my)**2)*gg(2,2)
              end do
            end do

            baxfactor    = 1._DP

           !- for OUTPUT hst/*.mtr.* -
            metric_l( 1,iz) = zz(iz)   ! [ 1]
            metric_l( 2,iz) = theta    ! [ 2]
            metric_l( 3,iz) = omg(iz)  ! [ 3]
            metric_l( 4,iz) = domgdx   ! [ 4]
            metric_l( 5,iz) = domgdy   ! [ 5]
            metric_l( 6,iz) = domgdz   ! [ 6]
            metric_l( 7,iz) = gg(1,1)  ! [ 7]
            metric_l( 8,iz) = gg(1,2)  ! [ 8]
            metric_l( 9,iz) = gg(1,3)  ! [ 9]
            metric_l(10,iz) = gg(2,2)  ! [10]
            metric_l(11,iz) = gg(2,3)  ! [11]
            metric_l(12,iz) = gg(3,3)  ! [12]
            metric_l(13,iz) = rootg(iz)! [13]
           !-------------------------


          else

            write( olog, * ) " # wrong choice of the equilibrium "
            call flush(olog)
            call MPI_Finalize(ierr_mpi)
            stop

          end if


          do im = 0, nm
            do my = ist_y, iend_y
              do mx = -nx, nx
                kmo           = sqrt( 2._DP * ksq(mx,my,iz) * mu(im) / omg(iz) ) &
                               * dsqrt( tau(ranks)*Anum(ranks) ) / Znum(ranks)
                call math_j0( kmo, j0(mx,my,iz,im) )
                call math_j1( kmo, j1(mx,my,iz,im) )
                call math_j2( kmo, j2(mx,my,iz,im) )
              end do
            end do
          end do


          do my = ist_y, iend_y
            do mx = -nx, nx
              bb     = ksq(mx,my,iz) / omg(iz)**2 &
                        * tau(ranks)*Anum(ranks)/(Znum(ranks)**2)
              call math_g0( bb, g0(mx,my,iz) )
            end do
          end do


!!!! debug (Jan 2012)
!          write( olog, fmt="(1p,10e15.7)" )    &
!            zz(iz), omg(iz), mir(iz,0), dpara(iz), jcob(iz), &
!            ksq(1,2,iz), kvs(1,2,iz,1,0), kvd(1,2,iz,1,0), j0(1,2,iz,0), g0(1,2,iz)
!!!! debug (Jan 2012)


        end do   ! iz loop ends

!- OUTPUT ascii data hst/*.mtr.* -
        call MPI_gather(metric_l(1,-nz), num_omtr*2*nz, MPI_DOUBLE_PRECISION,        &
                        metric_g(1,-global_nz), num_omtr*2*nz, MPI_DOUBLE_PRECISION, &
                        0, zsp_comm_world, ierr_mpi)
        if ( rankg == 0 ) then
          do iz = -global_nz, global_nz-1
            write( omtr, fmt="(f15.8,SP,256E24.14e3)") metric_g(:,iz)
          end do
          call flush(omtr)
        end if
!---------------------------------

! --- operator settings ---


        cfsrf   = 0._DP
        cfsrf_l = 0._DP
        do iz = -nz, nz-1
!          cfsrf_l   = cfsrf_l + 1._DP / omg(iz)
          cfsrf_l   = cfsrf_l + rootg(iz)
                                            ! normalization coefficient for 
                                            ! the surface average
        end do

        call MPI_Allreduce( cfsrf_l, cfsrf, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, zsp_comm_world, ierr_mpi )


                                            ! --- debug
                                            !    write( olog, * ) " *** z, omg "
                                            !    do iz = -nz, nz-1
                                            !      write( olog, * ) zz(iz), omg(iz)
                                            !    end do
                                            !    write( olog, * )


        if ( vel_rank == 0 ) then
          do iz = -nz, nz-1
            !dvp(iz)   = vp(iz,1)
            dvp(iz)  = sqrt( 2._DP * (0.5_DP * dm**2) * omg(iz) )
          end do
        end if

        call MPI_Bcast( dvp, 2*nz, MPI_DOUBLE_PRECISION, 0, &
                        vel_comm_world, ierr_mpi )


        do my = ist_y_g, iend_y_g
          ck(my-ist_y_g)   = exp( ui * 2._DP * pi * del_c &
                                     * real( n_tht * my, kind=DP ) )
          dj(my-ist_y_g)   = - m_j * n_tht * my
                                            !  del_c = q_0*n_alp-int(q_0*n_alp)
                                            !  m_j   = 2*n_alp*q_d
        end do


        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              fmx(iz,iv,im)   = exp( - 0.5_DP * vl(iv)**2 - omg(iz) * mu(im) ) &
                              / sqrt( twopi**3 )
            end do
          end do
        end do

        allocate( ww(-nx:nx,0:ny,-nz:nz-1) )

! --- GK polarization factor for efield calculation 
        fct_poisson(:,:,:) = 0._DP
        fct_e_energy(:,:,:) = 0._DP

        ww(:,:,:) = 0._DP
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx

              if ( rankw == 0 .and. mx == 0 .and. my == 0 ) then !- (0,0) mode

                fct_poisson(mx,my,iz) = 0._DP
                fct_e_energy(mx,my,iz) = 0._DP

              else

                ww(mx,my,iz) = lambda_i * ksq(mx,my,iz)
                do is = 0, ns-1
                  bb   = ksq(mx,my,iz) / omg(iz)**2 &
                          * tau(is)*Anum(is)/(Znum(is)**2)
                  call math_g0( bb, gg0 )
                  ww(mx,my,iz) = ww(mx,my,iz)  &
                               + Znum(is) * fcs(is) / tau(is) * ( 1._DP - gg0 )
                end do
                fct_poisson(mx,my,iz) = 1._DP / ww(mx,my,iz)
                fct_e_energy(mx,my,iz) = ww(mx,my,iz)

              end if

            end do
          end do
        end do


! --- ZF-factor for adiabatic model
        if ( ns == 1 ) then

          ww(:,:,:) = 0._DP
          do iz = -nz, nz-1
            my = 0
              do mx = -nx, nx
                ww(mx,my,iz) = ( 1._DP - g0(mx,my,iz) )       &
                             / ( 1._DP - g0(mx,my,iz) + tau(0)*tau_ad )
              end do
          end do

          call intgrl_fsrf ( ww, fctgt )

          if ( rankw == 0 )  then
            fctgt(0)   = ( 1._DP - g0(0,0,0) ) / ( 1._DP - g0(0,0,0) + tau(0)*tau_ad )
                                              ! g0(0,0,iz) has no z dependence
          endif

        endif

        deallocate( ww )

        allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
        allocate( nw(-nx:nx,0:ny,-nz:nz-1) )
        wf(:,:,:,:,:) = ( 0._DP, 0._DP )
        nw(:,:,:) = ( 0._DP, 0._DP )

! --- GK polarization factor for mfield calculation 
        fct_ampere(:,:,:) = 0._DP
        fct_m_energy(:,:,:) = 0._DP

        if ( beta .ne. 0._DP ) then
       
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx, nx
                    wf(mx,my,iz,iv,im) = Znum(ranks) * fcs(ranks) / Anum(ranks)  &
                                       * vl(iv)**2 * j0(mx,my,iz,im)**2 * fmx(iz,iv,im)
                  end do
                end do
              end do
            end do
          end do
  
          call intgrl_v0_moment_ms ( wf, nw )
  
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                fct_ampere(mx,my,iz) = 1._DP / real( ksq(mx,my,iz) + beta * nw(mx,my,iz), kind=DP )
                fct_m_energy(mx,my,iz) = ksq(mx,my,iz) / beta
              end do
            end do
          end do
  
          if ( rankw == 0 ) then
            do iz = -nz, nz-1
              fct_ampere(0,0,iz) = 0._DP
              fct_m_energy(0,0,iz) = 0._DP
            end do
          end if

        end if

        deallocate( wf )
        deallocate( nw )


! --- set collision frequencies and v-space functions for multi-species GK collision
        read(inml,nml=nu_ref)
        call colli_set_param(q_0, eps_r, nust)
        !if (trim(time_advnc) == "imp_colli" .or. trim(time_advnc) == "auto_init") then
        if (trim(col_type) == "full" .or. trim(col_type) == "lorentz" .or. trim(time_advnc) == "imp_colli") then
          call colliimp_set_param
        end if
        !!! call colliimp_set_param
        call dtc_init( lx, ly, vmax )

        gg_g(1,1,:) = metric_l( 7,:) !%%% For shearflow_lagrange
        gg_g(1,2,:) = metric_l( 8,:) !%%% For shearflow_lagrange
        gg_g(1,3,:) = metric_l( 9,:) !%%% For shearflow_lagrange
        gg_g(2,2,:) = metric_l(10,:) !%%% For shearflow_lagrange
        gg_g(2,3,:) = metric_l(11,:) !%%% For shearflow_lagrange
        gg_g(3,3,:) = metric_l(12,:) !%%% For shearflow_lagrange
        gg_g(2,1,:) = gg_g(1,2,:)    !%%% For shearflow_lagrange
        gg_g(3,1,:) = gg_g(1,3,:)    !%%% For shearflow_lagrange
        gg_g(3,2,:) = gg_g(2,3,:)    !%%% For shearflow_lagrange
        s_hat_g = s_hat              !%%% For shearflow_lagrange->running
        eps_r_g = eps_r              !%%% For shearflow_lagrange->running

        write( olog, * ) " # Collision parameters"
        write( olog, * ) ""
        write( olog, * ) " # Nref [m^-3]  = ", Nref
        write( olog, * ) " # Lref [m]     = ", Lref
        write( olog, * ) " # Tref [keV]   = ", Tref
        write( olog, * ) " # col_type     = ", col_type
        write( olog, * ) " # iFLR         = ", iFLR
        write( olog, * ) " # icheck       = ", icheck
        write( olog, * ) 
        write( olog, * ) " # Normalized collisionality: nu*"
        do is1 = 0, ns-1
        do is2 = 0, ns-1
        write( olog, * ) " # a, b, nu*_ab = ", is1, is2, nust(is1,is2)
        end do
        end do
        write( olog, * ) 
        if ( trim(col_type) == "LB" ) then
          write( olog, * ) " # col.-freq. bias factor for LB, nu = ", nu(:)
          write( olog, * ) 
        end if

        Zeff = 0._DP
        do is1 = 1, ns-1
          Zeff = Zeff + fcs(is1)*Znum(is1)
        end do
        write( olog, * ) " # Zeff         = ", Zeff
        write( olog, * ) 



  END SUBROUTINE set_cnfig


!--------------------------------------
  SUBROUTINE set_value( ff, phi, Al, hh, time )
!--------------------------------------

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh

    real(kind=DP), intent(out) :: time


    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf, dh, cf, ef
    real(kind=DP),    dimension(:),         allocatable :: rr
    character(15) :: colliflag
    integer :: input_status
    integer :: mx, my, iz, iv, im, nx_init, gmy


      ff(:,:,:,:,:) = ( 0._DP, 0._DP )
      phi(:,:,:)    = ( 0._DP, 0._DP )
      Al(:,:,:)     = ( 0._DP, 0._DP )
      hh(:,:,:,:,:) = ( 0._DP, 0._DP )


      if( inum == 1 ) then


        time     = 0._DP

        allocate( rr((2*nx+1)*(global_ny+1)) )

        if (init_random) then
          call math_random ( rr )
        else
          rr(:) = 0._DP
        end if

        if ( nx0 > nx ) then 
          nx_init = nx
        else        
          nx_init = nx0
        end if

!!!        my = 0 for R-H test
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist1_y, iend_y
            ! do my = ist_y, iend_y
                gmy = my + (ny+1)*rankw
                do mx = -nx_init, nx_init
                  ff(mx,my,iz,iv,im)   = dns1(ranks) * fmx(iz,iv,im)  &
                          * ( 1._DP + vl(iv) + zz(iz) )**2            &
                          * exp( -zz(iz)**2 / (0.2_DP*pi)**2 ) &
                          * exp( ui * twopi * rr(mx+nx+1+(2*nx+1)*gmy) )
                end do
              end do
            end do
          end do
        end do

        if ( rankw == 0 ) then
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                my = 0
                  do mx = 1, nx
                    ff(-mx,my,iz,iv,im) = conjg( ff(mx,my,iz,iv,im) )
                  end do
                ff(0,0,iz,iv,im) = ( 0._DP, 0._DP )
              end do
            end do
          end do
        end if

        deallocate( rr )


        if ( ch_res ) call set_ch_resolution ( ff, time )


      else


        allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )

        time   = - 1._DP

        do 
          !fj start 202010
          !read( unit=icnt, iostat=input_status ) time, wf
          call fileio_read_cnt( wf, time, input_status )
          !fj end 202010

          if ( input_status < 0 ) then
            write( olog, * ) &
               " # end of file of unit=30 is detected --> stop"
            call flush(olog)
            call MPI_Abort(MPI_COMM_WORLD, ierr_mpi)
            stop
          end if

          if ( input_status > 0 ) then
            write( olog, * ) &
               " # input error of unit=30 is detected --> stop"
            call flush(olog)
            call MPI_Abort(MPI_COMM_WORLD, ierr_mpi)
            stop
          end if

          if ( time > eps ) exit
        end do


        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ff(mx,my,iz,iv,im) = wf(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        end do

        deallocate( wf )


          write( olog, * ) ""
          write( olog, * ) " # simulation is re-started at t = ", time
          write( olog, * ) ""


      end if

      if (gamma_e /= 0._DP) then           !%%% For shearflow_lagrange
        !!!call shearflow_lagrange_init(time) !%%% For shearflow_lagrange
        call shearflow_running_init(time) !%%% For shearflow_lagrange->running
      end if                               !%%% For shearflow_lagrange

      call bndry_zvm_bound_f( ff )

      call fld_esfield( ff, phi )
      if ( beta .ne. 0._DP ) then
        call fld_emfield_ff( ff, Al )
      end if
      call fld_ff2hh( ff, Al, hh )

      call tips_reality( hh )

      allocate( dh(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( cf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( ef(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )

      colliflag = "collisional"
      call caldlt_rev( colliflag, ff, phi, Al, hh, dh, cf, ef )

      deallocate( dh )
      deallocate( cf )
      deallocate( ef )

!! --- for debug
!      call MPI_Finalize(ierr_mpi)
!      stop


  END SUBROUTINE set_value


!--------------------------------------
  SUBROUTINE set_ch_resolution ( ff, time )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    real(kind=DP), intent(out) :: time


    !--- Set perpendicular resolution employed in the input file. ---
    !!! NOTE !!!
    !    Resolutions in (z,v,m,s) should be kept the same.
    !    Since lx and ly should be larger than or equal to lx0, ly0,
    !    xfold and yfold are set to be integers.
    !!!!!!!!!!!!
      integer, parameter :: xfold = 1   ! kxmin0 = xfold * kxmin
      integer, parameter :: yfold = 1   ! kymin0 = yfold * kymin
      integer, parameter :: nx0 = 95, global_ny0 = 1, nprocw0 = 1
      integer, parameter :: ny0 = global_ny0 / nprocw0
      real(kind=DP) :: amplify = 1._DP
    !----------------------------------------------------------------


    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    integer :: input_status
    integer :: ny_size0, nwk0, irw, ist_y_g0, iend_y_g0
    integer :: mx, my, iz, iv, im, mxw, myw
    character(6)   :: crank
    character(1)   :: srank



      allocate( wf(-nx0:nx0,0:ny0,-nz:nz-1,1:2*nv,0:nm) )

      do irw = 0, nprocw0-1

        ny_size0 = global_ny0 + 1 
        if( mod(ny_size0,nprocw0) == 0 )  then
          nwk0    = ny_size0 / nprocw0
        else
          nwk0    = ny_size0 / nprocw0 + 1
        endif
        !--- global index range ---------------- 
        ist_y_g0  = nwk0*irw
        iend_y_g0 = min( nwk0*(irw+1)-1, (ny_size0-1) )

        write( crank, fmt="(i6.6)" ) irw + nprocw0*rankz + nprocw0*nprocz*rankv  &
            + nprocw0*nprocz*nprocv*rankm + nprocw0*nprocz*nprocv*nprocm*ranks
        write( srank, fmt="(i1.1)" ) ranks

        !fj start 202010
        !open( icnt, file=trim(f_cnt)//crank//".cnt.000", &
        !      form="unformatted", status="old", action="read" )
        !read( unit=icnt, iostat=input_status ) time, wf
        !rewind( icnt )
        !close( icnt )
        call fileio_open_icnt( trim(f_cnt) )
        call fileio_read_cnt( wf, time, input_status )
        call fileio_close_icnt
        !fj end 202010

        if ( ist_y_g <= iend_y_g0 * yfold .and. ist_y_g0 * yfold <= iend_y_g ) then
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y_g0, iend_y_g0
                  myw = my * yfold
                  if ( ist_y_g <= myw .and. myw <= iend_y_g ) then
                    do mx = -nx0, nx0
                      mxw = mx * xfold
                      if ( -nx <= mxw .and. mxw <= nx ) then
                        ff(mxw,myw-ist_y_g,iz,iv,im) = wf(mx,my-ist_y_g0,iz,iv,im) * amplify
                      end if
                    end do
                  end if
                end do
              end do
            end do
          end do
        end if

      end do

      deallocate( wf )

      write( olog, * ) ""
      write( olog, * ) " # simulation is re-started at t = ", time
      write( olog, * ) " # perpendicular resolutions are changed:"
      write( olog, * ) " # lx is ", xfold, "times larger" 
      write( olog, * ) " # ly is ", yfold, "times larger" 
      write( olog, * ) " # nx =", nx0, "to", nx
      write( olog, * ) " # global_ny =", global_ny0, "to", global_ny
      write( olog, * ) " # nprocw from", nprocw0, "to", nprocw
      write( olog, * ) ""


  END SUBROUTINE set_ch_resolution


END MODULE GKV_set
