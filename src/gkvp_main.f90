PROGRAM GKV_main
!-------------------------------------------------------------------------------
!
!    GKV+: nonlinear gyrokinetic Vlasov code in a flux tube geometry
!
!    Hierarchy of the modules (The lower should be complied earlier)
!    ------------------------
!        main
!         |
!        set, out
!         |
!        advnc, dtc, trans
!         |
!        colli, colliimp, exb, shearflow
!         |
!        bndry, fft, fld, zfilter
!         |
!        clock, intgrl, tips, freq, igs, vmecbzx, fileio
!         |
!        mpienv, math
!         |
!        header
!
!    Update history of gkvp_main.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - Adapt to modification of freq module.
!      gkvp_f0.52 (S. Maeyama, Sep 2018)
!        - Updated for implicit collision solver.
!      gkvp_f0.40 (M. Nakata, June 2014)
!        - Updated for realistic tokamak equilibrium, 
!          multi-species collision 
!      gkvp_f0.30 (S. Maeyama, March 2013)
!        - Updated for electromagnetic, multi-species,
!          MHD equilibrium, 5D-parallelization
!      gkvp_r0.3 (T.-H. Watanabe, Jun 2011)
!        - GKV is rearranged to Fortran90 module style. 
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_set,   only: set_init, set_close
  use GKV_clock, only: clock_timer, clock_sta, clock_end, clock_reset
  use GKV_out,   only: out_cntrl, out_contnu
  use GKV_dtc,   only: dtc_cntrl, flag_time_advnc, flag_time_split
  use GKV_fld,   only: fld_esfield
  use GKV_advnc, only: advnc_rkgsteps_rev
  use GKV_colliimp, only: colliimp_colli
  use GKV_fft,   only: fft_pre
  use GKV_freq,  only: freq_set, freq_conv
  use GKV_tips,  only: tips_flush
  use GKV_shearflow,  only: shearflow_kxmap, shearflow_lagrange_remesh, &
                            shearflow_running_lap

  implicit none

  complex(kind=DP), &
    dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

  complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1)       :: Al, phi

  complex(kind=DP), &
    dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh

  real(kind=DP) :: time
  character(15) :: colliflag
  integer :: loop, iflg, cflg


    call mpienv_init( nprocw, nprocz, nprocv, nprocm, nprocs )

    call clock_timer( 0, iflg )
                                           call clock_sta(1)
                                         ! call fapp_start("pre",1,1)
    call fft_pre( )
    call set_init( ff, phi, Al, hh, time )
      write( olog, * ) " # simulation is started at t = ", time

    if ( calc_type == "lin_freq" ) call freq_set( time )
    call out_cntrl( ff, phi, Al, hh, time, 0 )

    if ( adapt_dt ) call dtc_cntrl( ff, phi, Al, hh, time, 0 )
                                         ! call fapp_stop("pre",1,1)
                                           call clock_end(1)
                                           call clock_reset
    
    loop   = 0
    cflg   = 0
    call flush(olog)
                                           call clock_sta(2)
                                         !!call PAT_region_begin(1,"timesteploop",ierr_mpi)
                                         ! call fipp_start
                                         ! call fapp_start("timesteploop",2,1)

    do

      if ( time > tend - eps ) exit

      time   = time + dt
      loop   = loop + 1

      if (flag_time_advnc == 0) then ! 4th-order RKG explicit time integration

        colliflag = "collisional"
        call advnc_rkgsteps_rev( colliflag, ff, phi, Al, hh )

      else if (flag_time_advnc == 1) then ! 2nd-order operator split with implicit collision

        if (flag_time_split == 0) then
                                           call clock_sta(17)
                                         ! call fapp_start("colliimp",17,1)
          call colliimp_colli( 0.5_DP*dt, ff, phi, Al, hh )
                                         ! call fapp_stop("colliimp",17,1)
                                           call clock_end(17)
          colliflag = "collisionless"
          call advnc_rkgsteps_rev( colliflag, ff, phi, Al, hh )
          flag_time_split = 1
                                           !%%% For debug %%%
                                           ! if (rankg==0) write(olog,*) &
                                           ! loop, time, "half-step"
                                           !%%%%%%%%%%%%%%%%%
        else if (flag_time_split == 1) then
                                           call clock_sta(17)
                                         ! call fapp_start("colliimp",17,1)
          call colliimp_colli( dt, ff, phi, Al, hh )
                                         ! call fapp_stop("colliimp",17,1)
                                           call clock_end(17)
          colliflag = "collisionless"
          call advnc_rkgsteps_rev( colliflag, ff, phi, Al, hh )
                                           !%%% For debug %%%
                                           ! if (rankg==0) write(olog,*) &
                                           ! loop, time, "one-step"
                                           !%%%%%%%%%%%%%%%%%
        end if

      end if

      if (gamma_e /= 0._DP) then
        !call shearflow_kxmap( time, ff, phi, Al, hh )
        !if (time > tlim_exb - eps .AND. cflg == 0 ) then 
        !  write( olog, * ) ""
        !  write( olog, * ) " ########## CAUTION! ############"
        !  write( olog, * ) " # time variable exceeds the time-limit: tlim_exb = ", tlim_exb
        !  write( olog, * ) " # --> GKV is still running, but you need to check the results after tlim_exb."
        !  write( olog, * ) " ########## CAUTION! ############"
        !  write( olog, * ) ""
        !  cflg = 1
        !  !!! exit 
        !end if
        !!!call shearflow_lagrange_remesh(time,ff,phi,Al,hh) !%%% For shearflow_lagrange
        call shearflow_running_lap(time)                     !%%% For shearflow_lagrange->running
      end if

                                           call clock_sta(10)
                                         ! call fapp_start("output",10,1)
      call out_cntrl( ff, phi, Al, hh, time, 1 )
      if ( adapt_dt ) call dtc_cntrl( ff, phi, Al, hh, time, 1 )
      if ( calc_type == "lin_freq" ) then
        if ( all(freq_conv) ) then
          write( olog, * ) " # Growth rate and frequency are well converged."
          exit
        end if
      end if
                                         ! call fapp_stop("output",10,1)
                                           call clock_end(10)

! --- output continu file every 10000 steps
      if (mod(loop+10000,10000) == 0 ) then 
                                           call clock_sta(16)
                                         ! call fapp_start("checkp",16,1)
        write( olog, * ) "# check-point at time = ", time
        call out_contnu ( ff, time )
        call tips_flush
                                         ! call fapp_stop("checkp",16,1)
                                           call clock_end(16)
      end if
! ---
      call clock_timer( 1, iflg )
      
      if( iflg == 1 ) exit

    end do
                                         ! call fapp_stop("timesteploop",2,1)
                                         ! call fipp_stop
                                         !!call PAT_region_end(1,ierr_mpi)
                                           call clock_end(2)

                                           call clock_sta(3)
                                         ! call fapp_start("post",3,1)
    call out_cntrl( ff, phi, Al, hh, time, 2 )
      write( olog, * ) " # simulation is stopped at t = ", time
                                         ! call fapp_stop("post",3,1)
                                           call clock_end(3)
    call clock_timer( 2, iflg )

    call set_close

    call MPI_Finalize ( ierr_mpi )

  stop


END PROGRAM GKV_main
