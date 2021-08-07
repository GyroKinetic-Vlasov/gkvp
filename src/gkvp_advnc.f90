MODULE GKV_advnc
!-------------------------------------------------------------------------------
!
!    Calculate df/dt and time advance by Runge-Kutta-Gill method
!
!    Update history of gkvp_advnc.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - Unitialized access for padding iend_y<my is removed.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_fld,   only: fld_esfield, fld_emfield_hh, fld_hh2ff
  use GKV_exb,   only: exb_NL_term
  use GKV_colli, only: colli_LB!, colli_full
  use GKV_colliimp, only: colliimp_calc_colli_full
  use GKV_bndry, only: bndry_bound_e,  &
                       bndry_zv_buffin, bndry_zv_sendrecv, bndry_zv_buffout
  use GKV_clock, only: clock_sta, clock_end
  use GKV_zfilter, only: zfilter
  use GKV_tips,  only: tips_reality
  use GKV_shearflow, only: shearflow_lagrange_update, &
                           shearflow_lagrange_dealias, kx_lagrange, &
                           shearflow_running_update

  implicit none

  private

  integer, save :: nchunk_zv = 1, nchunk_yzv = 1, nchunk_yz = 1

  public   advnc_rkgsteps_rev, caldlt_rev


CONTAINS


!--------------------------------------
  SUBROUTINE advnc_rkgsteps_rev( colliflag, ff, phi, Al, hh )
!--------------------------------------
!     time integration of GK equation using RKG method

    character(15), intent(in) :: colliflag ! = "collisional"
                                           ! = "collisionless"
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh

    complex(kind=DP), save,  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: qh
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: dh, cf, ef
    integer :: mx, my, iz, iv, im, istep
    integer, save :: iflg
    data iflg / 0 /
!$  integer :: nthreads, omp_get_num_threads

!      if ( iflg == 0 ) then
        iflg = 1
!$OMP parallel
        do im = 0, nm
!$OMP do
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = 0, ny
                do mx = -nx, nx
                  qh(mx,my,iz,iv,im) = ( 0._DP, 0._DP )
                end do
              end do
            end do
          end do
!$OMP end do nowait
        end do
!$OMP end parallel

!$OMP parallel default(shared)
!$OMP master
!$    nthreads = omp_get_num_threads()
!$    if (nthreads > 1) then
!$      nchunk_zv = ((2*nz)*(2*nv)-1) / (nthreads-1) + 1
!$      nchunk_yzv = ((iend_y-ist_y+1)*(2*nz)*(2*nv)-1) / (nthreads-1) + 1
!$      nchunk_yz = ((iend_y-ist_y+1)*(2*nz)-1) / (nthreads-1) + 1
!$    end if
!$OMP end master
!$OMP end parallel
!      end if

      allocate( dh(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( cf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( ef(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )

      do istep = 1, 4

        if (istep == 2 .or. istep == 4) then          !%%% For shearflow_lagrange
          !!!call shearflow_lagrange_update(0.5_DP * dt) !%%% For shearflow_lagrange
          call shearflow_running_update(0.5_DP * dt) !%%% For shearflow_lagrange->running
        end if                                        !%%% For shearflow_lagrange

        call caldlt_rev( colliflag, ff, phi, Al, hh, dh, cf, ef )

                                           call clock_sta(11)
                                         ! call fapp_start("rkg",11,1)
        call rkg( hh, dh, qh, istep )
                                         ! call fapp_stop("rkg",11,1)
                                           call clock_end(11)

        call tips_reality ( hh )

        !!!call shearflow_lagrange_dealias(ff, phi, Al, hh) !%%% For shearflow_lagrange

                                           call clock_sta(12)
                                         ! call fapp_start("esfield",12,1)
        if ( beta > 0._DP ) then
          call fld_emfield_hh( hh, Al )
        end if
        call fld_hh2ff( hh, Al, ff )
        call fld_esfield( ff, phi )
                                         ! call fapp_stop("esfield",12,1)
                                           call clock_end(12)
      end do

      deallocate( dh )
      deallocate( cf )
      deallocate( ef )


  END SUBROUTINE advnc_rkgsteps_rev


!--------------------------------------
  SUBROUTINE rkg( hh, dh, qh, istep )
!--------------------------------------
!     Runge-Kutta-Gill

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh, qh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh
    integer, intent(in) :: istep

    real(kind=DP) :: c1, c2, cq, c0
    integer :: mx, my, iz, iv, im


      if      ( istep == 1 ) then
        c1   =  0.5_DP
        c2   = -1._DP
        cq   = -2._DP
        c0   =  1._DP
      else if ( istep == 2 ) then
        c1   =  1._DP - sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 3 ) then
        c1   =  1._DP + sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 4 ) then
        c1   =  1._DP / 6._DP
        c2   = -1._DP / 3._DP
        cq   =  0._DP
        c0   =  0._DP
      end if

!$OMP parallel do collapse(3)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                hh(mx,my,iz,iv,im) = hh(mx,my,iz,iv,im)           &
                                   + c1 * dt * dh(mx,my,iz,iv,im) &
                                   + c2 * qh(mx,my,iz,iv,im)
                qh(mx,my,iz,iv,im) = cq * qh(mx,my,iz,iv,im) &
                                   + c0 * dt * dh(mx,my,iz,iv,im)
              end do
            end do
          end do
        end do
      end do


  END SUBROUTINE rkg


!--------------------------------------
  SUBROUTINE caldlt_rev( colliflag, ff, phi, Al, hh, dh, cf, ef )
!--------------------------------------
!     increment of delta-f within a time step

    character(15), intent(in) :: colliflag ! = "collisional"
                                           ! = "collisionless"
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: dh, cf, ef

    complex(kind=DP), dimension(:,:,:,:), allocatable :: psi, chi
    integer :: mx, my, iz, iv, im

      allocate( psi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) )
      allocate( chi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) )

!$OMP parallel default(none) &
!$OMP shared(psi,chi,phi,Al,j0,ist_y,iend_y) &
!$OMP private(mx,my,iz,im)
!$OMP do collapse(2)
      do im = 0, nm
        do iz = -nz, nz-1
          do my = iend_y, ny
            psi(:,my,iz,im) = (0._DP, 0._DP)
            chi(:,my,iz,im) = (0._DP, 0._DP)
          end do
          do my = ist_y, iend_y
            do mx = -nx, nx
              psi(mx,my,iz,im) = j0(mx,my,iz,im) * phi(mx,my,iz)
              chi(mx,my,iz,im) = j0(mx,my,iz,im) * Al(mx,my,iz)
            end do
          end do
        end do
      end do
!$OMP end do
!$OMP end parallel

                                           call clock_sta(13)
                                         ! call fapp_start("literm",13,1)
     !%%% Linear collisionless term %%%
      call caldlt_linear( ff, psi, chi, dh )

     !%%% Collision term %%%
      if ( trim(colliflag) == "collisional" ) then

        if ( trim(col_type) == "LB" ) then
          call colli_LB( ff, phi, cf )
        else if ( trim(col_type) == "full" .or. &
                  trim(col_type) == "lorentz" ) then
          !call colli_full( ff, phi, cf )
          call colliimp_calc_colli_full( ff, phi, cf )
        else 
          write(olog,*) "## Illegal choice for col_type!! ---> stop"
          call flush(olog)
          call MPI_Finalize(ierr_mpi)
          stop
        end if

      else if ( trim(colliflag) == "collisionless" ) then

!!$OMP parallel workshare
!        cf(:,:,:,:,:) = (0._DP, 0._DP)
!!$OMP end parallel workshare

      else 

        write(olog,*) "## Illegal choice for colliflag:", colliflag
        call flush(olog)
        call MPI_Finalize(ierr_mpi)
        stop

      end if
                                         ! call fapp_stop("literm",13,1)
                                           call clock_end(13)

                                           call clock_sta(14)
                                         ! call fapp_start("nlterm",14,1)
     !%%% Nonlinear term %%%
      call exb_NL_term( hh, psi, chi, ef )
                                         ! call fapp_stop("nlterm",14,1)
                                           call clock_end(14)

     !%%% dh/dt = (Linear collisionless) + (Collision) - (Nonlinear) %%%
      if ( trim(colliflag) == "collisional" ) then
!$OMP parallel
        do im = 0, nm
!$OMP do
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              !do my = 0, ny
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dh(mx,my,iz,iv,im) = dh(mx,my,iz,iv,im) &
                                     + cf(mx,my,iz,iv,im) &
                                     - ef(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
!$OMP end do nowait
        end do
!$OMP end parallel
      else if ( trim(colliflag) == "collisionless" ) then
!$OMP parallel
        do im = 0, nm
!$OMP do
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              !do my = 0, ny
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dh(mx,my,iz,iv,im) = dh(mx,my,iz,iv,im) &
                                     - ef(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
!$OMP end do nowait
        end do
!$OMP end parallel
      end if
                                           call clock_sta(15)
                                         ! call fapp_start("zfilter",15,1)
      if ( trim(z_filt) == "on" ) then
        call zfilter( dh )
      end if                                
                                         ! call fapp_stop("zfilter",15,1)
                                           call clock_end(15)

      deallocate( psi )
      deallocate( chi )

  END SUBROUTINE caldlt_rev


!--------------------------------------
  SUBROUTINE caldlt_linear( ff, psi, chi, dh )
!--------------------------------------
!     increment of delta-f within a time step

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh

    complex(kind=DP), dimension(:,:,:,:), allocatable :: &
                         zb1be, zb1te, zb2be, zb2te, zb1bo, zb1to, zb2bo, zb2to
    complex(kind=DP), dimension(:,:,:,:), allocatable :: vb1e, vb2e, vb1o, vb2o
    integer :: im


      allocate( zb1be(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
      allocate( zb1te(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
      allocate( zb2be(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
      allocate( zb2te(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
      allocate( zb1bo(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
      allocate( zb1to(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
      allocate( zb2bo(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
      allocate( zb2to(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
      allocate( vb1e(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) )
      allocate( vb2e(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) )
      allocate( vb1o(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) )
      allocate( vb2o(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) )

      call bndry_bound_e( psi )

      call literm_k_rev( ff, psi, chi, dh )

!$OMP parallel default(none) &
!$OMP shared(ff,psi,chi,dh) &
!$OMP shared(zb1be,zb1te,zb2be,zb2te,zb1bo,zb1to,zb2bo,zb2to,vb1e,vb2e,vb1o,vb2o) &
!$OMP private(im)

!!%%% Without overlap %%%
!      do im = 0, nm
!        call bndry_zv_buffin( ff(:,:,:,:,im), zb1be, zb1te, vb1e )
!!$OMP barrier
!!$OMP master
!        call bndry_zv_sendrecv( zb1be, zb1te, zb2be, zb2te, vb1e, vb2e )
!!$OMP end master
!!$OMP barrier
!        call bndry_zv_buffout( zb2be, zb2te, vb2e, ff(:,:,:,:,im) )
!!$OMP barrier
!        call literm_zv( ff(:,:,:,:,im), psi(:,:,:,im), im, dh(:,:,:,:,im) )
!!$OMP barrier
!      end do
!!%%%%%%%%%%%%%%%%%%%%%%%


!!%%% With overlap %%%
      do im = 0, nm+3
        if (mod(im,2) == 0) then ! even
!$OMP master
          if (0+1<=im .and. im<=nm+1) call bndry_zv_sendrecv( zb1bo, zb1to, zb2bo, zb2to, vb1o, vb2o )
!$OMP end master
          if (0  <=im .and. im<=nm  ) call bndry_zv_buffin( ff(:,:,:,:,im), zb1be, zb1te, vb1e )
          if (0+2<=im .and. im<=nm+2) call bndry_zv_buffout( zb2be, zb2te, vb2e, ff(:,:,:,:,im-2) )
          if (0+3<=im .and. im<=nm+3) call literm_zv( ff(:,:,:,:,im-3), psi(:,:,:,im-3), im-3, dh(:,:,:,:,im-3) )
        else                     ! odd
!$OMP master
          if (0+1<=im .and. im<=nm+1) call bndry_zv_sendrecv( zb1be, zb1te, zb2be, zb2te, vb1e, vb2e )
!$OMP end master
          if (0  <=im .and. im<=nm  ) call bndry_zv_buffin( ff(:,:,:,:,im), zb1bo, zb1to, vb1o )
          if (0+2<=im .and. im<=nm+2) call bndry_zv_buffout( zb2bo, zb2to, vb2o, ff(:,:,:,:,im-2) )
          if (0+3<=im .and. im<=nm+3) call literm_zv( ff(:,:,:,:,im-3), psi(:,:,:,im-3), im-3, dh(:,:,:,:,im-3) )
        end if
!$OMP barrier
      end do
!!%%%%%%%%%%%%%%%%%%%%

!$OMP end parallel

      deallocate( zb1be )
      deallocate( zb1te )
      deallocate( zb2be )
      deallocate( zb2te )
      deallocate( zb1bo )
      deallocate( zb1to )
      deallocate( zb2bo )
      deallocate( zb2to )
      deallocate( vb1e )
      deallocate( vb2e )
      deallocate( vb1o )
      deallocate( vb2o )


  END SUBROUTINE caldlt_linear


!--------------------------------------
  SUBROUTINE literm_k_rev ( ff, psi, chi, lf )
!--------------------------------------
!     z-derivative of ff

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: lf

    real(kind=DP) :: cs1, cs2, kvd, kvs
    integer  ::  mx, my, iz, iv, im


                                           call clock_sta(1320)
                                         ! call fapp_start("literm_perp",1320,1)

      cs1    = sgn(ranks) * Znum(ranks) / tau(ranks)
      cs2    = sqrt( tau(ranks) / Anum(ranks) )

!$OMP parallel do collapse(3) private(kvd,kvs)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                !kvd = kx(mx) * vdx(iz,iv,im) + ky(my) * vdy(iz,iv,im)
                kvd = kx_lagrange(mx,my) * vdx(iz,iv,im) + ky(my) * vdy(iz,iv,im) !%%% For shearflow_lagrange
                kvs = ky(my) * vsy(iz,iv,im)
                lf(mx,my,iz,iv,im) =                  &
                   - ui * kvd * ff(mx,my,iz,iv,im)    &
                   - cs1 * fmx(iz,iv,im) * (          &
                       + ui * kvd * psi(mx,my,iz,im)  &
                       - ui * kvs                     &
                            * ( psi(mx,my,iz,im) - cs2 * vl(iv) * chi(mx,my,iz,im) ) )
              end do
            end do
          end do
        end do
      end do
                                         ! call fapp_stop("literm_perp",1320,1)
                                           call clock_end(1320)


  END SUBROUTINE literm_k_rev


!--------------------------------------
  SUBROUTINE literm_zv ( ff, psi, im, lf )
!--------------------------------------
!     (z,v)-derivative of ff

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb)            :: psi
    integer, intent(in) :: im
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)             :: lf

    real(kind=DP), dimension(-nz:nz-1) :: cefz, cefz2
    real(kind=DP) :: cefv, cs1
    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1330)
                                         ! call fapp_start("literm_para",1330,1)
!$OMP end master

      cs1    = sgn(ranks) * Znum(ranks) / tau(ranks)
      do iz = -nz, nz-1
        cefz(iz)   = 1._DP / ( 12._DP * dpara(iz) ) * sqrt( tau(ranks) / Anum(ranks) )
        cefz2(iz)  = 1._DP / ( 60._DP * dpara(iz) ) * sqrt( tau(ranks) / Anum(ranks) )
      end do
      cefv   = 1._DP / ( 12._DP * dv ) * sqrt( tau(ranks) / Anum(ranks) )

      if (trim(z_calc) == "cf4") then

!$OMP do collapse(2) schedule(dynamic,nchunk_zv)
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                 !!!- vl(iv) * cefz(iz) * (              &
                 - (vl(iv) * cefz(iz) - gamma_e /(s_hat_g*12._DP * (zz(0)-zz(-1))))* ( & !%%% For shearflow_lagrange->running
                     -         ff(mx,my,iz+2,iv)      &
                     + 8._DP * ff(mx,my,iz+1,iv)      &
                     - 8._DP * ff(mx,my,iz-1,iv)      &
                     +         ff(mx,my,iz-2,iv) )    &
                 + mir(iz,im) * cefv * (              &
                     -         ff(mx,my,iz,iv+2)      &
                     + 8._DP * ff(mx,my,iz,iv+1)      &
                     - 8._DP * ff(mx,my,iz,iv-1)      &
                     +         ff(mx,my,iz,iv-2) )    &
                 - cs1 * fmx(iz,iv,im) * (            &
                       vl(iv) * cefz(iz) * (          &
                         -         psi(mx,my,iz+2)    &
                         + 8._DP * psi(mx,my,iz+1)    &
                         - 8._DP * psi(mx,my,iz-1)    &
                         +         psi(mx,my,iz-2) ) )&
                 - art_diff * (                       &
                     +         ff(mx,my,iz+2,iv)      &
                     - 4._DP * ff(mx,my,iz+1,iv)      &
                     + 6._DP * ff(mx,my,iz  ,iv)      &
                     - 4._DP * ff(mx,my,iz-1,iv)      &
                     +         ff(mx,my,iz-2,iv) )
            end do
          end do
        end do
      end do
!$OMP end do nowait

      else if (trim(z_calc) == "up5") then

        do iv = 1, 2*nv
          if ( vl(iv) > 0._DP ) then
!$OMP do collapse(2) schedule(dynamic,nchunk_yz)
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                     - vl(iv) * cefz2(iz) * (             &
                         - 3._DP * ff(mx,my,iz+2,iv)      &
                         +30._DP * ff(mx,my,iz+1,iv)      &
                         +20._DP * ff(mx,my,iz  ,iv)      &
                         -60._DP * ff(mx,my,iz-1,iv)      &
                         +15._DP * ff(mx,my,iz-2,iv)      &
                         - 2._DP * ff(mx,my,iz-3,iv) )    &
                     + mir(iz,im) * cefv * (              &
                         -         ff(mx,my,iz,iv+2)      &
                         + 8._DP * ff(mx,my,iz,iv+1)      &
                         - 8._DP * ff(mx,my,iz,iv-1)      &
                         +         ff(mx,my,iz,iv-2) )    &
                     - cs1 * fmx(iz,iv,im) * (            &
                           vl(iv) * cefz(iz) * (          &
                             -         psi(mx,my,iz+2)    &
                             + 8._DP * psi(mx,my,iz+1)    &
                             - 8._DP * psi(mx,my,iz-1)    &
                             +         psi(mx,my,iz-2) ) )
                end do
              end do
            end do
!$OMP end do nowait
          else
!$OMP do collapse(2) schedule(dynamic,nchunk_yz)
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                     - vl(iv) * cefz2(iz) * (             &
                         + 2._DP * ff(mx,my,iz+3,iv)      &
                         -15._DP * ff(mx,my,iz+2,iv)      &
                         +60._DP * ff(mx,my,iz+1,iv)      &
                         -20._DP * ff(mx,my,iz  ,iv)      &
                         -30._DP * ff(mx,my,iz-1,iv)      &
                         + 3._DP * ff(mx,my,iz-2,iv) )    &
                     + mir(iz,im) * cefv * (              &
                         -         ff(mx,my,iz,iv+2)      &
                         + 8._DP * ff(mx,my,iz,iv+1)      &
                         - 8._DP * ff(mx,my,iz,iv-1)      &
                         +         ff(mx,my,iz,iv-2) )    &
                     - cs1 * fmx(iz,iv,im) * (            &
                           vl(iv) * cefz(iz) * (          &
                             -         psi(mx,my,iz+2)    &
                             + 8._DP * psi(mx,my,iz+1)    &
                             - 8._DP * psi(mx,my,iz-1)    &
                             +         psi(mx,my,iz-2) ) )
                end do
              end do
            end do
!$OMP end do nowait
          end if
        end do

      else

        write(olog,*) "## Illegal choice for z_calc!! ---> stop"
        call flush(olog)
        call MPI_Finalize(ierr_mpi)
        stop

      end if

!$OMP master
                                         ! call fapp_stop("literm_para",1330,1)
                                           call clock_end(1330)
!$OMP end master


  END SUBROUTINE literm_zv


END MODULE GKV_advnc
