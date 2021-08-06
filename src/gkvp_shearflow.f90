MODULE GKV_shearflow
!-------------------------------------------------------------------------------
!
!    Shearflow convection term  
!
!    Update history of gkvp_shearflow.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!      gkvp_f0.55 (M. Nakata, Dec 2018)   
!        - First implementation
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_math,   only: math_j0, math_j1, math_j2, math_g0, math_random
  use GKV_intgrl, only: intgrl_fsrf, intgrl_v0_moment_ms
  use GKV_fld,    only: fld_esfield, fld_emfield_ff, fld_ff2hh
  use GKV_tips,  only: tips_reality
  use GKV_colliimp, only: nbuff, gmu, gj0, gj1

  implicit none

  private

  public   shearflow_kxmap, &
           shearflow_lagrange_init, shearflow_lagrange_update,    &
           shearflow_lagrange_remesh, shearflow_lagrange_dealias, &
           kx_lagrange, gkx_lagrange,                             &
           shearflow_running_init, shearflow_running_update,      &
           shearflow_running_lap, zz_lagrange
           
  !%%% - Definition
  !%%%   kx: radial wavenumber in moving frame and remesh
  !%%%   kx_lagrange: radial wavenumber in lab frame
  !%%% - j0, j1, j2, g0, ksq, fct_poisson, fct_e_energy, 
  !%%%   fct_ampere, fct_m_energy, fctgt in gkvp_header.f90 are 
  !%%%   kx-dependent geometric quantities in lab frame, updated in this module.
  real(kind=DP), dimension(-nx:nx,0:ny) :: kx_lagrange
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: gkx_lagrange

  real(kind=DP), dimension(0:global_ny) :: gky
  real(kind=DP) :: dt_remesh, t_remesh
  real(kind=DP), dimension(0:ny) :: mx_dealias_left, mx_dealias_right

  real(kind=DP), dimension(-nz:nz-1) :: zz_lagrange
  real(kind=DP) :: dt_lap, t_lap


CONTAINS


!--------------------------------------
  SUBROUTINE shearflow_kxmap( time, ff, phi, Al, hh )
!--------------------------------------
!     discrete advection in kx direction due to the mean radial flow shear

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    real(kind=DP), intent(in)           :: time

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff_tmp

    integer :: mx, my, iz, iv, im, mx_new, gmy, tloop
    integer, dimension(0:ny) :: my_map, loop_mapping

      tloop = nint(time/dt)
      my_map(:) = -1
      if (rankw == 0) loop_mapping(0) = 1

      do my = ist1_y, iend_y
 
        gmy = my + (ny+1)*rankw
        !!!loop_mapping(my) = nint(kxmin_g/(kymin_g*abs(gamma_e)*dt))/gmy
        loop_mapping(my) = nint(kxmin_g/(kymin_g*gmy*abs(gamma_e)*dt))

        if (mod(tloop+loop_mapping(my),loop_mapping(my)) == 0 ) then 
          my_map(my) = my
        else 
          my_map(my) = -1
        end if 

      end do

      if (maxval(my_map) < 0 ) then 
        return
      else 

        if ( gamma_e > 0._DP ) then

!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,my_map,ff,ff_tmp) &
!$OMP private(mx,my,iz,iv,im,mx_new)
          do im = 0-nvb, nm+nvb
            do iv = 1-nvb, 2*nv+nvb
              do iz = -nz-nzb, nz-1+nzb
                do my = ist_y, iend_y

                  if (my == my_map(my)) then       
                    do mx = -nx+1, nx
                      mx_new = mx - 1  
                      ff_tmp(mx_new,my,iz,iv,im) = ff(mx,my,iz,iv,im)
                    end do 
                    ff_tmp(nx,my,iz,iv,im) = (0._DP, 0._DP)
                  end if

                end do 
              end do 
            end do 
          end do 

        else if (gamma_e < 0._DP) then  

!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,my_map,ff,ff_tmp) &
!$OMP private(mx,my,iz,iv,im,mx_new)
          do im = 0-nvb, nm+nvb
            do iv = 1-nvb, 2*nv+nvb
              do iz = -nz-nzb, nz-1+nzb
                do my = ist_y, iend_y

                  if (my == my_map(my)) then       
                    do mx = -nx, nx-1
                      mx_new = mx + 1  
                      ff_tmp(mx_new,my,iz,iv,im) = ff(mx,my,iz,iv,im)
                    end do 
                    ff_tmp(-nx,my,iz,iv,im) = (0._DP, 0._DP)
                  end if

                end do 
              end do 
            end do 
          end do 
      
        end if 

!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,my_map,ff,ff_tmp) &
!$OMP private(mx,my,iz,iv,im,mx_new)
        do im = 0-nvb, nm+nvb
          do iv = 1-nvb, 2*nv+nvb
            do iz = -nz-nzb, nz-1+nzb
              do my = ist_y, iend_y

                if (my == my_map(my)) then       
                  do mx = -nx, nx
                    ff(mx,my,iz,iv,im) = ff_tmp(mx,my,iz,iv,im)
                  end do
                end if

              end do
            end do
          end do
        end do

        call fld_esfield( ff, phi )
        if ( beta .ne. 0._DP ) then
          call fld_emfield_ff( ff, Al )
        end if
        call fld_ff2hh( ff, Al, hh )

        call tips_reality( hh )

      end if

  END SUBROUTINE shearflow_kxmap


!--------------------------------------
  SUBROUTINE shearflow_lagrange_init(time)
!--------------------------------------
  
    real(kind=DP), intent(in) :: time
    integer :: my, gmy, num_remesh

      dt_remesh = kxmin_g / kymin_g / gamma_e
      num_remesh = int((time + 0.5_DP*dt_remesh + eps)/dt_remesh)
      t_remesh  = ( 0.5d0 + num_remesh ) * dt_remesh - eps

      do my = ist_y, iend_y
        ! kx: radial wavenumber in moving frame & remesh
        ! kx_lagrange: radial wavenumber in lab frame
        kx_lagrange(:,my) = kx(:)                   &
                          - ky(my) * gamma_e * time &
                          + ky(my) * kxmin_g / kymin_g * num_remesh
      end do
      call shearflow_lagrange_kxdep_geometry
  
      do my = ist_y, iend_y
        gmy = my + (ny+1)*rankw
        mx_dealias_left(my) = -nx+0.5*gmy-1                     &
                            + ky(my) * gamma_e * time / kxmin_g &
                            - ky(my) / kymin_g * num_remesh
        mx_dealias_right(my) = nx-0.5*gmy+1                      &
                             + ky(my) * gamma_e * time / kxmin_g &
                             - ky(my) / kymin_g * num_remesh
      end do

!--- For exb_NL_term
      do my = 0, global_ny
        gky(my) = kymin_g * real(my, kind=DP)
        gkx_lagrange(:,my) = kx(:)                    &
                           - gky(my) * gamma_e * time &
                           + gky(my) * kxmin_g / kymin_g * num_remesh
      end do


  END SUBROUTINE shearflow_lagrange_init
  
  
!--------------------------------------
  SUBROUTINE shearflow_lagrange_update(dt_shear)
!--------------------------------------
  
    real(kind=DP), intent(in) :: dt_shear
    integer :: my
  
      do my = ist_y, iend_y
        kx_lagrange(:,my) = kx_lagrange(:,my) - ky(my) * gamma_e * dt_shear
      end do
      call shearflow_lagrange_kxdep_geometry
  
      do my = ist_y, iend_y
        mx_dealias_left(my) = mx_dealias_left(my) + ky(my) * gamma_e * dt_shear / kxmin_g
        mx_dealias_right(my) = mx_dealias_right(my) + ky(my) * gamma_e * dt_shear / kxmin_g
      end do
  
!--- For exb_NL_term
      do my = 0, global_ny
        gkx_lagrange(:,my) = gkx_lagrange(:,my) - gky(my) * gamma_e * dt_shear
      end do

  END SUBROUTINE shearflow_lagrange_update


!--------------------------------------
  SUBROUTINE shearflow_lagrange_kxdep_geometry
!--------------------------------------
! - j0, j1, j2, g0, ksq, fct_poisson, fct_e_energy, 
!   fct_ampere, fct_m_energy, fctgt in gkvp_header.f90 are 
!   kx-dependent geometric quantities in lab frame, updated in this module.

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: nw
    real(kind=DP), dimension(:,:,:), allocatable :: ww
    real(kind=DP) :: gg0, bb, kmo
    integer :: mx, my, iz, iv, im, is
    integer :: ibuff, mxy
 
!$OMP parallel do default(none) &
!$OMP shared(ist_y,iend_y,ksq,kx_lagrange,gg_g,ky) &
!$OMP private(mx,my,iz)
      do iz = -nz, nz-1 
        do my = ist_y, iend_y
          do mx = -nx, nx
            ksq(mx,my,iz) = (kx_lagrange(mx,my)**2)*gg_g(1,1,iz)         &
                          + 2._DP*kx_lagrange(mx,my)*ky(my)*gg_g(1,2,iz) &
                          + (ky(my)**2)*gg_g(2,2,iz)
          end do
        end do
      end do
  
!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,ranks,ksq,mu,omg,tau,Anum,Znum,j0,j1,j2) &
!$OMP private(mx,my,iz,im,kmo)
      do im = 0, nm
        do iz = -nz, nz-1 
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
      end do
  
!$OMP parallel do default(none) &
!$OMP shared(ist_y,iend_y,ranks,ksq,omg,tau,Anum,Znum,g0) &
!$OMP private(mx,my,iz,bb)
      do iz = -nz, nz-1 
        do my = ist_y, iend_y
          do mx = -nx, nx
            bb     = ksq(mx,my,iz) / omg(iz)**2 &
                      * tau(ranks)*Anum(ranks)/(Znum(ranks)**2)
            call math_g0( bb, g0(mx,my,iz) )
          end do
        end do
      end do

! --- GK polarization factor for efield calculation 
      allocate( ww(-nx:nx,0:ny,-nz:nz-1) )
!$OMP parallel default(none) &
!$OMP shared(ist_y,iend_y,rankw,lambda_i,ksq,omg,tau,Anum,Znum,fcs,ww,fct_poisson,fct_e_energy) &
!$OMP private(mx,my,iz,is,bb,gg0)
!$OMP workshare
      ww(:,:,:) = 0._DP
      fct_poisson(:,:,:) = 0._DP
      fct_e_energy(:,:,:) = 0._DP
!$OMP end workshare
!$OMP do
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
!$OMP end do nowait
!$OMP end parallel

! --- ZF-factor for adiabatic model
      if ( ns == 1 ) then
!$OMP parallel workshare
        ww(:,:,:) = 0._DP
!$OMP end parallel workshare
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

! --- GK polarization factor for mfield calculation 
      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( nw(-nx:nx,0:ny,-nz:nz-1) )
!$OMP parallel workshare
      wf(:,:,:,:,:) = ( 0._DP, 0._DP )
      nw(:,:,:) = ( 0._DP, 0._DP )
      fct_ampere(:,:,:) = 0._DP
      fct_m_energy(:,:,:) = 0._DP
!$OMP end parallel workshare
      if ( beta .ne. 0._DP ) then
!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,ranks,Znum,fcs,Anum,vl,j0,fmx,wf) &
!$OMP private(mx,my,iz,iv,im)
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
!$OMP parallel do default(none) &
!$OMP shared(ist_y,iend_y,ksq,beta,nw,fct_ampere,fct_m_energy) &
!$OMP private(mx,my,iz)
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

!--- Bessel functions for colliimp module ---
!$OMP parallel do default(none) &
!$OMP shared(ist_y,iend_y,spc_rank,ksq,gmu,omg,tau,Anum,Znum,gj0,gj1) &
!$OMP private(mx,my,iz,is,im,mxy,ibuff,kmo)
      do ibuff = 0, nbuff-1
        mxy = ibuff + nbuff * spc_rank
        if (mxy <= (2*nx+1)*(ny+1)-1) then
          mx = mod(mxy,2*nx+1) - nx
          my = mxy / (2*nx+1)
          do iz = -nz, nz-1
            do is = 0, ns-1
              do im = 0, global_nm
                kmo = sqrt( 2._DP * ksq(mx,my,iz) * gmu(im) / omg(iz) ) &
                    * sqrt( tau(is)*Anum(is) ) / Znum(is)
                call math_j0( kmo, gj0(im,is,iz,ibuff) )
                call math_j1( kmo, gj1(im,is,iz,ibuff) )
              end do
            end do
          end do
        else
          gj0(:,:,:,ibuff) = 0._DP
          gj1(:,:,:,ibuff) = 0._DP
        end if
      end do

  END SUBROUTINE shearflow_lagrange_kxdep_geometry

  
!--------------------------------------
  SUBROUTINE shearflow_lagrange_remesh(time,ff,phi,Al,hh)
!--------------------------------------
  
    real(kind=DP), intent(in) :: time
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    integer :: mx, my, gmy, iz, iv, im
  
      if (time > t_remesh) then
        if (rankg==0) write(*,*) "remesh,t=",time,t_remesh

        do my = ist_y, iend_y
          gmy = my + (ny+1)*rankw
          if (gmy == 0) then
            !ff(-nx:nx,my,:,:,:) = ff(-nx:nx,my,:,:,:)
            !hh(-nx:nx,my,:,:,:) = hh(-nx:nx,my,:,:,:)
            !phi(-nx:nx,my,:) = phi(-nx:nx,my,:)
            !Al(-nx:nx,my,:) = Al(-nx:nx,my,:)
          else if (gmy <= 2*nx) then
            if (gamma_e > 0._DP) then
!$OMP parallel default(none) &
!$OMP shared(my,gmy,ff,hh) &
!$OMP private(iz,iv,im)
!$OMP do collapse(2)
              do im = 0-nvb, nm+nvb
                do iv = 1-nvb, 2*nv+nvb
                  do iz = -nz-nzb, nz-1+nzb
                    ff(-nx:nx-gmy,my,iz,iv,im) = ff(-nx+gmy:nx,my,iz,iv,im)
                    ff(nx-gmy+1:nx,my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP do collapse(2)
              do im = 0, nm
                do iv = 1, 2*nv
                  do iz = -nz, nz-1
                    hh(-nx:nx-gmy,my,iz,iv,im) = hh(-nx+gmy:nx,my,iz,iv,im)
                    hh(nx-gmy+1:nx,my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP end parallel
              phi(-nx:nx-gmy,my,:) = phi(-nx+gmy:nx,my,:)
              phi(nx-gmy+1:nx,my,:) = (0._DP, 0._DP)
              Al(-nx:nx-gmy,my,:) = Al(-nx+gmy:nx,my,:)
              Al(nx-gmy+1:nx,my,:) = (0._DP, 0._DP)
            else
!$OMP parallel default(none) &
!$OMP shared(my,gmy,ff,hh) &
!$OMP private(iz,iv,im)
!$OMP do collapse(2)
              do im = 0-nvb, nm+nvb
                do iv = 1-nvb, 2*nv+nvb
                  do iz = -nz-nzb, nz-1+nzb
                    ff(-nx+gmy:nx,my,iz,iv,im) = ff(-nx:nx-gmy,my,iz,iv,im)
                    ff(-nx:-nx+gmy-1,my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP do collapse(2)
              do im = 0, nm
                do iv = 1, 2*nv
                  do iz = -nz, nz-1
                    hh(-nx+gmy:nx,my,iz,iv,im) = hh(-nx:nx-gmy,my,iz,iv,im)
                    hh(-nx:-nx+gmy-1,my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP end parallel
              phi(-nx+gmy:nx,my,:) = phi(-nx:nx-gmy,my,:)
              phi(-nx:-nx+gmy-1,my,:) = (0._DP, 0._DP)
              Al(-nx+gmy:nx,my,:) = Al(-nx:nx-gmy,my,:)
              Al(-nx:-nx+gmy-1,my,:) = (0._DP, 0._DP)
            end if
          else
!$OMP parallel default(none) &
!$OMP shared(my,gmy,ff,hh) &
!$OMP private(iz,iv,im)
!$OMP do collapse(2)
              do im = 0-nvb, nm+nvb
                do iv = 1-nvb, 2*nv+nvb
                  do iz = -nz-nzb, nz-1+nzb
                    ff(:,my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP do collapse(2)
              do im = 0, nm
                do iv = 1, 2*nv
                  do iz = -nz, nz-1
                    hh(:,my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP end parallel
            phi(:,my,:) = (0._DP, 0._DP)
            Al(:,my,:) = (0._DP, 0._DP)
          end if
        end do
  
        do my = ist_y, iend_y
          kx_lagrange(:,my) = kx_lagrange(:,my) + ky(my) * kxmin_g / kymin_g
        end do
        call shearflow_lagrange_kxdep_geometry
  
        do my = ist_y, iend_y
          mx_dealias_left(my) = mx_dealias_left(my) - ky(my) / kymin_g
          mx_dealias_right(my) = mx_dealias_right(my) - ky(my) / kymin_g
        end do

!--- For exb_NL_term
        do my = 0, global_ny
          gkx_lagrange(:,my) = gkx_lagrange(:,my) + gky(my) * kxmin_g / kymin_g
        end do

        t_remesh = t_remesh + dt_remesh
      
      end if
  
  END SUBROUTINE shearflow_lagrange_remesh
  
  
!--------------------------------------
  SUBROUTINE shearflow_lagrange_dealias(ff,phi,Al,hh)
!--------------------------------------
  
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    integer :: my, gmy, iz, iv, im, mxl(0:ny), mxr(0:ny)
  
      mxl(:) = int(mx_dealias_left(:)) - 1
      mxr(:) = int(mx_dealias_right(:)) + 1
      do my = ist_y, iend_y
        gmy = my + (ny+1)*rankw
        !write(*,*) mxl(my), mx_dealias_left(my), mxr(my), mx_dealias_right(my)
        if (gmy <= nx) then
          if (-nx <= mxl(my)) then
!$OMP parallel default(none) &
!$OMP shared(my,mxl,ff,hh) &
!$OMP private(iz,iv,im)
!$OMP do collapse(2)
              do im = 0-nvb, nm+nvb
                do iv = 1-nvb, 2*nv+nvb
                  do iz = -nz-nzb, nz-1+nzb
                    ff(-nx:mxl(my),my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP do collapse(2)
              do im = 0, nm
                do iv = 1, 2*nv
                  do iz = -nz, nz-1
                    hh(-nx:mxl(my),my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP end parallel
            phi(-nx:mxl(my),my,:) = (0._DP, 0._DP)
            Al(-nx:mxl(my),my,:) = (0._DP, 0._DP)
          end if
          if (mxr(my) <= nx) then
!$OMP parallel default(none) &
!$OMP shared(my,mxr,ff,hh) &
!$OMP private(iz,iv,im)
!$OMP do collapse(2)
              do im = 0-nvb, nm+nvb
                do iv = 1-nvb, 2*nv+nvb
                  do iz = -nz-nzb, nz-1+nzb
                    ff(mxr(my):nx,my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP do collapse(2)
              do im = 0, nm
                do iv = 1, 2*nv
                  do iz = -nz, nz-1
                    hh(mxr(my):nx,my,iz,iv,im) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP end parallel
            phi(mxr(my):nx,my,:) = (0._DP, 0._DP)
            Al(mxr(my):nx,my,:) = (0._DP, 0._DP)
          end if
        else
!$OMP parallel default(none) &
!$OMP shared(my,ff,hh) &
!$OMP private(iz,iv,im)
!$OMP do collapse(2)
              do im = 0-nvb, nm+nvb
                do iv = 1-nvb, 2*nv+nvb
                  do iz = -nz-nzb, nz-1+nzb
                    ff(:,my,:,:,:) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP do collapse(2)
              do im = 0, nm
                do iv = 1, 2*nv
                  do iz = -nz, nz-1
                    hh(:,my,:,:,:) = (0._DP, 0._DP)
                  end do
                end do
              end do
!$OMP end do nowait
!$OMP end parallel
            phi(:,my,:) = (0._DP, 0._DP)
            Al(:,my,:) = (0._DP, 0._DP)
        end if
      end do
  
  END SUBROUTINE shearflow_lagrange_dealias


!--------------------------------------
  SUBROUTINE shearflow_running_init(time)
!--------------------------------------
  
    real(kind=DP), intent(in) :: time
    integer :: my, num_lap

      dt_lap = 2._DP * pi * s_hat_g / gamma_e
      num_lap = int((time + 0.5_DP*dt_lap + eps)/dt_lap)
      t_lap  = ( 0.5d0 + num_lap ) * dt_lap - eps

      do my = ist_y, iend_y
        ! kx: radial wavenumber in moving frame & lap
        ! kx_lagrange: radial wavenumber in lab frame
        kx_lagrange(:,my) = kx(:)                   &
                          - ky(my) * gamma_e * time &
                          + sign(2._DP * pi * num_lap * s_hat_g * ky(my), gamma_e*s_hat_g)
      end do
      zz_lagrange(:) = zz(:)                    &
                     + gamma_e / s_hat_g * time &
                     - sign(2._DP * pi * num_lap, gamma_e*s_hat_g)
      call shearflow_running_kxzzdep_geometry
  
!--- For exb_NL_term
      do my = 0, global_ny
        gky(my) = kymin_g * real(my, kind=DP)
        gkx_lagrange(:,my) = kx(:)                    &
                           - gky(my) * gamma_e * time &
                           + sign(2._DP * pi * num_lap * s_hat_g * gky(my), gamma_e*s_hat_g)
      end do


  END SUBROUTINE shearflow_running_init
  
  
!--------------------------------------
  SUBROUTINE shearflow_running_update(dt_shear)
!--------------------------------------
  
    real(kind=DP), intent(in) :: dt_shear
    integer :: my
  
      do my = ist_y, iend_y
        kx_lagrange(:,my) = kx_lagrange(:,my) - ky(my) * gamma_e * dt_shear
      end do
      zz_lagrange(:) = zz_lagrange(:) + gamma_e / s_hat_g * dt_shear
      call shearflow_running_kxzzdep_geometry
  
!--- For exb_NL_term
      do my = 0, global_ny
        gkx_lagrange(:,my) = gkx_lagrange(:,my) - gky(my) * gamma_e * dt_shear
      end do

  END SUBROUTINE shearflow_running_update


!--------------------------------------
  SUBROUTINE shearflow_running_kxzzdep_geometry
!--------------------------------------
! - j0, j1, j2, g0, ksq, fct_poisson, fct_e_energy, 
!   fct_ampere, fct_m_energy, fctgt in gkvp_header.f90 are 
!   kx-dependent geometric quantities in lab frame, updated in this module.

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: nw
    real(kind=DP), dimension(:,:,:), allocatable :: ww
    real(kind=DP) :: gg0, bb, kmo, cfsrf_l
    integer :: mx, my, iz, iv, im, is
    integer :: ibuff, mxy
    real(kind=DP) :: q_bar, r_major, alpha_MHD, p_total, dp_totaldx, beta_total
    real(kind=DP) :: theta, domgdz, domgdx, domgdy, kkx, kky, dz, dm

    !%%% Valid only for s-alpha %%%
      dz = zz(0)-zz(-1)
      do iz = -nz, nz-1
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
        theta = zz_lagrange(iz)
        omg(iz)   = 1._DP - eps_r_g * cos( theta )        ! s-alpha with eps-expansion
        rootg(iz) = q_0*r_major/omg(iz)
        dpara(iz) = dz* q_0 * r_major
        domgdz  = eps_r_g * sin( theta )
        domgdx  = -cos( theta ) / r_major
        domgdy  = 0._DP
        gg_g(1,1,iz) = 1._DP
        gg_g(1,2,iz) = s_hat_g*zz_lagrange(iz) - alpha_MHD*sin(zz_lagrange(iz)) ! with Shafranov shift
        gg_g(1,3,iz) = 0._DP
        gg_g(2,1,iz) = gg_g(1,2,iz)
        gg_g(2,2,iz) = 1._DP + (s_hat_g*zz_lagrange(iz) - alpha_MHD*sin(zz_lagrange(iz)))**2 ! with Shafranov shift
        gg_g(2,3,iz) = 1._DP/(r_major*eps_r_g)
        gg_g(3,1,iz) = gg_g(1,3,iz)
        gg_g(3,2,iz) = gg_g(2,3,iz)
        gg_g(3,3,iz) = 1._DP/((r_major*eps_r_g)**2)
        kkx = -r_major * (q_0/q_bar) &
                       * ( gg_g(1,1,iz)*gg_g(2,3,iz) - gg_g(1,2,iz)*gg_g(1,3,iz) )*domgdz
        kky =  r_major * (q_bar/q_0) &
                       * ( domgdx - ( gg_g(1,2,iz)*gg_g(2,3,iz) - gg_g(2,2,iz)*gg_g(1,3,iz) )*domgdz )
        do im = 0, nm
          vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )
          mir(iz,im) = mu(im) * (q_0/q_bar) * domgdz / ( omg(iz)*rootg(iz) )
          do iv = 1, 2*nv
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
          end do
        end do   ! im loop ends
      end do
    !%%%%

      cfsrf   = 0._DP
      cfsrf_l = 0._DP
      do iz = -nz, nz-1
        cfsrf_l   = cfsrf_l + rootg(iz)
      end do
      call MPI_Allreduce( cfsrf_l, cfsrf, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, zsp_comm_world, ierr_mpi )

      dm = sqrt(2._DP*mu(1)) - sqrt(2._DP*mu(0))
      if ( vel_rank == 0 ) then
        do iz = -nz, nz-1
          dvp(iz)  = sqrt( 2._DP * (0.5_DP * dm**2) * omg(iz) )
        end do
      end if
      call MPI_Bcast( dvp, 2*nz, MPI_DOUBLE_PRECISION, 0, &
                      vel_comm_world, ierr_mpi )
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            fmx(iz,iv,im)   = exp( - 0.5_DP * vl(iv)**2 - omg(iz) * mu(im) ) &
                            / sqrt( twopi**3 )
          end do
        end do
      end do
 
 
!$OMP parallel do default(none) &
!$OMP shared(ist_y,iend_y,ksq,kx_lagrange,gg_g,ky) &
!$OMP private(mx,my,iz)
      do iz = -nz, nz-1 
        do my = ist_y, iend_y
          do mx = -nx, nx
            ksq(mx,my,iz) = (kx_lagrange(mx,my)**2)*gg_g(1,1,iz)         &
                          + 2._DP*kx_lagrange(mx,my)*ky(my)*gg_g(1,2,iz) &
                          + (ky(my)**2)*gg_g(2,2,iz)
          end do
        end do
      end do
  
!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,ranks,ksq,mu,omg,tau,Anum,Znum,j0,j1,j2) &
!$OMP private(mx,my,iz,im,kmo)
      do im = 0, nm
        do iz = -nz, nz-1 
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
      end do
  
!$OMP parallel do default(none) &
!$OMP shared(ist_y,iend_y,ranks,ksq,omg,tau,Anum,Znum,g0) &
!$OMP private(mx,my,iz,bb)
      do iz = -nz, nz-1 
        do my = ist_y, iend_y
          do mx = -nx, nx
            bb     = ksq(mx,my,iz) / omg(iz)**2 &
                      * tau(ranks)*Anum(ranks)/(Znum(ranks)**2)
            call math_g0( bb, g0(mx,my,iz) )
          end do
        end do
      end do

! --- GK polarization factor for efield calculation 
      allocate( ww(-nx:nx,0:ny,-nz:nz-1) )
!$OMP parallel default(none) &
!$OMP shared(ist_y,iend_y,rankw,lambda_i,ksq,omg,tau,Anum,Znum,fcs,ww,fct_poisson,fct_e_energy) &
!$OMP private(mx,my,iz,is,bb,gg0)
!$OMP workshare
      ww(:,:,:) = 0._DP
      fct_poisson(:,:,:) = 0._DP
      fct_e_energy(:,:,:) = 0._DP
!$OMP end workshare
!$OMP do
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
!$OMP end do nowait
!$OMP end parallel

! --- ZF-factor for adiabatic model
      if ( ns == 1 ) then
!$OMP parallel workshare
        ww(:,:,:) = 0._DP
!$OMP end parallel workshare
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

! --- GK polarization factor for mfield calculation 
      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( nw(-nx:nx,0:ny,-nz:nz-1) )
!$OMP parallel workshare
      wf(:,:,:,:,:) = ( 0._DP, 0._DP )
      nw(:,:,:) = ( 0._DP, 0._DP )
      fct_ampere(:,:,:) = 0._DP
      fct_m_energy(:,:,:) = 0._DP
!$OMP end parallel workshare
      if ( beta .ne. 0._DP ) then
!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,ranks,Znum,fcs,Anum,vl,j0,fmx,wf) &
!$OMP private(mx,my,iz,iv,im)
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
!$OMP parallel do default(none) &
!$OMP shared(ist_y,iend_y,ksq,beta,nw,fct_ampere,fct_m_energy) &
!$OMP private(mx,my,iz)
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

!--- Bessel functions for colliimp module ---
!$OMP parallel do default(none) &
!$OMP shared(ist_y,iend_y,spc_rank,ksq,gmu,omg,tau,Anum,Znum,gj0,gj1) &
!$OMP private(mx,my,iz,is,im,mxy,ibuff,kmo)
      do ibuff = 0, nbuff-1
        mxy = ibuff + nbuff * spc_rank
        if (mxy <= (2*nx+1)*(ny+1)-1) then
          mx = mod(mxy,2*nx+1) - nx
          my = mxy / (2*nx+1)
          do iz = -nz, nz-1
            do is = 0, ns-1
              do im = 0, global_nm
                kmo = sqrt( 2._DP * ksq(mx,my,iz) * gmu(im) / omg(iz) ) &
                    * sqrt( tau(is)*Anum(is) ) / Znum(is)
                call math_j0( kmo, gj0(im,is,iz,ibuff) )
                call math_j1( kmo, gj1(im,is,iz,ibuff) )
              end do
            end do
          end do
        else
          gj0(:,:,:,ibuff) = 0._DP
          gj1(:,:,:,ibuff) = 0._DP
        end if
      end do

  END SUBROUTINE shearflow_running_kxzzdep_geometry

  
!--------------------------------------
  SUBROUTINE shearflow_running_lap(time)
!--------------------------------------
  
    real(kind=DP), intent(in) :: time
    integer :: my
  
      if (time > t_lap) then
        if (rankg==0) write(*,*) "lap,t=",time,t_lap

        do my = ist_y, iend_y
          kx_lagrange(:,my) = kx_lagrange(:,my) &
                            + sign(2._DP * pi * s_hat_g * ky(my), gamma_e*s_hat_g)
        end do
        zz_lagrange(:) = zz_lagrange(:) &
                       - sign(2._DP * pi, gamma_e*s_hat_g)
        !!!call shearflow_running_kxzzdep_geometry
  
!--- For exb_NL_term
        do my = 0, global_ny
          gkx_lagrange(:,my) = gkx_lagrange(:,my) &
                             + sign(2._DP * pi * s_hat_g * gky(my), gamma_e*s_hat_g)
        end do

        t_lap = t_lap + dt_lap
      
      end if
  
  END SUBROUTINE shearflow_running_lap
  
  
END MODULE GKV_shearflow
