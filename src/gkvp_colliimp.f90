MODULE GKV_colliimp
!-------------------------------------------------------------------------------
!
!    Collision term by implicit solver
!
!    Update history of gkvp_colliimp.f90
!    --------------
!      gkvp_f0.61 (S. Maeyama, Mar 2021)
!        - Treat tracer particles (fcs=0), that has no field-particle collision.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - Bias factor nu is modified.
!        - Initialization of padding iend_y<my is added.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_clock, only : clock_sta, clock_end
  use GKV_math, only : math_j0, math_j1, math_j2
  use GKV_fld, only : fld_esfield, fld_emfield_hh, fld_hh2ff


  implicit none

  private

  integer, parameter :: nprocvms = nprocv * nprocm * nprocs
  integer, parameter :: nbuff = ((2*nx+1)*(ny+1) - 1)/nprocvms + 1
                        !%%% NOTE %%%
                        ! if (mod((2*nx+1)*(ny+1),nprocvms)==0) then
                        !   nbuff = ((2*nx+1)*(ny+1))/nprocvms 
                        ! else
                        !   nbuff = ((2*nx+1)*(ny+1))/nprocvms + 1
                        ! end if
                        !%%%%%%%%%%%%
  real(kind=DP), save :: gvl(1:2*global_nv)
  real(kind=DP), save :: gmu(0:global_nm)
  real(kind=DP), save :: gvp(0:global_nm,-nz:nz-1)
  real(kind=DP), save :: gfmx(1:2*global_nv,0:global_nm,-nz:nz-1)
  real(kind=DP), save, &
    dimension(0:global_nm,0:ns-1,-nz:nz-1,0:nbuff-1) :: gj0, gj1
  real(kind=DP), save, &
    dimension(1:2*global_nv,0:global_nm,0:ns-1,-nz:nz-1) :: gnu_ds, gnu_ps, &
                                                            gnu_hs, gnu_gs
  real(kind=DP), save, &
    dimension(1:6,1:2*global_nv,0:global_nm,0:ns-1,0:ns-1,-nz:nz-1) &
                                             :: gvfunc, gx_tst, gy_fld

  integer, parameter :: iter_max = 100
  real(kind=DP), parameter :: res_error_max = 1.d-12

  public   colliimp_set_param, colliimp_colli, colliimp_calc_colli_full, &
           gvl, gvp, gnu_ds, gnu_ps, gnu_hs, gnu_gs, &
           nbuff, gmu, gj0, gj1 !%%% For shearflow_lagrange


CONTAINS


!--------------------------------------
  SUBROUTINE colliimp_set_param
!--------------------------------------
!   Set parameters for GK collision term

    real(kind=DP), dimension(:,:,:,:,:), allocatable :: gnu_d, gnu_p, gnu_h, gnu_g
    real(kind=DP) :: dm, kmo, gxxa, cph, dph, cgg, &
                     gc_t01, gc_t02, cintgrl
    integer :: mx, my, iz, iv, im, is, mxy, ibuff, ia, ib
                                      !%%% For debug %%%
                                      ! integer :: iproc
                                      !%%%%%%%%%%%%%%%%%
                                      call clock_sta(1700)

      allocate( gnu_d(1:2*global_nv,0:global_nm,0:ns-1,0:ns-1,-nz:nz-1) )
      allocate( gnu_p(1:2*global_nv,0:global_nm,0:ns-1,0:ns-1,-nz:nz-1) )
      allocate( gnu_h(1:2*global_nv,0:global_nm,0:ns-1,0:ns-1,-nz:nz-1) )
      allocate( gnu_g(1:2*global_nv,0:global_nm,0:ns-1,0:ns-1,-nz:nz-1) )

      if ( ns == 1 ) then
        write(olog,*) "# Adiabatic model (ns==1) is not supported in imp_colli"
        call flush(olog)
        call MPI_finalize(ierr_mpi)
        stop
      end if


      do iv = 1, 2*global_nv
        gvl(iv) = dv * ( real( iv - nv * nprocv - 1, kind=DP ) + 0.5_DP )
      end do

      dm = vmax / real( nprocm * ( nm+1 ) - 1, kind=DP )
      do im = 0, global_nm
        gmu(im) = 0.5_DP * ( dm * real( im, kind=DP ) )**2
      end do

      do iz = -nz, nz-1
        do im = 0, global_nm
          gvp(im,iz)  = sqrt( 2._DP * gmu(im) * omg(iz) )
        end do
      end do

      do iz = -nz, nz-1
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            gfmx(iv,im,iz) = exp( - 0.5_DP * gvl(iv)**2 - omg(iz) * gmu(im) ) &
                           / sqrt( twopi**3 )
          end do
        end do
      end do

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

      do iz = -nz, nz-1
        do ib = 0, ns-1
          do ia = 0, ns-1
            do im = 0, global_nm
              do iv = 1, 2*global_nv

                gxxa = dsqrt(gvl(iv)**2 + gvp(im,iz)**2) / dsqrt(2._DP)
                cph = derf(calpha(ia,ib) * gxxa)
                dph = 2._DP / dsqrt(pi) * dexp(- calpha(ia,ib)**2 * gxxa**2)
                cgg = (cph - calpha(ia,ib) * gxxa * dph)/(calpha(ia,ib)**2 * gxxa**2) * 0.5_DP

                gnu_d(iv,im,ia,ib,iz) = 0.75_DP*dsqrt(pi)*ctauiv(ia,ib)*(cph-cgg)/gxxa**3
                gnu_p(iv,im,ia,ib,iz) = 1.50_DP*dsqrt(pi)*ctauiv(ia,ib)*(  cgg  )/gxxa**3
               !gnu_h(iv,im,ia,ib,iz) = 0.75_DP*dsqrt(pi)*ctauiv(ia,ib)*calpha(ia,ib)*dph/gxxa**2
                gnu_h(iv,im,ia,ib,iz) = 0.75_DP*dsqrt(pi)*ctauiv(ia,ib)*calpha(ia,ib)*dph
                gnu_g(iv,im,ia,ib,iz) = gnu_p(iv,im,ia,ib,iz)*gxxa**2*(1._DP-calpha(ia,ib)**2)

                gc_t01 = - (1._DP + calpha(ia,ib)**2)*gfmx(iv,im,iz)*gnu_p(iv,im,ia,ib,iz) * gxxa**2*gvl(iv)
                gc_t02 = - 1.5_DP*dsqrt(pi)*ctauiv(ia,ib)*gfmx(iv,im,iz) &
                       * ( cph - calpha(ia,ib)*gxxa*(1._DP + calpha(ia,ib)**2)*dph ) / calpha(ia,ib)**2 / gxxa
                cintgrl = 2._DP * pi * gvp(im,iz) * dv * dvp(iz)

                gvfunc(1,iv,im,ia,ib,iz) = gc_t01 / gfmx(iv,im,iz)                          * cintgrl
                gvfunc(2,iv,im,ia,ib,iz) = gc_t01 / gfmx(iv,im,iz) * (gvp(im,iz) / gvl(iv)) * cintgrl
                gvfunc(3,iv,im,ia,ib,iz) = gc_t02 / gfmx(iv,im,iz)                          * cintgrl
                gvfunc(4,iv,im,ia,ib,iz) = gvl(iv)                                          * cintgrl
                gvfunc(5,iv,im,ia,ib,iz) = gvp(im,iz)                                       * cintgrl
                gvfunc(6,iv,im,ia,ib,iz) = (gxxa**2 - 1.5_DP)                               * cintgrl 

                gx_tst(1,iv,im,ia,ib,iz) = (ctheta(ia,ib) - 1._DP) * gfmx(iv,im,iz) * gvl(iv)
                gx_tst(2,iv,im,ia,ib,iz) = gx_tst(1,iv,im,ia,ib,iz) * gvp(im,iz) / gvl(iv) 
                gx_tst(3,iv,im,ia,ib,iz) = gx_tst(1,iv,im,ia,ib,iz) * (gxxa**2/1.5_DP - 1._DP) / gvl(iv)
                gx_tst(4,iv,im,ia,ib,iz) = (ctheta(ia,ib) - 1._DP) * (gc_t01                              &
                                            - (ctheta(ia,ib) - 1._DP) * calpha(ia,ib) * ctauiv(ia,ib)     &
                                              * gfmx(iv,im,iz) * gvl(iv) / dsqrt(1._DP + calpha(ia,ib)**2) )
                gx_tst(5,iv,im,ia,ib,iz) =  gx_tst(4,iv,im,ia,ib,iz) * gvp(im,iz) / gvl(iv)  
                gx_tst(6,iv,im,ia,ib,iz) = (ctheta(ia,ib) - 1._DP) * (gc_t02 * 2._DP/3._DP            &
                                            - (ctheta(ia,ib) - 1._DP) * calpha(ia,ib) * ctauiv(ia,ib) &
                                              * gfmx(iv,im,iz) * (gxxa**2/1.5_DP - 1._DP) * 2._DP     &
                                                / (1._DP + calpha(ia,ib)**2)**1.5 )

                if (fcs(ia) == 0.d0 .or. fcs(ib) == 0.d0) then !-care for tracer particle(dens=0)-
                  gy_fld(1:6,iv,im,ia,ib,iz) = 0._DP
                else
                  gy_fld(1,iv,im,ia,ib,iz) = - (fcs(ib)/Znum(ib)) / (fcs(ia)/Znum(ia)) * calpha(ia,ib) * Anum(ia)  & 
                                               * tau(ib) * ctheta(ia,ib) * ctheta(ib,ia) / tau(ia) / cgamma(ia,ib) &
                                                        * ( gc_t01 - cxi(ia,ib) * gfmx(iv,im,iz) * gvl(iv) ) 
                  gy_fld(2,iv,im,ia,ib,iz) = gy_fld(1,iv,im,ia,ib,iz) * gvp(im,iz) / gvl(iv)  
                  gy_fld(3,iv,im,ia,ib,iz) = - (fcs(ib)/Znum(ib)) / (fcs(ia)/Znum(ia))                          & 
                                                        * tau(ib) * ctheta(ia,ib) * ctheta(ib,ia) / ceta(ia,ib) &
                                                        * ( gc_t02                                          &
                                                            - cxi(ia,ib)/(1._DP+calpha(ia,ib)**2)*gfmx(iv,im,iz)&
                                                              *(2._DP*gxxa**2 - 3._DP) ) 
                  gy_fld(4,iv,im,ia,ib,iz) = - gy_fld(1,iv,im,ia,ib,iz)*cxi(ib,ia) 
                  gy_fld(5,iv,im,ia,ib,iz) = - gy_fld(2,iv,im,ia,ib,iz)*cxi(ib,ia) 
                  gy_fld(6,iv,im,ia,ib,iz) = - gy_fld(3,iv,im,ia,ib,iz)*2._DP*cxi(ib,ia)/(1._DP+calpha(ib,ia)**2) 
                end if

              end do
            end do
          end do
        end do
      end do

      do iz = -nz, nz-1
        do ia = 0, ns-1
          gnu_ds(:,:,ia,iz) = 0._DP
          gnu_ps(:,:,ia,iz) = 0._DP
          gnu_hs(:,:,ia,iz) = 0._DP
          gnu_gs(:,:,ia,iz) = 0._DP
          do ib = 0, ns-1
            gnu_ds(:,:,ia,iz) = gnu_ds(:,:,ia,iz) + nu(ia) * gnu_d(:,:,ia,ib,iz)
            gnu_ps(:,:,ia,iz) = gnu_ps(:,:,ia,iz) + nu(ia) * gnu_p(:,:,ia,ib,iz)
            gnu_hs(:,:,ia,iz) = gnu_hs(:,:,ia,iz) + nu(ia) * gnu_h(:,:,ia,ib,iz)
            gnu_gs(:,:,ia,iz) = gnu_gs(:,:,ia,iz) + nu(ia) * gnu_g(:,:,ia,ib,iz)
          end do              !- nu(ia) is a bias factor given in namelist
        end do
      end do
                                      !%%% For debug %%%
                                      ! mx = 1; my = 1; iz = 0; is = 0
                                      ! if (rankw == 0 .and. rankz == 0 .and. &
                                      !     ranks == is) then
                                      !   do im = 0, nm
                                      !   do iv = 1, 2*nv
                                      !     write(910000+rankg,*)  &
                                      !         vl(iv), vp(iz,im), &
                                      !         fmx(iz,iv,im), &
                                      !         j0(mx,my,iz,im), &
                                      !         nu_ds(iz,iv,im)
                                      !   end do
                                      !   write(910000+rankg,*)
                                      !   end do
                                      ! end if
                                      !
                                      ! mx = 1; my = 1; iz = 0; is = 0
                                      ! if (rankw == 0 .and. rankz == 0) then
                                      ! call mxmy2ibuffiproc(mx,my,ibuff,iproc)
                                      ! if (spc_rank == iproc) then
                                      !   do im = 0, global_nm
                                      !   do iv = 1, 2*global_nv
                                      !     write(91,*) gvl(iv), gvp(im,iz), &
                                      !                 gfmx(iv,im,iz),      &
                                      !                 gj0(im,is,iz,ibuff), &
                                      !                 gnu_ds(iv,im,is,iz)
                                      !   end do
                                      !   write(91,*)
                                      !   end do
                                      ! end if
                                      ! end if
                                      !%%%%%%%%%%%%%%%%%
      deallocate( gnu_d )
      deallocate( gnu_p )
      deallocate( gnu_h )
      deallocate( gnu_g )
                                      call clock_end(1700)


  END SUBROUTINE colliimp_set_param


  SUBROUTINE ibuffiproc2mxmy(ibuff, iproc, mx, my)
    integer, intent(in) :: ibuff, iproc
    integer, intent(out) :: mx, my
    integer :: mxy
      mxy = ibuff + nbuff * iproc
      !if (mxy <= (2*nx+1)*(ny+1)-1) then
        mx = mod(mxy,2*nx+1) - nx
        my = mxy / (2*nx+1)
      !end if
  END SUBROUTINE ibuffiproc2mxmy


  SUBROUTINE mxmy2ibuffiproc(mx, my, ibuff, iproc)
    integer, intent(in) :: mx, my
    integer, intent(out) :: ibuff, iproc
    integer :: mxy
      mxy = mx + nx + (2*nx+1) * my
      ibuff = mod(mxy,nbuff)
      iproc = mxy / nbuff
  END SUBROUTINE mxmy2ibuffiproc


!!!!--------------------------------------
!!!  SUBROUTINE colliimp_colli( ldt, ff, phi, Al, hh )
!!!!--------------------------------------
!!!!   Collsion operator calculation interface
!!!
!!!    real(kind=DP), intent(in) :: ldt
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1)             :: phi, Al
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
!!!
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: sender, recver
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wh
!!!    integer :: iz, ibuff
!!!
!!!                                      !%%% For debug %%%
!!!                                      ! complex(kind=DP) :: wphi, wAl
!!!                                      ! integer :: mx, my, iproc
!!!                                      !%%%%%%%%%%%%%%%%%
!!!      allocate( sender(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) )
!!!      allocate( recver(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) )
!!!      allocate( wh(1:2*global_nv,0:global_nm,0:ns-1,-nz:nz-1,0:nbuff-1) )
!!!
!!!      call vms2xy_pack(hh, sender)
!!!      call vms2xy_transpose(sender, recver)
!!!      call vms2xy_unpack(recver, wh)
!!!      
!!!                                      call clock_sta(1720)
!!!!$OMP parallel default(none) &
!!!!$OMP shared(ldt,wh) &
!!!!$OMP private(iz,ibuff)
!!!!$OMP do collapse(2)
!!!      do ibuff = 0, nbuff-1
!!!        do iz = -nz, nz-1
!!!          call implicit_collision_solver(ldt, wh(:,:,:,iz,ibuff), iz, ibuff)
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!!$OMP end parallel
!!!                                      call clock_end(1720)
!!!
!!!      call xy2vms_pack(wh, sender)
!!!      call xy2vms_transpose(sender, recver)
!!!      call xy2vms_unpack(recver, hh)
!!!
!!!      if ( beta .ne. 0._DP ) then
!!!        call fld_emfield_hh( hh, Al )
!!!      end if
!!!      call fld_hh2ff( hh, Al, ff )
!!!      call fld_esfield( ff, phi )
!!!                                      !%%% For debug %%%
!!!                                      ! mx = 1; my = 1
!!!                                      ! call mxmy2ibuffiproc(mx,my,ibuff,iproc)
!!!                                      ! if (rankw == 0 .and. &
!!!                                      !     spc_rank == iproc) then
!!!                                      !   do iz = -nz, nz-1
!!!                                      !     call emfield_hh(wh(:,:,:,iz,ibuff),&
!!!                                      !                     wAl, iz, ibuff)
!!!                                      !     call esfield(wh(:,:,:,iz,ibuff), &
!!!                                      !                  wphi, iz, ibuff)
!!!                                      !     write(10000+rankg,*) zz(iz), &
!!!                                      !            dble(wphi),           &
!!!                                      !            aimag(wphi),          &
!!!                                      !            dble(phi(mx,my,iz)),  &
!!!                                      !            aimag(phi(mx,my,iz)), &
!!!                                      !            dble(wAl),            &
!!!                                      !            aimag(wAl),           &
!!!                                      !            dble(Al(mx,my,iz)),   &
!!!                                      !            aimag(Al(mx,my,iz))
!!!                                      !   end do
!!!                                      ! end if
!!!                                      ! call MPI_finalize(ierr_mpi)
!!!                                      ! stop
!!!                                      !%%%%%%%%%%%%%%%%%
!!!      deallocate( sender )
!!!      deallocate( recver )
!!!      deallocate( wh )
!!!
!!!  END SUBROUTINE colliimp_colli
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE vms2xy_pack(hh, sender)
!!!!--------------------------------------
!!!!   Data pack for transpose from (x,y*,z*,v*,m*,s*) decomposition
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
!!!    complex(kind=DP), intent(out), &
!!!      dimension(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) :: sender
!!!
!!!    integer :: mx, my, iz, iv, im, ibuff, iproc
!!!
!!!
!!!                                      call clock_sta(1710)
!!!!$OMP parallel default(none) &
!!!!$OMP shared(hh,sender) &
!!!!$OMP private(mx,my,iz,iv,im,ibuff,iproc)
!!!!$OMP workshare
!!!      sender(:,:,:,:,:) = (0._DP, 0._DP)
!!!!$OMP end workshare
!!!!$OMP do collapse(2)
!!!      do im = 0, nm
!!!        do iv = 1, 2*nv
!!!          do iz = -nz, nz-1
!!!            do my = 0, ny
!!!              do mx = -nx, nx
!!!                call mxmy2ibuffiproc(mx, my, ibuff, iproc)
!!!                sender(iv,im,iz,ibuff,iproc) = hh(mx,my,iz,iv,im)
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!!$OMP end parallel
!!!                                      call clock_end(1710)
!!!                                      !%%% For debug %%%
!!!                                      ! if (rankg == 0) then
!!!                                      ! do my = 0, ny
!!!                                      ! do mx = -nx, nx
!!!                                      !   mxy = mx+nx + (2*nx+1)*my
!!!                                      !   ibuff = mod(mxy, nbuff)
!!!                                      !   iproc = mxy / nbuff
!!!                                      !   write(99,*) mx, my, mxy, ibuff, iproc
!!!                                      ! end do
!!!                                      ! end do
!!!                                      ! end if
!!!                                      !%%%%%%%%%%%%%%%%%
!!!
!!!  END SUBROUTINE vms2xy_pack
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE vms2xy_transpose(sender, recver)
!!!!--------------------------------------
!!!!   Transpose from (x,y*,z*,v*,m*,s*) to ((x,y*)***,z*,v,m,s) decomposition
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) :: sender
!!!    complex(kind=DP), intent(out), &
!!!      dimension(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) :: recver
!!!
!!!
!!!                                      call clock_sta(1711)
!!!      call mpi_alltoall( sender,                      &
!!!                         (2*nv)*(nm+1)*(2*nz)*nbuff,  &
!!!                         MPI_DOUBLE_COMPLEX,          &
!!!                         recver,                      &
!!!                         (2*nv)*(nm+1)*(2*nz)*nbuff,  &
!!!                         MPI_DOUBLE_COMPLEX,          &
!!!                         spc_comm_world,              &
!!!                         ierr_mpi     )
!!!                                      call clock_end(1711)
!!!
!!!  END SUBROUTINE vms2xy_transpose
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE vms2xy_unpack(recver, wh)
!!!!--------------------------------------
!!!!   Data unpack for transpose to ((x,y*)***,z*,v,m,s) decomposition
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) :: recver
!!!    complex(kind=DP), intent(out), &
!!!      dimension(1:2*global_nv,0:global_nm,0:ns-1,-nz:nz-1,0:nbuff-1) :: wh
!!!
!!!    integer :: iz, iv, im, giv, gim, ibuff, irs, irm, irv, iproc
!!!
!!!
!!!                                      call clock_sta(1712)
!!!!$OMP parallel default(none) &
!!!!$OMP shared(recver,wh) &
!!!!$OMP private(iz,iv,im,ibuff,irs,irm,irv,iproc,gim,giv)
!!!!$OMP workshare
!!!      wh(:,:,:,:,:) = (0._DP, 0._DP)
!!!!$OMP end workshare
!!!!$OMP do collapse(5)
!!!      do ibuff = 0, nbuff-1
!!!        do iz = -nz, nz-1
!!!          do irs = 0, nprocs-1
!!!            do irm = 0, nprocm-1
!!!              do irv = 0, nprocv-1
!!!                iproc = irv + nprocv*irm + nprocv*nprocm*irs
!!!                do im = 0, nm
!!!                  gim = im + (nm+1)*irm
!!!                  do iv = 1, 2*nv
!!!                    giv = iv + (2*nv)*irv
!!!                    wh(giv,gim,irs,iz,ibuff) = recver(iv,im,iz,ibuff,iproc)
!!!                  end do
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!!$OMP end parallel
!!!                                      call clock_end(1712)
!!!                                      !%%% For debug %%%
!!!                                      ! if (rankg == 0) then
!!!                                      ! do irs = 0, nprocs-1
!!!                                      ! do irm = 0, nprocm-1
!!!                                      ! do irv = 0, nprocv-1
!!!                                      !   iproc = irv + nprocv*irm &
!!!                                      !         + nprocv*nprocm*irs
!!!                                      ! do im = 0, nm
!!!                                      !   gim = im + (nm+1)*irm
!!!                                      ! do iv = 1, 2*nv
!!!                                      !   giv = iv + (2*nv)*irv
!!!                                      !   write(97,*) irs, irm, irv, iproc, &
!!!                                      !               giv, gim
!!!                                      ! end do
!!!                                      ! end do
!!!                                      ! end do
!!!                                      ! end do
!!!                                      ! end do
!!!                                      ! end if
!!!                                      !%%%%%%%%%%%%%%%%%
!!!
!!!  END SUBROUTINE vms2xy_unpack
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE xy2vms_pack(wh, sender)
!!!!--------------------------------------
!!!!   Data pack for transpose from ((x,y*)***,z*,v,m,s) decomposition
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(1:2*global_nv,0:global_nm,0:ns-1,-nz:nz-1,0:nbuff-1) :: wh
!!!    complex(kind=DP), intent(out), &
!!!      dimension(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) :: sender
!!!
!!!    integer :: iz, iv, im, giv, gim, ibuff, irs, irm, irv, iproc
!!!
!!!
!!!                                      call clock_sta(1730)
!!!!$OMP parallel default(none) &
!!!!$OMP shared(wh,sender) &
!!!!$OMP private(iz,iv,im,ibuff,irs,irm,irv,iproc,gim,giv)
!!!!$OMP workshare
!!!      sender(:,:,:,:,:) = (0._DP, 0._DP)
!!!!$OMP end workshare
!!!!$OMP do collapse(5)
!!!      do ibuff = 0, nbuff-1
!!!        do iz = -nz, nz-1
!!!          do irs = 0, nprocs-1
!!!            do irm = 0, nprocm-1
!!!              do irv = 0, nprocv-1
!!!                iproc = irv + nprocv*irm + nprocv*nprocm*irs
!!!                do im = 0, nm
!!!                  gim = im + (nm+1)*irm
!!!                  do iv = 1, 2*nv
!!!                    giv = iv + (2*nv)*irv
!!!                    sender(iv,im,iz,ibuff,iproc) = wh(giv,gim,irs,iz,ibuff)
!!!                  end do
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!!$OMP end parallel
!!!                                      call clock_end(1730)
!!! 
!!!  END SUBROUTINE xy2vms_pack
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE xy2vms_transpose(sender, recver)
!!!!--------------------------------------
!!!!   Transpose from ((x,y*)***,z*,v,m,s) to (x,y*,z*,v*,m*,s*) decomposition
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(1:2*nv,0:nm,0:nbuff-1,-nz:nz-1,0:nprocvms-1) :: sender
!!!    complex(kind=DP), intent(out), &
!!!      dimension(1:2*nv,0:nm,0:nbuff-1,-nz:nz-1,0:nprocvms-1) :: recver
!!!
!!!
!!!                                      call clock_sta(1731)
!!!      call mpi_alltoall( sender,                      &
!!!                         (2*nv)*(nm+1)*(2*nz)*nbuff,  &
!!!                         MPI_DOUBLE_COMPLEX,          &
!!!                         recver,                      &
!!!                         (2*nv)*(nm+1)*(2*nz)*nbuff,  &
!!!                         MPI_DOUBLE_COMPLEX,          &
!!!                         spc_comm_world,              &
!!!                         ierr_mpi     )
!!!                                      call clock_end(1731)
!!!
!!!  END SUBROUTINE xy2vms_transpose
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE xy2vms_unpack(recver, hh)
!!!!--------------------------------------
!!!!   Data unpack for transpose to (x,y*,z*,v*,m*,s*) decomposition
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) :: recver
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
!!!
!!!    integer :: mx, my, iz, iv, im, ibuff, iproc
!!!
!!!
!!!                                      call clock_sta(1732)
!!!!$OMP parallel default(none) &
!!!!$OMP shared(recver,hh) &
!!!!$OMP private(mx,my,iz,iv,im,ibuff,iproc)
!!!!$OMP do collapse(2)
!!!      do im = 0, nm
!!!        do iv = 1, 2*nv
!!!          do iz = -nz, nz-1
!!!            do my = 0, ny
!!!              do mx = -nx, nx
!!!                call mxmy2ibuffiproc(mx, my, ibuff, iproc)
!!!                hh(mx,my,iz,iv,im) = recver(iv,im,iz,ibuff,iproc)
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!!$OMP end parallel
!!!                                      call clock_end(1732)
!!!
!!!  END SUBROUTINE xy2vms_unpack


!--------------------------------------
  SUBROUTINE colliimp_colli( ldt, ff, phi, Al, hh )
!--------------------------------------
!   Collsion operator calculation interface

    real(kind=DP), intent(in) :: ldt
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)             :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh

    complex(kind=DP), dimension(:,:,:,:), allocatable ::  &
           send1e, recv1e, send2e, recv2e, send1o, recv1o, send2o, recv2o
    complex(kind=DP), dimension(:,:,:,:), allocatable :: wh1e, wh2e, wh1o, wh2o
    integer :: iv, im, iz, is, ibuff, iproc

      allocate( send1e(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( recv1e(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( send2e(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( recv2e(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( wh1e(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) )
      allocate( wh2e(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) )
      allocate( send1o(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( recv1o(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( send2o(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( recv2o(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( wh1o(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) )
      allocate( wh2o(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) )

!$OMP parallel default(none) &
!$OMP shared(ldt,hh) &
!$OMP shared(send1e,recv1e,send2e,recv2e,wh1e,wh2e) &
!$OMP shared(send1o,recv1o,send2o,recv2o,wh1o,wh2o) &
!$OMP private(iv,im,iz,is,ibuff,iproc)
      do iproc = 0, nprocvms-1
!$OMP do
        do ibuff = 0, nbuff-1
          do im = 0, nm
            do iv = 1, 2*nv
              send1e(iv,im,ibuff,iproc) = (0._DP, 0._DP)
              send2e(iv,im,ibuff,iproc) = (0._DP, 0._DP)
              send1o(iv,im,ibuff,iproc) = (0._DP, 0._DP)
              send2o(iv,im,ibuff,iproc) = (0._DP, 0._DP)
            end do
          end do
        end do
!$OMP end do nowait
      end do
      do ibuff = 0, nbuff-1
!$OMP do
        do is = 0, ns-1
          do im = 0, global_nm
            do iv = 1, 2*global_nv
              wh1e(iv,im,is,ibuff) = (0._DP, 0._DP)
              wh1o(iv,im,is,ibuff) = (0._DP, 0._DP)
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP barrier

!!%%% Without overlap %%%
!      do iz = -nz, nz-1
!        call vms2xy_pack_iz(iz, hh, send1e)
!!$OMP barrier
!!$OMP master
!        call vms2xy_transpose_iz(send1e, recv1e)
!!$OMP end master
!!$OMP barrier
!        call vms2xy_unpack_iz(recv1e, wh1e)
!!$OMP barrier
!        call solver_wrapper_iz(ldt, iz, wh1e, wh2e)
!!$OMP barrier
!        call xy2vms_pack_iz(wh2e, send2e)
!!$OMP barrier
!!$OMP master
!        call xy2vms_transpose_iz(send2e, recv2e)
!!$OMP end master
!!$OMP barrier
!        call xy2vms_unpack_iz(iz, recv2e, hh)
!!$OMP barrier
!      end do
!!%%%%%%%%%%%%%%%%%%%%%%%

!%%% With overlap %%%
      do iz = -nz, nz-1+6

        if (mod(iz+nz,2) == 0) then ! even

!$OMP master
          if (-nz+1<=iz .and. iz<=nz-1+1) call vms2xy_transpose_iz(send1o, recv1o)
          if (-nz+5<=iz .and. iz<=nz-1+5) call xy2vms_transpose_iz(send2o, recv2o)
!$OMP end master
          if (-nz  <=iz .and. iz<=nz-1  ) call vms2xy_pack_iz(iz, hh, send1e)
          if (-nz+2<=iz .and. iz<=nz-1+2) call vms2xy_unpack_iz(recv1e, wh1e)
          if (-nz+3<=iz .and. iz<=nz-1+3) call solver_wrapper_iz(ldt, iz-3, wh1o, wh2o)
          if (-nz+4<=iz .and. iz<=nz-1+4) call xy2vms_pack_iz(wh2e, send2e)
          if (-nz+6<=iz .and. iz<=nz-1+6) call xy2vms_unpack_iz(iz-6, recv2e, hh)

        else                        ! odd

!$OMP master
          if (-nz+1<=iz .and. iz<=nz-1+1) call vms2xy_transpose_iz(send1e, recv1e)
          if (-nz+5<=iz .and. iz<=nz-1+5) call xy2vms_transpose_iz(send2e, recv2e)
!$OMP end master
          if (-nz  <=iz .and. iz<=nz-1  ) call vms2xy_pack_iz(iz, hh, send1o)
          if (-nz+2<=iz .and. iz<=nz-1+2) call vms2xy_unpack_iz(recv1o, wh1o)
          if (-nz+3<=iz .and. iz<=nz-1+3) call solver_wrapper_iz(ldt, iz-3, wh1e, wh2e)
          if (-nz+4<=iz .and. iz<=nz-1+4) call xy2vms_pack_iz(wh2o, send2o)
          if (-nz+6<=iz .and. iz<=nz-1+6) call xy2vms_unpack_iz(iz-6, recv2o, hh)

        end if
!$OMP barrier

      end do
!%%%%%%%%%%%%%%%%%%%%

!$OMP end parallel

      if ( beta .ne. 0._DP ) then
        call fld_emfield_hh( hh, Al )
      end if
      call fld_hh2ff( hh, Al, ff )
      call fld_esfield( ff, phi )

      deallocate( send1e )
      deallocate( recv1e )
      deallocate( send2e )
      deallocate( recv2e )
      deallocate( wh1e )
      deallocate( wh2e )
      deallocate( send1o )
      deallocate( recv1o )
      deallocate( send2o )
      deallocate( recv2o )
      deallocate( wh1o )
      deallocate( wh2o )

  END SUBROUTINE colliimp_colli


!--------------------------------------
  SUBROUTINE vms2xy_pack_iz(iz, hh, sender)
!--------------------------------------
!   Data pack for transpose from (x,y*,z*,v*,m*,s*) decomposition

    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(inout), &
      dimension(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) :: sender

    integer :: mx, my, iv, im, ibuff, iproc

!$OMP master
                                      call clock_sta(1710)
!$OMP end master
      do im = 0, nm
!$OMP do schedule(dynamic)
        do iv = 1, 2*nv
            do my = 0, ny
              do mx = -nx, nx
                call mxmy2ibuffiproc(mx, my, ibuff, iproc)
                sender(iv,im,ibuff,iproc) = hh(mx,my,iz,iv,im)
              end do
            end do
        end do
!$OMP end do nowait
      end do
!$OMP master
                                      call clock_end(1710)
!$OMP end master

  END SUBROUTINE vms2xy_pack_iz


!--------------------------------------
  SUBROUTINE vms2xy_transpose_iz(sender, recver)
!--------------------------------------
!   Transpose from (x,y*,z*,v*,m*,s*) to ((x,y*)***,z*,v,m,s) decomposition

    complex(kind=DP), intent(in), &
      dimension(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) :: sender
    complex(kind=DP), intent(out), &
      dimension(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) :: recver

                                      call clock_sta(1711)
      call mpi_alltoall( sender,               &
                         (2*nv)*(nm+1)*nbuff,  &
                         MPI_DOUBLE_COMPLEX,   &
                         recver,               &
                         (2*nv)*(nm+1)*nbuff,  &
                         MPI_DOUBLE_COMPLEX,   &
                         spc_comm_world,       &
                         ierr_mpi     )
                                      call clock_end(1711)

  END SUBROUTINE vms2xy_transpose_iz


!--------------------------------------
  SUBROUTINE vms2xy_unpack_iz(recver, wh)
!--------------------------------------
!   Data unpack for transpose to ((x,y*)***,z*,v,m,s) decomposition

    complex(kind=DP), intent(in), &
      dimension(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) :: recver
    complex(kind=DP), intent(inout), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) :: wh

    integer :: iv, im, giv, gim, ibuff, irs, irm, irv, iproc

!$OMP master
                                      call clock_sta(1712)
!$OMP end master
!$OMP do collapse(4) schedule(dynamic)
      do ibuff = 0, nbuff-1
          do irs = 0, nprocs-1
            do irm = 0, nprocm-1
              do irv = 0, nprocv-1
                iproc = irv + nprocv*irm + nprocv*nprocm*irs
                do im = 0, nm
                  gim = im + (nm+1)*irm
                  do iv = 1, 2*nv
                    giv = iv + (2*nv)*irv
                    wh(giv,gim,irs,ibuff) = recver(iv,im,ibuff,iproc)
                  end do
                end do
              end do
            end do
          end do
      end do
!$OMP end do nowait
!$OMP master
                                      call clock_end(1712)
!$OMP end master

  END SUBROUTINE vms2xy_unpack_iz


!--------------------------------------
  SUBROUTINE xy2vms_pack_iz(wh, sender)
!--------------------------------------
!   Data pack for transpose from ((x,y*)***,z*,v,m,s) decomposition

    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) :: wh
    complex(kind=DP), intent(inout), &
      dimension(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) :: sender

    integer :: iv, im, giv, gim, ibuff, irs, irm, irv, iproc

!$OMP master
                                      call clock_sta(1730)
!$OMP end master
!$OMP do collapse(4) schedule(dynamic)
      do ibuff = 0, nbuff-1
          do irs = 0, nprocs-1
            do irm = 0, nprocm-1
              do irv = 0, nprocv-1
                iproc = irv + nprocv*irm + nprocv*nprocm*irs
                do im = 0, nm
                  gim = im + (nm+1)*irm
                  do iv = 1, 2*nv
                    giv = iv + (2*nv)*irv
                    sender(iv,im,ibuff,iproc) = wh(giv,gim,irs,ibuff)
                  end do
                end do
              end do
            end do
          end do
      end do
!$OMP end do nowait
!$OMP master
                                      call clock_end(1730)
!$OMP end master
 
  END SUBROUTINE xy2vms_pack_iz


!--------------------------------------
  SUBROUTINE xy2vms_transpose_iz(sender, recver)
!--------------------------------------
!   Transpose from ((x,y*)***,z*,v,m,s) to (x,y*,z*,v*,m*,s*) decomposition

    complex(kind=DP), intent(in), &
      dimension(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) :: sender
    complex(kind=DP), intent(out), &
      dimension(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) :: recver

                                      call clock_sta(1731)
      call mpi_alltoall( sender,               &
                         (2*nv)*(nm+1)*nbuff,  &
                         MPI_DOUBLE_COMPLEX,   &
                         recver,               &
                         (2*nv)*(nm+1)*nbuff,  &
                         MPI_DOUBLE_COMPLEX,   &
                         spc_comm_world,       &
                         ierr_mpi     )
                                      call clock_end(1731)

  END SUBROUTINE xy2vms_transpose_iz


!--------------------------------------
  SUBROUTINE xy2vms_unpack_iz(iz, recver, hh)
!--------------------------------------
!   Data unpack for transpose to (x,y*,z*,v*,m*,s*) decomposition

    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) :: recver
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh

    integer :: mx, my, iv, im, ibuff, iproc

!$OMP master
                                      call clock_sta(1732)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic)
      do im = 0, nm
        do iv = 1, 2*nv
            do my = 0, ny
              do mx = -nx, nx
                call mxmy2ibuffiproc(mx, my, ibuff, iproc)
                hh(mx,my,iz,iv,im) = recver(iv,im,ibuff,iproc)
              end do
            end do
        end do
      end do
!$OMP end do nowait
!$OMP master
                                      call clock_end(1732)
!$OMP end master

  END SUBROUTINE xy2vms_unpack_iz


!--------------------------------------
  SUBROUTINE solver_wrapper_iz(ldt, iz, whin, whout)
!--------------------------------------
!   Implicit collision solver wrapper for overlap

    real(kind=DP), intent(in) :: ldt
    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) :: whin
    complex(kind=DP), intent(inout), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) :: whout

    complex(kind=DP), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: wc3
    integer :: ibuff

!$OMP master
                                      call clock_sta(1720)
!$OMP end master
!$OMP do schedule(dynamic)
      do ibuff = 0, nbuff-1
        wc3(:,:,:) = whin(:,:,:,ibuff)
        call implicit_collision_solver(ldt, wc3, iz, ibuff)
        whout(:,:,:,ibuff) = wc3(:,:,:)
      end do
!$OMP end do nowait
!$OMP master
                                      call clock_end(1720)
!$OMP end master

  END SUBROUTINE solver_wrapper_iz


!--------------------------------------
  SUBROUTINE implicit_collision_solver(ldt, wh, iz, ibuff)
!-------------------------------------------------------------------------------
!
!    2nd C.N. : 2nd-order Crank-Nicolson method
!
!        (h(t+dt) - h(t))/dt = C((g(t+dt) + g(t))/2)
!      where g is a function of h, i.e., g(t)=g(h(t)).
!      Since C is a linearized collision operator, it is regard as A*x = b,
!        A = 1 - dt/2*Cg,
!        x = h(t+dt),
!        b = g(h(t)) + dt/2*C(g(h(t))).
!
!    1st B.E. : 1st-order backward Euler method
!
!        (f(t+dt) - f(t))/dt = C(f(t+dt))
!      It is regard as A*x = b,
!        A = 1 - dt*Cg,
!        x = f(t+dt),
!        b = g(h(t)).
!
!    Difference between 2nd C.N and 1st B.E. appears in "calc_ap" and "calc_b".
!
!-------------------------------------------------------------------------------

    real(kind=DP), intent(in) :: ldt
    complex(kind=DP), intent(inout), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: wh
    integer, intent(in) :: iz, ibuff
  
    complex(kind=DP), dimension(1:2*global_nv,0:global_nm,0:ns-1) :: b, r
    real(kind=DP) :: res_error
    integer :: iter
  
      !if (maxval(abs(wh(:,:,:))) < 1.d-100) then
      if (maxval(abs(wh(:,:,:))) < 1.d-40) then
        return
      else
        res_error = 1._DP
      end if

      call calc_b(ldt, wh, b, iz, ibuff)
      do iter = 0, iter_max-1
        if (res_error < res_error_max) exit
        call bi_cgstab(ldt, iter, wh, b, r, iz, ibuff)
       !call gcr      (ldt, iter, wh, b, r, iz, ibuff)
       !call gcr1     (ldt, iter, wh, b, r, iz, ibuff)
       !call gcrm     (ldt, iter, wh, b, r, iz, ibuff)
        call calc_error(r, b, res_error)
      end do
      if (iter == iter_max) then
        write(olog,*) "# Reach maximum iteration : ", iz, ibuff, iter, res_error
      end if
                                      !%%% For debug %%%
                                      ! write(olog,*) "#iter= ", iter, iz, ibuff
                                      !%%%%%%%%%%%%%%%%%
  
  END SUBROUTINE implicit_collision_solver


!--------------------------------------
  SUBROUTINE bi_cgstab(ldt, iter, x, b, r, iz, ibuff)
!-------------------------------------------------------------------------------
!
!    Bi-conjugate gradient stabilization method (A*x = b for Non-Hermitian A)
!
!-------------------------------------------------------------------------------

    real(kind=DP), intent(in) :: ldt
    integer, intent(in) :: iter
    complex(kind=DP), intent(inout), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: x, r
    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: b
    integer, intent(in) :: iz, ibuff

    complex(kind=DP), save, &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: p, r0
!$OMP threadprivate(p,r0)
    complex(kind=DP), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: s, ap, as
    complex(kind=DP) :: alpha, beta, omega, sum1, sum2, sum3
  
      if (iter == 0) then                           ! r_0 = b - A*x_0, p_0 = r_0
        call calc_r(ldt, x, b, r, iz, ibuff)
        p(:,:,:) = r(:,:,:)
        r0(:,:,:) = r(:,:,:)
      end if
  
      call calc_ap(ldt, p, ap, iz, ibuff)
      sum1 = sum(r0(:,:,:) * r(:,:,:))
      sum2 = sum(r0(:,:,:) * ap(:,:,:))
      alpha = sum1 / sum2                      ! alpha = (r_0,r_k) / (r_0,A*p_k)
  
      s(:,:,:) = r(:,:,:) - alpha * ap(:,:,:)          ! s_k = r_k - alpha*A*p_k
  
      call calc_ap(ldt, s, as, iz, ibuff)
      sum2 = sum(s(:,:,:) * as(:,:,:))
      sum3 = sum(as(:,:,:)**2)
      omega = sum2 / sum3                  ! omega = (A*s_k,s_k) / (A*s_k,A*s_k)
  
      x(:,:,:) = x(:,:,:) + alpha * p(:,:,:) + omega * s(:,:,:)
                                           ! x_k+1 = x_k * alpha*p_k + omega*s_k
      r(:,:,:) = s(:,:,:) - omega * as(:,:,:)        ! r_k+1 = s_k - omega*A*s_k
  
      sum3 = sum(r0(:,:,:) * r(:,:,:))
      beta = alpha * sum3 / (omega * sum1)
                              ! beta = (alpha / omega) * (r_0,r_k+1) / (r_0,r_k)
  
      p(:,:,:) = r(:,:,:) + beta * (p(:,:,:) - omega * ap(:,:,:))
                                    ! p_k+1 = r_k+1 + beta * (p_k - omega*A*p_k)
  
  END SUBROUTINE bi_cgstab


!!--------------------------------------
!  SUBROUTINE gcr(ldt, iter, x, b, r, iz, ibuff)
!!-------------------------------------------------------------------------------
!!
!!    Generalized conjugate gradient method (A*x = b for Non-Hermitian A)
!!
!!-------------------------------------------------------------------------------
!
!    real(kind=DP), intent(in) :: ldt
!    integer, intent(in) :: iter
!    complex(kind=DP), intent(inout), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: x, r
!    complex(kind=DP), intent(in), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: b
!    integer, intent(in) :: iz, ibuff
!
!    complex(kind=DP), save, &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1,0:iter_max-1) :: p, ap
!!$OMP threadprivate(p,ap)
!    complex(kind=DP), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: ar
!    complex(kind=DP) :: alpha, beta, sum1
!    complex(kind=DP), save, &
!      dimension(0:iter_max-1) :: sum2
!!$OMP threadprivate(sum2)
!    integer :: i
!  
!      if (iter == 0) then                           ! r_0 = b - A*x_0, p_0 = r_0
!        call calc_r(ldt, x, b, r, iz, ibuff)
!        p(:,:,:,iter) = r(:,:,:)
!        call calc_ap(ldt, p(:,:,:,iter), ap(:,:,:,iter), iz, ibuff)
!      end if
!  
!      sum1 = sum(ap(:,:,:,iter) * r(:,:,:))
!      sum2(iter) = sum(ap(:,:,:,iter) * ap(:,:,:,iter))
!      alpha = sum1 / sum2(iter)            ! alpha = (A*p_k,r_k) / (A*p_k,A*p_k)
!  
!      x(:,:,:) = x(:,:,:) + alpha * p(:,:,:,iter)    ! x_k+1 = x_k + alpha*p_k
!      r(:,:,:) = r(:,:,:) - alpha * ap(:,:,:,iter)   ! r_k+1 = r_k - alpha*A*p_k
!  
!      call calc_ap(ldt, r, ar, iz, ibuff)
!
!      p(:,:,:,iter+1) = r(:,:,:)
!      ap(:,:,:,iter+1) = ar(:,:,:)
!      do i = 0, iter
!        sum1 = sum(ap(:,:,:,i) * ar(:,:,:))
!        beta = - sum1 / sum2(i)     ! beta_i = - (A*p_i,A*r_k+1) / (A*p_i,A*p_i)
!        p(:,:,:,iter+1) = p(:,:,:,iter+1) + beta * p(:,:,:,i)
!                                            ! p_k+1 = r_k+1 + Sum_i^k beta_i*p_i
!        ap(:,:,:,iter+1) = ap(:,:,:,iter+1) + beta * ap(:,:,:,i)
!                                      ! A*p_k+1 = A*r_k+1 + Sum_i^k beta_i*A*p_i
!      end do
!  
!  END SUBROUTINE gcr
!
!
!!--------------------------------------
!  SUBROUTINE gcr1(ldt, iter, x, b, r, iz, ibuff)
!!-------------------------------------------------------------------------------
!!
!!    Generalized conjugate gradient method (A*x = b for Non-Hermitian A)
!!
!!-------------------------------------------------------------------------------
!
!    real(kind=DP), intent(in) :: ldt
!    integer, intent(in) :: iter
!    complex(kind=DP), intent(inout), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: x, r
!    complex(kind=DP), intent(in), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: b
!    integer, intent(in) :: iz, ibuff
!
!    complex(kind=DP), save, &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: p, ap
!!$OMP threadprivate(p,ap)
!    complex(kind=DP), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: ar
!    complex(kind=DP) :: alpha, beta, sum1, sum2
!  
!      if (iter == 0) then                           ! r_0 = b - A*x_0, p_0 = r_0
!        call calc_r(ldt, x, b, r, iz, ibuff)
!        p(:,:,:) = r(:,:,:)
!        call calc_ap(ldt, p, ap, iz, ibuff)
!      end if
!  
!      sum1 = sum(ap(:,:,:) * r(:,:,:))
!      sum2 = sum(ap(:,:,:) * ap(:,:,:))
!      alpha = sum1 / sum2                  ! alpha = (A*p_k,r_k) / (A*p_k,A*p_k)
!  
!      x(:,:,:) = x(:,:,:) + alpha * p(:,:,:)         ! x_k+1 = x_k + alpha*p_k
!      r(:,:,:) = r(:,:,:) - alpha * ap(:,:,:)        ! r_k+1 = r_k - alpha*A*p_k
!  
!      call calc_ap(ldt, r, ar, iz, ibuff)
!      sum1 = sum(ap(:,:,:) * ar(:,:,:))
!      beta = - sum1 / sum2            ! beta = - (A*p_k,A*r_k+1) / (A*p_k,A*p_k)
!  
!      p(:,:,:) = r(:,:,:) + beta * p(:,:,:)           ! p_k+1 = r_k+1 + beta*p_k
!      ap(:,:,:) = ar(:,:,:) + beta * ap(:,:,:)  ! A*p_k+1 = A*r_k+1 + beta*A*p_k
!  
!  END SUBROUTINE gcr1
!
!
!!--------------------------------------
!  SUBROUTINE gcrm(ldt, iter, x, b, r, iz, ibuff)
!!-------------------------------------------------------------------------------
!!
!!    Generalized conjugate gradient method (A*x = b for Non-Hermitian A)
!!
!!-------------------------------------------------------------------------------
!
!    real(kind=DP), intent(in) :: ldt
!    integer, intent(in) :: iter
!    complex(kind=DP), intent(inout), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: x, r
!    complex(kind=DP), intent(in), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: b
!    integer, intent(in) :: iz, ibuff
!
!    integer, parameter :: m = 10
!    complex(kind=DP), save, &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1,0:m-1) :: p, ap
!!$OMP threadprivate(p,ap)
!    complex(kind=DP), &
!      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: ar
!    complex(kind=DP) :: alpha, beta, sum1
!    complex(kind=DP), save, &
!      dimension(0:m-1) :: sum2
!!$OMP threadprivate(sum2)
!    integer :: i, k
!  
!      k = mod(iter,m)
!      if (k == 0) then                    ! r_0 = b - A*x_0, p_0 = r_0
!        call calc_r(ldt, x, b, r, iz, ibuff)
!        p(:,:,:,k) = r(:,:,:)
!        call calc_ap(ldt, p(:,:,:,k), ap(:,:,:,k), iz, ibuff)
!      end if
!  
!      sum1 = sum(ap(:,:,:,k) * r(:,:,:))
!      sum2(k) = sum(ap(:,:,:,k) * ap(:,:,:,k))
!      alpha = sum1 / sum2(k)            ! alpha = (A*p_k,r_k) / (A*p_k,A*p_k)
!  
!      x(:,:,:) = x(:,:,:) + alpha * p(:,:,:,k)    ! x_k+1 = x_k + alpha*p_k
!      r(:,:,:) = r(:,:,:) - alpha * ap(:,:,:,k)   ! r_k+1 = r_k - alpha*A*p_k
!  
!      call calc_ap(ldt, r, ar, iz, ibuff)
!
!      p(:,:,:,k+1) = r(:,:,:)
!      ap(:,:,:,k+1) = ar(:,:,:)
!      do i = 0, k
!        sum1 = sum(ap(:,:,:,i) * ar(:,:,:))
!        beta = - sum1 / sum2(i)     ! beta_i = - (A*p_i,A*r_k+1) / (A*p_i,A*p_i)
!        p(:,:,:,k+1) = p(:,:,:,k+1) + beta * p(:,:,:,i)
!                                            ! p_k+1 = r_k+1 + Sum_i^k beta_i*p_i
!        ap(:,:,:,k+1) = ap(:,:,:,k+1) + beta * ap(:,:,:,i)
!                                      ! A*p_k+1 = A*r_k+1 + Sum_i^k beta_i*A*p_i
!      end do
!  
!  END SUBROUTINE gcrm


!--------------------------------------
  SUBROUTINE calc_b(ldt, wh, b, iz, ibuff)
!-------------------------------------------------------------------------------
!
!    Calculate b
!
!-------------------------------------------------------------------------------

    real(kind=DP), intent(in) :: ldt
    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: wh
    complex(kind=DP), intent(out), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: b
    integer, intent(in) :: iz, ibuff
 
    complex(kind=DP), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: ww
    complex(kind=DP) :: phi, Al

      call emfield_hh(wh, Al, iz, ibuff)
      call esfield(wh, phi, iz, ibuff)

      if (iFLR == 0) then
        call setww_ff(wh, Al, ww, iz, ibuff) ! Gyrocenter distribution
      else
        call setww_gg(wh, phi, Al, ww, iz, ibuff) ! Non-adiabatic distribution
      end if

      if (trim(col_type) == "LB") then
        call collision_LB(ww, b, iz, ibuff)
      else if (trim(col_type) == "lorentz") then
        call collision_lorentz(ww, b, iz, ibuff)
      else if (trim(col_type) == "full") then
        call collision_full(ww, b, iz, ibuff)
      else
        write(olog,*) "# colliimp module does not support col_type=", col_type
        stop
      end if

      b(:,:,:) = wh(:,:,:) + 0.5_DP * ldt * b(:,:,:) ! 2nd C.N.
     !b(:,:,:) = wh(:,:,:)                           ! 1st B.E.
  
  END SUBROUTINE calc_b


!--------------------------------------
  SUBROUTINE calc_ap(ldt, p, ap, iz, ibuff)
!-------------------------------------------------------------------------------
!
!    Calculate A*p^k
!
!-------------------------------------------------------------------------------

    real(kind=DP), intent(in) :: ldt
    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: p
    complex(kind=DP), intent(out), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: ap
    integer, intent(in) :: iz, ibuff
 
    complex(kind=DP), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: ww
    complex(kind=DP) :: phi, Al

      call emfield_hh(p, Al, iz, ibuff)
      call esfield(p, phi, iz, ibuff)

      if (iFLR == 0) then
        call setww_ff(p, Al, ww, iz, ibuff) ! Gyrocenter distribution
      else
        call setww_gg(p, phi, Al, ww, iz, ibuff) ! Non-adiabatic distribution
      end if

      if (trim(col_type) == "LB") then
        call collision_LB(ww, ap, iz, ibuff)
      else if (trim(col_type) == "lorentz") then
        call collision_lorentz(ww, ap, iz, ibuff)
      else if (trim(col_type) == "full") then
        call collision_full(ww, ap, iz, ibuff)
      else
        write(olog,*) "# colliimp module does not support col_type=", col_type
        stop
      end if

      ap(:,:,:) = p(:,:,:) - 0.5_DP * ldt * ap(:,:,:) ! 2nd C.N.
     !ap(:,:,:) = p(:,:,:) - ldt * ap(:,:,:)          ! 1st B.E.
  
  END SUBROUTINE calc_ap


!--------------------------------------
  SUBROUTINE calc_r(ldt, x, b, r, iz, ibuff)
!-------------------------------------------------------------------------------
!
!    Calculate residual  r^n = b - A*x^n
!
!-------------------------------------------------------------------------------

    real(kind=DP), intent(in) :: ldt
    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: x, b
    complex(kind=DP), intent(out), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: r
    integer, intent(in) :: iz, ibuff
  
      call calc_ap(ldt, x, r, iz, ibuff)
      r(:,:,:) = b(:,:,:) - r(:,:,:)
  
  END SUBROUTINE calc_r
  
  
!--------------------------------------
  SUBROUTINE calc_error(r, b, error)
!-------------------------------------------------------------------------------
!
!    Evaluate relative error of A*x^n = b
!
!      When analytic solution x^* is known,
!
!        Relative error = ||x^n - x^*|| / ||x^*||
!
!      When x^* is not available, error is often replaced by followings:
!
!        Relative residual error = ||r^n|| / ||b||
!        Relative iteration error = ||x^n - x^n-1|| / ||x^n||
!
!      where residual is  r^n = b - A*x^n.
!
!      Several types of norm of the vector ||x|| are available,
!
!              1-norm  ||x|| = sum_i |x_i|
!              2-norm  ||x|| = sum_i |x_i|^2
!              p-norm  ||x|| = sum_i |x_i|^p
!        Maximum norm  ||x|| = max(|x_i|)
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: r, b
    real(kind=DP), intent(out) :: error

     !error = sum(abs(r(:,:,:))) / sum(abs(b(:,:,:)))
     !error = sum(abs(r(:,:,:))**2) / sum(abs(b(:,:,:))**2)
      error = maxval(abs(r(:,:,:))) / maxval(abs(b(:,:,:)))

  END SUBROUTINE calc_error


!--------------------------------------
  SUBROUTINE emfield_hh(hh, Al, iz, ibuff)
!--------------------------------------
!   Magnetic field calculation

    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: hh
    complex(kind=DP), intent(out) :: Al
    integer, intent(in) :: iz, ibuff
 
    complex(kind=DP) :: jpara, wf, wfvp, wfvp1
    real(kind=DP) :: cs2, cintgrl
    integer :: iv, im, is, mx, my

      jpara = (0._DP, 0._DP)
      do is = 0, ns-1
        cs2 = sgn(is) * fcs(is) * sqrt(tau(is) / Anum(is))

        do im = 1, global_nm
          cintgrl = 2._DP * pi * gvp(im,iz) * dv * dvp(iz)
          do iv = 1, 2*global_nv
            wf = hh(iv,im,is) * gj0(im,is,iz,ibuff) * cs2 * gvl(iv)
            jpara = jpara + wf * cintgrl
          end do
        end do

       !- edge compensation -
        im = 1
          do iv = 1, 2*global_nv
            wfvp  = gvp(im  ,iz) * hh(iv,im  ,is) * gj0(im  ,is,iz,ibuff) &
                  * cs2 * gvl(iv)
            wfvp1 = gvp(im+1,iz) * hh(iv,im+1,is) * gj0(im+1,is,iz,ibuff) &
                  * cs2 * gvl(iv)
            jpara = jpara - ( - wfvp / 12._DP       &
                            + ( wfvp1               &
                              - wfvp * 2._DP        &
                              ) * 11._DP / 720._DP  &
                            ) * (2._DP * pi * dv * dvp(iz))
          end do

      end do

      call ibuffiproc2mxmy(ibuff, spc_rank, mx, my)
      Al = jpara * beta * fct_ampere(mx,my,iz)

  END SUBROUTINE emfield_hh


!--------------------------------------
  SUBROUTINE esfield(hh, phi, iz, ibuff)
!--------------------------------------
!   Electrostatic field calculation

    implicit none
    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: hh
    complex(kind=DP), intent(out) :: phi
    integer, intent(in) :: iz, ibuff
 
    complex(kind=DP) :: rho, wf, wfvp, wfvp1
    real(kind=DP) :: cs1, cintgrl
    integer :: iv, im, is, mx, my

      rho = (0._DP, 0._DP)
      do is = 0, ns-1
        cs1 = sgn(is) * fcs(is)

        do im = 1, global_nm
          cintgrl = 2._DP * pi * gvp(im,iz) * dv * dvp(iz)
          do iv = 1, 2*global_nv
            wf = hh(iv,im,is) * gj0(im,is,iz,ibuff) * cs1
            rho = rho + wf * cintgrl
          end do
        end do

       !- edge compensation -
        im = 1
          do iv = 1, 2*global_nv
            wfvp  = gvp(im  ,iz) * hh(iv,im  ,is) * gj0(im  ,is,iz,ibuff) * cs1
            wfvp1 = gvp(im+1,iz) * hh(iv,im+1,is) * gj0(im+1,is,iz,ibuff) * cs1
            rho = rho - ( - wfvp / 12._DP       &
                        + ( wfvp1               &
                          - wfvp * 2._DP        &
                          ) * 11._DP / 720._DP  &
                        ) * (2._DP * pi * dv * dvp(iz))
          end do

      end do

      call ibuffiproc2mxmy(ibuff, spc_rank, mx, my)
      phi = rho * fct_poisson(mx,my,iz)

  END SUBROUTINE esfield


!--------------------------------------
  SUBROUTINE setww_ff(hh, Al, ff, iz, ibuff)
!--------------------------------------
!   Set gyrocenter distribution ff

    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: hh
    complex(kind=DP), intent(in) :: Al
    complex(kind=DP), intent(out), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: ff
    integer, intent(in) :: iz, ibuff
 
    real(kind=DP) :: cs2
    integer :: iv, im, is, ivb

      do is = 0, ns-1
        cs2 = sgn(is) * Znum(is) / sqrt( Anum(is) * tau(is) )
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            ff(iv,im,is) = hh(iv,im,is) - gfmx(iv,im,iz) * gj0(im,is,iz,ibuff) &
                                                          * cs2 * gvl(iv) * Al
          end do
        end do
      end do

    !- Boundary condition -
      do ivb = 1, nvb
        ff(1-ivb,:,:) = (0._DP, 0._DP)
        ff(2*global_nv+ivb,:,:) = (0._DP, 0._DP)
        ff(:,global_nm+ivb,:) = (0._DP, 0._DP)
        ff(:,-ivb,:) = ff(:,ivb,:)
      end do
    !-

  END SUBROUTINE setww_ff


!--------------------------------------
  SUBROUTINE setww_gg(hh, phi, Al, gg, iz, ibuff)
!--------------------------------------
!   Set non-adiabatic distribution gg

    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: hh
    complex(kind=DP), intent(in) :: phi, Al
    complex(kind=DP), intent(out), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: gg
    integer, intent(in) :: iz, ibuff
 
    real(kind=DP) :: cs1, cs2
    integer :: iv, im, is, ivb

      do is = 0, ns-1
        cs1 = sgn(is) * Znum(is) / tau(is)
        cs2 = sgn(is) * Znum(is) / sqrt( Anum(is) * tau(is) )
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            gg(iv,im,is) = hh(iv,im,is) + gfmx(iv,im,iz) * gj0(im,is,iz,ibuff) &
                                            * (cs1 * phi - cs2 * gvl(iv) * Al)
          end do
        end do
      end do

    !- Boundary condition -
      do ivb = 1, nvb
        gg(1-ivb,:,:) = (0._DP, 0._DP)
        gg(2*global_nv+ivb,:,:) = (0._DP, 0._DP)
        gg(:,global_nm+ivb,:) = (0._DP, 0._DP)
        gg(:,-ivb,:) = gg(:,ivb,:)
      end do
    !-

  END SUBROUTINE setww_gg


!--------------------------------------
  SUBROUTINE collision_LB(wf, cf, iz, ibuff)
!--------------------------------------
!   Lenard-Bernstein operator

    complex(kind=DP), intent(in), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: wf
    complex(kind=DP), intent(out), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: cf
    integer, intent(in) :: iz, ibuff
  
    real(kind=DP) :: nu_s, cv1, cv2, cm1, cm2, cflr
    integer :: iv, im, is, mxy, mx, my

      mxy = ibuff + nbuff * spc_rank
      if (mxy <= (2*nx+1)*(ny+1)-1) then

        mx = mod(mxy,2*nx+1) - nx
        my = mxy / (2*nx+1)
        do is = 0, ns-1
          nu_s = nu(is) * 3._DP * sqrt(pi) * ctauiv(is,is) / 4._DP
               !- nu(is) is a bias factor given in namelist
          cv1 = nu_s / ( 12._DP * dv )
          cv2 = nu_s / ( 12._DP * dv * dv )
          cm1 = nu_s / ( 12._DP * dvp(iz) )
          cm2 = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
          cflr = nu_s * ksq(mx,my,iz) * Anum(is) * tau(is)  &
                 / ( Znum(is) * omg(iz) )**2 * real(iFLR, kind=DP)
          im = 0 
            do iv = 1, 2*global_nv
              cf(iv,im,is) =                                        &
                         ( -          wf(iv+2,im,is)                &
                           + 16._DP * wf(iv+1,im,is)                &
                           - 30._DP * wf(iv  ,im,is)                &
                           + 16._DP * wf(iv-1,im,is)                &
                           -          wf(iv-2,im,is)                &
                         ) * cv2                                    &
                       + ( -          wf(iv+2,im,is)                &
                           +  8._DP * wf(iv+1,im,is)                &
                           -  8._DP * wf(iv-1,im,is)                &
                           +          wf(iv-2,im,is)                &
                         ) * cv1 * gvl(iv)                          &
                       + ( -          wf(iv,im+2,is)                &
                           + 16._DP * wf(iv,im+1,is)                &
                           - 30._DP * wf(iv,im  ,is)                &
                           + 16._DP * wf(iv,im-1,is)                &
                           -          wf(iv,im-2,is)                &
                         ) * cm2 * 2._DP                            &
                       + nu_s * 3._DP * wf(iv,im,is)                &
                       - cflr * wf(iv,im,is)
            end do
          do im = 1, global_nm
            do iv = 1, 2*global_nv
              cf(iv,im,is) =                                            &
                         ( -          wf(iv+2,im,is)                    &
                           + 16._DP * wf(iv+1,im,is)                    &
                           - 30._DP * wf(iv  ,im,is)                    &
                           + 16._DP * wf(iv-1,im,is)                    &
                           -          wf(iv-2,im,is)                    &
                         ) * cv2                                        &
                       + ( -          wf(iv+2,im,is)                    &
                           +  8._DP * wf(iv+1,im,is)                    &
                           -  8._DP * wf(iv-1,im,is)                    &
                           +          wf(iv-2,im,is)                    &
                         ) * cv1 * gvl(iv)                              &
                       + ( -          wf(iv,im+2,is)                    &
                           + 16._DP * wf(iv,im+1,is)                    &
                           - 30._DP * wf(iv,im  ,is)                    &
                           + 16._DP * wf(iv,im-1,is)                    &
                           -          wf(iv,im-2,is)                    &
                         ) * cm2                                        &
                       + ( -          wf(iv,im+2,is)                    &
                           +  8._DP * wf(iv,im+1,is)                    &
                           -  8._DP * wf(iv,im-1,is)                    &
                           +          wf(iv,im-2,is)                    &
                         ) * cm1 * ( gvp(im,iz) + 1._DP / gvp(im,iz) )  &
                       + nu_s * 3._DP * wf(iv,im,is)                    &   
                       - cflr * wf(iv,im,is)
            end do
          end do
        end do

      else

        cf(:,:,:) = (0._DP, 0._DP)

      end if

  END SUBROUTINE collision_LB


!--------------------------------------
  SUBROUTINE collision_lorentz(wf, cf, iz, ibuff)
!--------------------------------------
!   Lorentz operator

    complex(kind=DP), intent(in), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: wf
    complex(kind=DP), intent(out), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: cf
    integer, intent(in) :: iz, ibuff
  
    real(kind=DP) :: cv1, cv2, cm1, cm2, cvm, cflr
    integer :: iv, im, is, mxy, mx, my

      mxy = ibuff + nbuff * spc_rank
      if (mxy <= (2*nx+1)*(ny+1)-1) then

        mx = mod(mxy,2*nx+1) - nx
        my = mxy / (2*nx+1)
        do is = 0, ns-1
          cv1 = 1._DP / (12._DP * dv)
          cv2 = 1._DP / (12._DP * dv**2)
          cm1 = 1._DP / (12._DP * dvp(iz))
          cm2 = 1._DP / (12._DP * dvp(iz)**2)
          cvm = 1._DP / (144._DP * dv * dvp(iz))
          cflr = ksq(mx,my,iz) * Anum(is) * tau(is)  &
                 / ( Znum(is) * omg(iz) )**2 * real(iFLR, kind=DP)

          im = 0
            do iv = 1, 2*global_nv
              cf(iv,im,is) = 0.5_DP * gnu_ds(iv,im,is,iz) * (                  &
                          - 2._DP * gvl(iv) * cv1 * (-         wf(iv+2,im,is)  &
                                                     + 8._DP * wf(iv+1,im,is)  &
                                                     - 8._DP * wf(iv-1,im,is)  &
                                                     +         wf(iv-2,im,is)) &
                       + 2._DP * gvl(iv)**2 * cm2 * (-         wf(iv,im+2,is)  &
                                                     +16._DP * wf(iv,im+1,is)  &
                                                     -30._DP * wf(iv,im  ,is)  &
                                                     +16._DP * wf(iv,im-1,is)  &
                                                     -         wf(iv,im-2,is)) &
                           )                                                   &
                           - 0.25_DP * gnu_ds(iv,im,is,iz)                     &
                             * (2._DP * gvl(iv)**2)                            &
                             * cflr * wf(iv,im,is)
            end do

          do im = 1, global_nm
            do iv = 1, 2*global_nv
              cf(iv,im,is) = 0.5_DP * gnu_ds(iv,im,is,iz) * (                  &
                              gvp(im,iz)**2 * cv2 * (-         wf(iv+2,im,is)  &
                                                     +16._DP * wf(iv+1,im,is)  &
                                                     -30._DP * wf(iv  ,im,is)  &
                                                     +16._DP * wf(iv-1,im,is)  &
                                                     -         wf(iv-2,im,is)) &
                          - 2._DP * gvl(iv) * cv1 * (-         wf(iv+2,im,is)  &
                                                     + 8._DP * wf(iv+1,im,is)  &
                                                     - 8._DP * wf(iv-1,im,is)  &
                                                     +         wf(iv-2,im,is)) &
                           - 2._DP * gvl(iv) * gvp(im,iz) &
                                          * cvm * (+         wf(iv+2,im+2,is)  &
                                                   - 8._DP * wf(iv+2,im+1,is)  &
                                                   + 8._DP * wf(iv+2,im-1,is)  &
                                                   -         wf(iv+2,im-2,is)  &
                                                   - 8._DP * wf(iv+1,im+2,is)  &
                                                   +64._DP * wf(iv+1,im+1,is)  &
                                                   -64._DP * wf(iv+1,im-1,is)  &
                                                   + 8._DP * wf(iv+1,im-2,is)  &
                                                   + 8._DP * wf(iv-1,im+2,is)  &
                                                   -64._DP * wf(iv-1,im+1,is)  &
                                                   +64._DP * wf(iv-1,im-1,is)  &
                                                   - 8._DP * wf(iv-1,im-2,is)  &
                                                   -         wf(iv-2,im+2,is)  &
                                                   + 8._DP * wf(iv-2,im+1,is)  &
                                                   - 8._DP * wf(iv-2,im-1,is)  &
                                                   +         wf(iv-2,im-2,is)) &
                           + ((gvl(iv)**2 - gvp(im,iz)**2) / gvp(im,iz)) &
                                            * cm1 * (-         wf(iv,im+2,is)  &
                                                     + 8._DP * wf(iv,im+1,is)  &
                                                     - 8._DP * wf(iv,im-1,is)  &
                                                     +         wf(iv,im-2,is)) &
                               + gvl(iv)**2 * cm2 * (-         wf(iv,im+2,is)  &
                                                     +16._DP * wf(iv,im+1,is)  &
                                                     -30._DP * wf(iv,im  ,is)  &
                                                     +16._DP * wf(iv,im-1,is)  &
                                                     -         wf(iv,im-2,is)) &
                           )                                                   &
                           - 0.25_DP * gnu_ds(iv,im,is,iz)                     &
                             * (2._DP * gvl(iv)**2 + gvp(im,iz)**2)            &
                             * cflr * wf(iv,im,is)
            end do
          end do
        end do

      else

        cf(:,:,:) = (0._DP, 0._DP)

      end if

  END SUBROUTINE collision_lorentz


!--------------------------------------
  SUBROUTINE collision_full(wf, cf, iz, ibuff)
!--------------------------------------
!   Sugama operator

    complex(kind=DP), intent(in), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: wf
    complex(kind=DP), intent(out), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: cf
    integer, intent(in) :: iz, ibuff

    complex(kind=DP), dimension(1:6,0:ns-1,0:ns-1) :: moment_ab

                             !%%% For debug %%%
                             ! complex(kind=DP), &
                             ! dimension(1:2*global_nv,0:global_nm,0:ns-1) &
                             !                              :: cft, cfd, cff
                             ! integer :: mx, my, iv, im
                             ! integer, save :: iflg_dbg
                             ! data iflg_dbg / 0 /
                             !%%%%%%%%%%%%%%%%%

      call collision_full_CT(wf, cf, iz, ibuff)

      call collision_full_calc_moment(wf, moment_ab, iz, ibuff)
      call collision_full_DT(moment_ab, cf, iz, ibuff)
      call collision_full_CF(moment_ab, cf, iz, ibuff)

                             !%%% For debug %%%
                             ! call ibuffiproc2mxmy&
                             !                   (ibuff,spc_rank,mx,my)
                             ! if (rankw == 0 .and. rankz == 0 .and. &
                             !     mx == 1 .and. my == 1 .and.       &
                             !     iz == 0 .and. iflg_dbg == 0) then
                             ! iflg_dbg = 1
                             ! call collision_full_ct(wf,cft,iz,ibuff)
                             ! call collision_full_dt(moment_ab,cfd,iz,ibuff)
                             ! call collision_full_cf(moment_ab,cff,iz,ibuff)
                             ! do im = 0, global_nm
                             ! do iv = 1, 2*global_nv
                             !   write(97,*) gvl(iv), gvp(im,iz), &
                             !              dble(wf(iv,im,0)),    &
                             !              aimag(wf(iv,im,0)),   &
                             !              dble(cft(iv,im,0)),   &
                             !              aimag(cft(iv,im,0)),  &
                             !              dble(cfd(iv,im,0)),   &
                             !              aimag(cfd(iv,im,0)),  &
                             !              dble(cff(iv,im,0)),   &
                             !              aimag(cff(iv,im,0))
                             ! end do
                             !   write(97,*)
                             ! end do
                             ! end if
                             !%%%%%%%%%%%%%%%%%

  END SUBROUTINE collision_full


!--------------------------------------
  SUBROUTINE collision_full_CT(wf, cf, iz, ibuff)
!--------------------------------------
!   Sugama operator, test particle part (scattering+slowing down+diffusion)

    complex(kind=DP), intent(in), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: wf
    complex(kind=DP), intent(out), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: cf
    integer, intent(in) :: iz, ibuff
  
    real(kind=DP) :: cv1, cv2, cm1, cm2, cvm, cflr
    integer :: iv, im, is, mxy, mx, my

      mxy = ibuff + nbuff * spc_rank
      if (mxy <= (2*nx+1)*(ny+1)-1) then

        mx = mod(mxy,2*nx+1) - nx
        my = mxy / (2*nx+1)
        do is = 0, ns-1

          cv1 = 1._DP / (12._DP * dv)
          cv2 = 1._DP / (12._DP * dv**2)
          cm1 = 1._DP / (12._DP * dvp(iz))
          cm2 = 1._DP / (12._DP * dvp(iz)**2)
          cvm = 1._DP / (144._DP * dv * dvp(iz))
          cflr = ksq(mx,my,iz) * Anum(is) * tau(is)  &
                 / ( Znum(is) * omg(iz) )**2 * real(iFLR, kind=DP)

          im = 0
            do iv = 1, 2*global_nv
              cf(iv,im,is) = ( -          wf(iv+2,im,is)                       &
                               + 16._DP * wf(iv+1,im,is)                       &
                               - 30._DP * wf(iv  ,im,is)                       &
                               + 16._DP * wf(iv-1,im,is)                       &
                               -          wf(iv-2,im,is)                       &
                             ) * cv2                                           &
                               * (   gnu_ps(iv,im,is,iz) * gvl(iv)**2          &
                                 ) * 0.5_DP                                    &
                           + ( -          wf(iv,im+2,is)                       &
                               + 16._DP * wf(iv,im+1,is)                       &
                               - 30._DP * wf(iv,im  ,is)                       &
                               + 16._DP * wf(iv,im-1,is)                       &
                               -          wf(iv,im-2,is)                       &
                             ) * cm2                                           &
                               * (   gnu_ds(iv,im,is,iz) * gvl(iv)**2          &
                                 )                                             &
                           + ( -          wf(iv+2,im,is)                       &
                               +  8._DP * wf(iv+1,im,is)                       &
                               -  8._DP * wf(iv-1,im,is)                       &
                               +          wf(iv-2,im,is)                       &
                             ) * cv1                                           &
                               * gnu_gs(iv,im,is,iz) * gvl(iv)                 &
                          !+ ( gnu_hs(iv,im,is,iz) * gxxa**2 * 2._DP           &
                           + ( gnu_hs(iv,im,is,iz) * 2._DP                     &
                                - 0.25_DP * cflr                               &
                                 * gnu_ds(iv,im,is,iz) * 2._DP * gvl(iv)**2    &
                             ) * wf(iv,im,is)
            end do

          do im = 1, global_nm
            do iv = 1, 2*global_nv
              cf(iv,im,is) = ( -          wf(iv+2,im,is)                       &
                               + 16._DP * wf(iv+1,im,is)                       &
                               - 30._DP * wf(iv  ,im,is)                       &
                               + 16._DP * wf(iv-1,im,is)                       &
                               -          wf(iv-2,im,is)                       &
                             ) * cv2                                           &
                               * (   gnu_ps(iv,im,is,iz) * gvl(iv)**2          &
                                   + gnu_ds(iv,im,is,iz) * gvp(im,iz)**2       &
                                 ) * 0.5_DP                                    &
                           + ( -          wf(iv,im+2,is)                       &
                               + 16._DP * wf(iv,im+1,is)                       &
                               - 30._DP * wf(iv,im  ,is)                       &
                               + 16._DP * wf(iv,im-1,is)                       &
                               -          wf(iv,im-2,is)                       &
                             ) * cm2                                           &
                               * (   gnu_ds(iv,im,is,iz) * gvl(iv)**2          &
                                   + gnu_ps(iv,im,is,iz) * gvp(im,iz)**2       &
                                 ) * 0.5_DP                                    &
                           + ( +          wf(iv+2,im+2,is)                     &
                               -  8._DP * wf(iv+2,im+1,is)                     &
                               +  8._DP * wf(iv+2,im-1,is)                     &
                               -          wf(iv+2,im-2,is)                     &
                               -  8._DP * wf(iv+1,im+2,is)                     &
                               + 64._DP * wf(iv+1,im+1,is)                     &
                               - 64._DP * wf(iv+1,im-1,is)                     &
                               +  8._DP * wf(iv+1,im-2,is)                     &
                               +  8._DP * wf(iv-1,im+2,is)                     &
                               - 64._DP * wf(iv-1,im+1,is)                     &
                               + 64._DP * wf(iv-1,im-1,is)                     &
                               -  8._DP * wf(iv-1,im-2,is)                     &
                               -          wf(iv-2,im+2,is)                     &
                               +  8._DP * wf(iv-2,im+1,is)                     &
                               -  8._DP * wf(iv-2,im-1,is)                     &
                               +          wf(iv-2,im-2,is)                     &
                             ) * cvm                                           &
                               * gvl(iv) * gvp(im,iz)                          &
                               * (   gnu_ps(iv,im,is,iz)                       &
                                   - gnu_ds(iv,im,is,iz) )                     &
                           + ( -          wf(iv+2,im,is)                       &
                               +  8._DP * wf(iv+1,im,is)                       &
                               -  8._DP * wf(iv-1,im,is)                       &
                               +          wf(iv-2,im,is)                       &
                             ) * cv1                                           &
                               * gnu_gs(iv,im,is,iz) * gvl(iv)                 &
                           + ( -          wf(iv,im+2,is)                       &
                               +  8._DP * wf(iv,im+1,is)                       &
                               -  8._DP * wf(iv,im-1,is)                       &
                               +          wf(iv,im-2,is)                       &
                             ) * cm1                                           &
                               * (   gnu_gs(iv,im,is,iz) * gvp(im,iz)          &
                                   + gnu_ds(iv,im,is,iz) * 0.5_DP              &
                                     * (gvl(iv)**2 / gvp(im,iz) + gvp(im,iz)) )&
                          !+ ( gnu_hs(iv,im,is,iz) * gxxa**2 * 2._DP           &
                           + ( gnu_hs(iv,im,is,iz) * 2._DP                     &
                                - 0.25_DP * cflr                               &
                                 * (  gnu_ds(iv,im,is,iz)                      &
                                       * (2._DP * gvl(iv)**2 + gvp(im,iz)**2)  &
                                     + gnu_ps(iv,im,is,iz) * gvp(im,iz)**2 )   &
                             ) * wf(iv,im,is)
            end do
          end do
        end do

      else

        cf(:,:,:) = (0._DP, 0._DP)

      end if

  END SUBROUTINE collision_full_CT


!--------------------------------------
  SUBROUTINE collision_full_calc_moment(wf, moment_ab, iz, ibuff)
!--------------------------------------
!   Sugama operator, calculate velocity moments

    complex(kind=DP), intent(in), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: wf
    complex(kind=DP), intent(out), &
      dimension(1:6,0:ns-1,0:ns-1) :: moment_ab
    integer, intent(in) :: iz, ibuff

    complex(kind=DP) :: wfvp, wfvp1
    integer :: iv, im, ia, ib
             
      moment_ab(:,:,:) = (0._DP, 0._DP)

      if ( iFLR == 1 ) then  ! full-GK

        do ib = 0, ns-1
          do ia = 0, ns-1
  
            do im = 1, global_nm
              do iv = 1, 2*global_nv
                  moment_ab(1,ia,ib) = moment_ab(1,ia,ib) + wf(iv,im,ia) * gj0(im,ia,iz,ibuff) * gvfunc(1,iv,im,ia,ib,iz)
                  moment_ab(2,ia,ib) = moment_ab(2,ia,ib) + wf(iv,im,ia) * gj1(im,ia,iz,ibuff) * gvfunc(2,iv,im,ia,ib,iz)
                  moment_ab(3,ia,ib) = moment_ab(3,ia,ib) + wf(iv,im,ia) * gj0(im,ia,iz,ibuff) * gvfunc(3,iv,im,ia,ib,iz)
                  moment_ab(4,ia,ib) = moment_ab(4,ia,ib) + wf(iv,im,ia) * gj0(im,ia,iz,ibuff) * gvfunc(4,iv,im,ia,ib,iz)
                  moment_ab(5,ia,ib) = moment_ab(5,ia,ib) + wf(iv,im,ia) * gj1(im,ia,iz,ibuff) * gvfunc(5,iv,im,ia,ib,iz)
                  moment_ab(6,ia,ib) = moment_ab(6,ia,ib) + wf(iv,im,ia) * gj0(im,ia,iz,ibuff) * gvfunc(6,iv,im,ia,ib,iz)
              end do
            end do
  
           !- edge compensation -
            im = 1
              do iv = 1, 2*global_nv
                wfvp  = wf(iv,im  ,ia) * gj0(im  ,ia,iz,ibuff) * gvfunc(1,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gj0(im+1,ia,iz,ibuff) * gvfunc(1,iv,im+1,ia,ib,iz)
                moment_ab(1,ia,ib) = moment_ab(1,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
                wfvp  = wf(iv,im  ,ia) * gj1(im  ,ia,iz,ibuff) * gvfunc(2,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gj1(im+1,ia,iz,ibuff) * gvfunc(2,iv,im+1,ia,ib,iz)
                moment_ab(2,ia,ib) = moment_ab(2,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
                wfvp  = wf(iv,im  ,ia) * gj0(im  ,ia,iz,ibuff) * gvfunc(3,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gj0(im+1,ia,iz,ibuff) * gvfunc(3,iv,im+1,ia,ib,iz)
                moment_ab(3,ia,ib) = moment_ab(3,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
                wfvp  = wf(iv,im  ,ia) * gj0(im  ,ia,iz,ibuff) * gvfunc(4,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gj0(im+1,ia,iz,ibuff) * gvfunc(4,iv,im+1,ia,ib,iz)
                moment_ab(4,ia,ib) = moment_ab(4,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
                wfvp  = wf(iv,im  ,ia) * gj1(im  ,ia,iz,ibuff) * gvfunc(5,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gj1(im+1,ia,iz,ibuff) * gvfunc(5,iv,im+1,ia,ib,iz)
                moment_ab(5,ia,ib) = moment_ab(5,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
                wfvp  = wf(iv,im  ,ia) * gj0(im  ,ia,iz,ibuff) * gvfunc(6,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gj0(im+1,ia,iz,ibuff) * gvfunc(6,iv,im+1,ia,ib,iz)
                moment_ab(6,ia,ib) = moment_ab(6,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
              end do
          
          end do
        end do

      else if ( iFLR == 0 ) then ! DK-limit

        do ib = 0, ns-1
          do ia = 0, ns-1
  
            do im = 1, global_nm
              do iv = 1, 2*global_nv
                  moment_ab(1,ia,ib) = moment_ab(1,ia,ib) + wf(iv,im,ia) * gvfunc(1,iv,im,ia,ib,iz)
                  moment_ab(3,ia,ib) = moment_ab(3,ia,ib) + wf(iv,im,ia) * gvfunc(3,iv,im,ia,ib,iz)
                  moment_ab(4,ia,ib) = moment_ab(4,ia,ib) + wf(iv,im,ia) * gvfunc(4,iv,im,ia,ib,iz)
                  moment_ab(6,ia,ib) = moment_ab(6,ia,ib) + wf(iv,im,ia) * gvfunc(6,iv,im,ia,ib,iz)
              end do
            end do
  
           !- edge compensation -
            im = 1
              do iv = 1, 2*global_nv
                wfvp  = wf(iv,im  ,ia) * gvfunc(1,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gvfunc(1,iv,im+1,ia,ib,iz)
                moment_ab(1,ia,ib) = moment_ab(1,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
                wfvp  = wf(iv,im  ,ia) * gvfunc(3,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gvfunc(3,iv,im+1,ia,ib,iz)
                moment_ab(3,ia,ib) = moment_ab(3,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
                wfvp  = wf(iv,im  ,ia) * gvfunc(4,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gvfunc(4,iv,im+1,ia,ib,iz)
                moment_ab(4,ia,ib) = moment_ab(4,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
                wfvp  = wf(iv,im  ,ia) * gvfunc(6,iv,im  ,ia,ib,iz)
                wfvp1 = wf(iv,im+1,ia) * gvfunc(6,iv,im+1,ia,ib,iz)
                moment_ab(6,ia,ib) = moment_ab(6,ia,ib)        &
                                   - ( - wfvp / 12._DP         &
                                       + ( wfvp1               &
                                         - wfvp * 2._DP        &
                                         ) * 11._DP / 720._DP )
              end do
          
          end do
        end do

      end if

  END SUBROUTINE collision_full_calc_moment


!--------------------------------------
  SUBROUTINE collision_full_DT(moment_ab, cf, iz, ibuff)
!--------------------------------------
!   Sugama operator, test particle part (non-isothermal terms)

    complex(kind=DP), intent(in), &
      dimension(1:6,0:ns-1,0:ns-1) :: moment_ab
    complex(kind=DP), intent(inout), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: cf
    integer, intent(in) :: iz, ibuff
  
    integer :: iv, im, ia, ib, mxy

      mxy = ibuff + nbuff * spc_rank
      if (mxy <= (2*nx+1)*(ny+1)-1) then

        if ( iFLR == 1 ) then  ! full-GK
   
          do ib = 0, ns-1
            do ia = 0, ns-1
              do im = 0, global_nm
                do iv = 1, 2*global_nv
                  cf(iv,im,ia) = cf(iv,im,ia)             &
                               + gx_tst(1,iv,im,ia,ib,iz) &
                                * moment_ab(1,ia,ib)      &
                                 * gj0(im,ia,iz,ibuff)    & 
                               + gx_tst(2,iv,im,ia,ib,iz) &
                                * moment_ab(2,ia,ib)      &
                                 * gj1(im,ia,iz,ibuff)    & 
                               + gx_tst(3,iv,im,ia,ib,iz) &
                                * moment_ab(3,ia,ib)      &
                                 * gj0(im,ia,iz,ibuff)    & 
                               + gx_tst(4,iv,im,ia,ib,iz) &
                                * moment_ab(4,ia,ib)      &
                                 * gj0(im,ia,iz,ibuff)    & 
                               + gx_tst(5,iv,im,ia,ib,iz) &
                                * moment_ab(5,ia,ib)      &
                                 * gj1(im,ia,iz,ibuff)    & 
                               + gx_tst(6,iv,im,ia,ib,iz) &
                                * moment_ab(6,ia,ib)      &
                                 * gj0(im,ia,iz,ibuff)                  
                end do
              end do
            end do
          end do
   
        else if ( iFLR == 0 ) then ! DK-limit
   
          do ib = 0, ns-1
            do ia = 0, ns-1
              do im = 0, global_nm
                do iv = 1, 2*global_nv
                  cf(iv,im,ia) = cf(iv,im,ia)             &
                               + gx_tst(1,iv,im,ia,ib,iz) &
                                * moment_ab(1,ia,ib)      &
                               + gx_tst(3,iv,im,ia,ib,iz) &
                                * moment_ab(3,ia,ib)      &
                               + gx_tst(4,iv,im,ia,ib,iz) &
                                * moment_ab(4,ia,ib)      &
                               + gx_tst(6,iv,im,ia,ib,iz) &
                                * moment_ab(6,ia,ib)
                end do
              end do
            end do
          end do
   
        end if

      else

        !cf(:,:,:) = (0._DP, 0._DP)

      end if

  END SUBROUTINE collision_full_DT


!--------------------------------------
  SUBROUTINE collision_full_CF(moment_ab, cf, iz, ibuff)
!--------------------------------------
!   Sugama operator, field particle part

    complex(kind=DP), intent(in), &
      dimension(1:6,0:ns-1,0:ns-1) :: moment_ab
    complex(kind=DP), intent(inout), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: cf
    integer, intent(in) :: iz, ibuff
  
    integer :: iv, im, ia, ib, mxy

      mxy = ibuff + nbuff * spc_rank
      if (mxy <= (2*nx+1)*(ny+1)-1) then

        if ( iFLR == 1 ) then  ! full-GK

          do ib = 0, ns-1
            do ia = 0, ns-1
              do im = 0, global_nm
                do iv = 1, 2*global_nv
                  cf(iv,im,ia) = cf(iv,im,ia)             &
                               + gy_fld(1,iv,im,ia,ib,iz) &
                                * moment_ab(1,ib,ia)      &
                                 * gj0(im,ia,iz,ibuff)    & 
                               + gy_fld(2,iv,im,ia,ib,iz) &
                                * moment_ab(2,ib,ia)      &
                                 * gj1(im,ia,iz,ibuff)    & 
                               + gy_fld(3,iv,im,ia,ib,iz) &
                                * moment_ab(3,ib,ia)      &
                                 * gj0(im,ia,iz,ibuff)    & 
                               + gy_fld(4,iv,im,ia,ib,iz) &
                                * moment_ab(4,ib,ia)      &
                                 * gj0(im,ia,iz,ibuff)    & 
                               + gy_fld(5,iv,im,ia,ib,iz) &
                                * moment_ab(5,ib,ia)      &
                                 * gj1(im,ia,iz,ibuff)    & 
                               + gy_fld(6,iv,im,ia,ib,iz) &
                                * moment_ab(6,ib,ia)      &
                                 * gj0(im,ia,iz,ibuff)                  
                end do
              end do
            end do
          end do

        else if ( iFLR == 0 ) then ! DK-limit
   
          do ib = 0, ns-1
            do ia = 0, ns-1
              do im = 0, global_nm
                do iv = 1, 2*global_nv
                  cf(iv,im,ia) = cf(iv,im,ia)             &
                               + gy_fld(1,iv,im,ia,ib,iz) &
                                * moment_ab(1,ib,ia)      &
                               + gy_fld(3,iv,im,ia,ib,iz) &
                                * moment_ab(3,ib,ia)      &
                               + gy_fld(4,iv,im,ia,ib,iz) &
                                * moment_ab(4,ib,ia)      &
                               + gy_fld(6,iv,im,ia,ib,iz) &
                                * moment_ab(6,ib,ia)
                end do
              end do
            end do
          end do
   
        end if

      else

        !cf(:,:,:) = (0._DP, 0._DP)

      end if

  END SUBROUTINE collision_full_CF


!!!!--------------------------------------
!!!  SUBROUTINE colliimp_calc_colli_full( ff, phi, cf )
!!!!--------------------------------------
!!!!   Collsion operator calculation interface
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1)             :: phi
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf
!!!
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: sender, recver
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
!!!    complex(kind=DP), dimension(:,:,:), allocatable :: w1, w2
!!!    integer :: mx, my, iz, iv, im, ibuff, ivb
!!!
!!!      allocate( sender(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) )
!!!      allocate( recver(1:2*nv,0:nm,-nz:nz-1,0:nbuff-1,0:nprocvms-1) )
!!!      allocate( wf(1:2*global_nv,0:global_nm,0:ns-1,-nz:nz-1,0:nbuff-1) )
!!!      allocate( w1(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) )
!!!      allocate( w2(1:2*global_nv,0:global_nm,0:ns-1) )
!!!
!!!!$OMP parallel do collapse(2)
!!!      do im = 0, nm
!!!        do iv = 1, 2*nv
!!!          do iz = -nz, nz-1
!!!            do my = ist_y, iend_y
!!!              do mx = -nx, nx
!!!                cf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) + sgn(ranks) * Znum(ranks) &
!!!                          * fmx(iz,iv,im) / tau(ranks) * j0(mx,my,iz,im) * phi(mx,my,iz) &
!!!                                                                     * real(iFLR, kind=DP)
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!
!!!      call vms2xy_pack(cf, sender)
!!!      call vms2xy_transpose(sender, recver)
!!!      call vms2xy_unpack(recver, wf)
!!!      
!!!                                      call clock_sta(1720)
!!!!$OMP parallel default(none) &
!!!!$OMP shared(wf) &
!!!!$OMP private(w1,w2,iz,ibuff)
!!!      w1(:,:,:) = (0._DP, 0._DP)
!!!!$OMP do collapse(2)
!!!      do ibuff = 0, nbuff-1
!!!        do iz = -nz, nz-1
!!!          w1(1:2*global_nv,0:global_nm,0:ns-1) = wf(:,:,:,iz,ibuff)
!!!         !- Boundary condition -
!!!          do ivb = 1, nvb
!!!           !w1(1-ivb,:,:) = (0._DP, 0._DP)
!!!           !w1(2*global_nv+ivb,:,:) = (0._DP, 0._DP)
!!!           !w1(:,global_nm+ivb,:) = (0._DP, 0._DP)
!!!            w1(:,-ivb,:) = w1(:,ivb,:)
!!!          end do
!!!         !-
!!!          call collision_full(w1, w2, iz, ibuff)
!!!          wf(:,:,:,iz,ibuff) = w2(:,:,:)
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!!$OMP end parallel
!!!                                      call clock_end(1720)
!!!
!!!      call xy2vms_pack(wf, sender)
!!!      call xy2vms_transpose(sender, recver)
!!!      call xy2vms_unpack(recver, cf)
!!!
!!!      deallocate( sender )
!!!      deallocate( recver )
!!!      deallocate( wf )
!!!      deallocate( w1 )
!!!      deallocate( w2 )
!!!
!!!  END SUBROUTINE colliimp_calc_colli_full


!--------------------------------------
  SUBROUTINE colliimp_calc_colli_full( ff, phi, cf )
!--------------------------------------
!   Collsion operator calculation interface

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)             :: phi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf

    complex(kind=DP), dimension(:,:,:,:), allocatable ::  &
           send1e, recv1e, send2e, recv2e, send1o, recv1o, send2o, recv2o
    complex(kind=DP), dimension(:,:,:,:), allocatable :: wh1e, wh2e, wh1o, wh2o
    integer :: mx, my, iz, iv, im, is, ibuff, iproc

      allocate( send1e(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( recv1e(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( send2e(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( recv2e(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( wh1e(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) )
      allocate( wh2e(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) )
      allocate( send1o(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( recv1o(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( send2o(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( recv2o(1:2*nv,0:nm,0:nbuff-1,0:nprocvms-1) )
      allocate( wh1o(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) )
      allocate( wh2o(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) )

!$OMP parallel default(none) &
!$OMP shared(ff,phi,cf,ist_y,iend_y,sgn,Znum,fmx,tau,j0,iFLR,ranks) &
!$OMP shared(send1e,recv1e,send2e,recv2e,wh1e,wh2e) &
!$OMP shared(send1o,recv1o,send2o,recv2o,wh1o,wh2o) &
!$OMP private(mx,my,iz,iv,im,is,ibuff,iproc)
      do iproc = 0, nprocvms-1
!$OMP do
        do ibuff = 0, nbuff-1
          do im = 0, nm
            do iv = 1, 2*nv
              send1e(iv,im,ibuff,iproc) = (0._DP, 0._DP)
              send2e(iv,im,ibuff,iproc) = (0._DP, 0._DP)
              send1o(iv,im,ibuff,iproc) = (0._DP, 0._DP)
              send2o(iv,im,ibuff,iproc) = (0._DP, 0._DP)
            end do
          end do
        end do
!$OMP end do nowait
      end do
      do ibuff = 0, nbuff-1
!$OMP do
        do is = 0, ns-1
          do im = 0, global_nm
            do iv = 1, 2*global_nv
              wh1e(iv,im,is,ibuff) = (0._DP, 0._DP)
              wh1o(iv,im,is,ibuff) = (0._DP, 0._DP)
            end do
          end do
        end do
!$OMP end do nowait
      end do

      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = iend_y, ny
              cf(:,my,iz,iv,im) = (0._DP, 0._DP)
            end do
            do my = ist_y, iend_y
              do mx = -nx, nx
                cf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) + sgn(ranks) * Znum(ranks) &
                          * fmx(iz,iv,im) / tau(ranks) * j0(mx,my,iz,im) * phi(mx,my,iz) &
                                                                     * real(iFLR, kind=DP)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP barrier

!!%%% Without overlap %%%
!      do iz = -nz, nz-1
!        call vms2xy_pack_iz(iz, cf, send1e)
!!$OMP barrier
!!$OMP master
!        call vms2xy_transpose_iz(send1e, recv1e)
!!$OMP end master
!!$OMP barrier
!        call vms2xy_unpack_iz(recv1e, wh1e)
!!$OMP barrier
!        call collision_full_wrapper_iz(iz, wh1e, wh2e)
!!$OMP barrier
!        call xy2vms_pack_iz(wh2e, send2e)
!!$OMP barrier
!!$OMP master
!        call xy2vms_transpose_iz(send2e, recv2e)
!!$OMP end master
!!$OMP barrier
!        call xy2vms_unpack_iz(iz, recv2e, cf)
!!$OMP barrier
!      end do
!!%%%%%%%%%%%%%%%%%%%%%%%

!%%% With overlap %%%
      do iz = -nz, nz-1+6

        if (mod(iz+nz,2) == 0) then ! even

!$OMP master
          if (-nz+1<=iz .and. iz<=nz-1+1) call vms2xy_transpose_iz(send1o, recv1o)
          if (-nz+5<=iz .and. iz<=nz-1+5) call xy2vms_transpose_iz(send2o, recv2o)
!$OMP end master
          if (-nz  <=iz .and. iz<=nz-1  ) call vms2xy_pack_iz(iz, cf, send1e)
          if (-nz+2<=iz .and. iz<=nz-1+2) call vms2xy_unpack_iz(recv1e, wh1e)
          if (-nz+3<=iz .and. iz<=nz-1+3) call collision_full_wrapper_iz(iz-3, wh1o, wh2o)
          if (-nz+4<=iz .and. iz<=nz-1+4) call xy2vms_pack_iz(wh2e, send2e)
          if (-nz+6<=iz .and. iz<=nz-1+6) call xy2vms_unpack_iz(iz-6, recv2e, cf)

        else                        ! odd

!$OMP master
          if (-nz+1<=iz .and. iz<=nz-1+1) call vms2xy_transpose_iz(send1e, recv1e)
          if (-nz+5<=iz .and. iz<=nz-1+5) call xy2vms_transpose_iz(send2e, recv2e)
!$OMP end master
          if (-nz  <=iz .and. iz<=nz-1  ) call vms2xy_pack_iz(iz, cf, send1o)
          if (-nz+2<=iz .and. iz<=nz-1+2) call vms2xy_unpack_iz(recv1o, wh1o)
          if (-nz+3<=iz .and. iz<=nz-1+3) call collision_full_wrapper_iz(iz-3, wh1e, wh2e)
          if (-nz+4<=iz .and. iz<=nz-1+4) call xy2vms_pack_iz(wh2o, send2o)
          if (-nz+6<=iz .and. iz<=nz-1+6) call xy2vms_unpack_iz(iz-6, recv2o, cf)

        end if
!$OMP barrier

      end do
!%%%%%%%%%%%%%%%%%%%%

!$OMP end parallel

      deallocate( send1e )
      deallocate( recv1e )
      deallocate( send2e )
      deallocate( recv2e )
      deallocate( wh1e )
      deallocate( wh2e )
      deallocate( send1o )
      deallocate( recv1o )
      deallocate( send2o )
      deallocate( recv2o )
      deallocate( wh1o )
      deallocate( wh2o )

  END SUBROUTINE colliimp_calc_colli_full


!--------------------------------------
  SUBROUTINE collision_full_wrapper_iz( iz, wfin, wfout )
!--------------------------------------
!   Collsion operator calculation interface

    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) :: wfin
    complex(kind=DP), intent(out), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1,0:nbuff-1) :: wfout

    complex(kind=DP), &
      dimension(1-nvb:2*global_nv+nvb,0-nvb:global_nm+nvb,0:ns-1) :: w1
    complex(kind=DP), &
      dimension(1:2*global_nv,0:global_nm,0:ns-1) :: w2
    integer :: ibuff, ivb

!$OMP master
                                      call clock_sta(1720)
!$OMP end master
      w1(:,:,:) = (0._DP, 0._DP)
!$OMP do schedule(dynamic)
      do ibuff = 0, nbuff-1
          w1(1:2*global_nv,0:global_nm,0:ns-1) = wfin(:,:,:,ibuff)
         !- Boundary condition -
          do ivb = 1, nvb
           !w1(1-ivb,:,:) = (0._DP, 0._DP)
           !w1(2*global_nv+ivb,:,:) = (0._DP, 0._DP)
           !w1(:,global_nm+ivb,:) = (0._DP, 0._DP)
            w1(:,-ivb,:) = w1(:,ivb,:)
          end do
         !-
          call collision_full(w1, w2, iz, ibuff)
          wfout(:,:,:,ibuff) = w2(:,:,:)
      end do
!$OMP end do nowait
!$OMP master
                                      call clock_end(1720)
!$OMP end master

  END SUBROUTINE collision_full_wrapper_iz


END MODULE GKV_colliimp
