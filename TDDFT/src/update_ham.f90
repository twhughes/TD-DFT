!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE update_hamiltonian(istepp)
  !-----------------------------------------------------------------------
  !
  ! ... Update the hamiltonian
  !
  USE kinds,         ONLY : dp
  USE ldaU,          ONLY : lda_plus_U
  USE scf,           ONLY : rho, rho_core, rhog_core, vltot, v, kedtau, vrs
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE io_global,     ONLY : stdout, ionode
  USE lsda_mod,      ONLY : nspin
  USE uspp,          ONLY : okvan, nkb
  USE dfunct,        ONLY : newd
  USE tddft_module   
  USE becmod,        ONLY : becp, allocate_bec_type, deallocate_bec_type, &
                            is_allocated_bec_type
  USE wvfct,         ONLY : nbnd
  USE ldaU,          ONLY : eth
  USE ener,          ONLY : etot, eband, etxc, ehart, etxcc, vtxc, ewld  
  USE cell_base,     ONLY : at, bg, alat, omega
  USE ions_base,     ONLY : zv, nat, nsp, ityp, tau
  USE gvect,         ONLY : ngm, gstart, g, gg, gcutm
  USE control_flags, ONLY : gamma_only
  USE vlocal,        ONLY : strf
  USE fde
  use mp_images, only : my_image_id
! modules needed for adding the add_efield_t
  USE noncollin_module, ONLY :  nspin_lsda
  implicit none
  interface
    subroutine calculate_eband(istepp)
      integer, intent(in), optional :: istepp
    end subroutine calculate_eband
  end interface

  integer, intent(in) :: istepp
  integer :: is
  real(dp) :: charge
  real(dp), external :: ewald

  call start_clock('updateH')
  
  if (is_allocated_bec_type(becp)) call deallocate_bec_type(becp)
  call sum_band_tddft(iunevcn,rho)
  if (do_fde) then
!      write(stdout,*)'coupled',coupled
!      call flush_unit(stdout)
      if (coupled) then
             call update_rho_fde(rho, .true.)
      else
             call update_rho_fde(rho, .false.)
      endif
  endif
           
  if (.not. is_allocated_bec_type(becp)) call allocate_bec_type(nkb, nbnd, becp)


  if (lda_plus_U) then
    call new_ns
    if (iverbosity > 20) call write_ns()
  end if
    

  ! calculate HXC-potential
  call v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, etotefield, charge, v )
   
! add local version of add_efield_t
! Caution: we're using our own tefield, which is different from the tefield in
! extfield module so when add_efield is called through v_of_rho, tefield is
! false even if it is true in tddft. Reason: first variable set to false and
! forces exit within add_efield.
   if (itdpulse>1.and.(dt*(istepp-1)<=ttend)) then
     DO is = 1, nspin_lsda
        CALL add_efield_t(v%of_r(1,is),   .true.)
     END DO
   endif 


  ! calculate total local potential (external + scf)
  call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    
 
  ! calculate eband and etot
 call calculate_eband()

  if (do_fde .and. .not. fde_init_rho) then
     ewld = ewald( alat, nat_fde, nsp, ityp_fde, zv, at, bg, tau_fde, &
                   omega, g, gg, ngm, gcutm, gstart, gamma_only, strf_fde )
  else
     ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                   omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  endif
! etot is calculated for each fragment separately, but the ehart is calculated
! using the total rho_fde. So in the uncoupled version, summing the energies up
! doesn't make sense because they rho_fde is NOT THE SAME on the different
! fragments and hence the total energy will be not correct. To fix this, we need
! to introduce a "total energy" for each fragment, where instead of summing the
! etot up over inter_fragment_comm tot etot_fde, the etot_fde of each fragment
! is partially updated with the etot of the fragment itself AND the ehart of the
! frozen fragments which also changes because it is calculated with the rho_fde.
  etot = eband + ( etxc - etxcc ) + ewld + ehart
  if (do_fde) call fde_energies()
!   if (ionode) then
!       if (do_fde) then
!          write(90+my_image_id,'(''ENERGY '',2X,I6,6E16.6)') istep, etot, eband, ehart, etxc, etxcc
!       else
!          write(90+my_image_id,'(''ENERGY '',2X,I6,5E16.6)') istep, etot, eband, ehart, etxc, etxcc
!       endif
!    endif

  ! calculate new D_nm matrix for ultrasoft pseudopotential
  if (okvan) then
    if (istepp == -1 .or. ( (nupdate_Dnm /= 0 .and. mod(istepp,nupdate_Dnm) == 0) ) ) then
      call newd()
!      if (iverbosity > 10) write(stdout,'(5X,''call newd'')')
    endif
  endif
    
  call deallocate_bec_type(becp)
  call stop_clock('updateH')
    
END SUBROUTINE update_hamiltonian

 !-----------------------------------------------------------------------
 SUBROUTINE calculate_eband(istepp)
  !-----------------------------------------------------------------------
  !
  USE kinds,                       ONLY : dp
!  USE wavefunctions_module,        ONLY : evc
  USE io_global,                   ONLY : stdout
  USE io_files,                    ONLY : nwordwfc, iunwfc, iunigk
  USE klist,                       ONLY : nks, wk
  USE wvfct,                       ONLY : nbnd, npwx, npw, wg, current_k, et, igk, g2kin
  USE lsda_mod,                    ONLY : current_spin, isk, nspin
  USE becmod,                      ONLY : becp, calbec, allocate_bec_type, &
                                          is_allocated_bec_type, deallocate_bec_type
  USE mp,                          ONLY : mp_sum, mp_barrier
  USE mp_pools,                    ONLY : intra_pool_comm, inter_pool_comm
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE uspp,                        ONLY : nkb, vkb
  USE pwcom
  USE tddft_module
  use constants,                   only : rytoev
  ! ... Calculate band energy
  !
  implicit none
  complex(dp), allocatable :: hpsi(:,:)
  integer :: ik, ibnd
  complex(dp), external :: zdotc
  integer, intent(in), optional :: istepp
  complex(dp), allocatable :: evc_temp(:,:)
 
  if (.not. is_allocated_bec_type(becp)) call allocate_bec_type(nkb, nbnd, becp)
  allocate(hpsi(npwx,nbnd))
  allocate(evc_temp(npwx,nbnd))
  eband = 0.d0
  ! loop over k-points     
  if (nks > 1) rewind (iunigk)
  do ik = 1, nks
    current_k = ik
    current_spin = isk(ik)
   
!*******************************************************************
!! NEW PART, taken from PW/src
!    npw = ngk(ik)
!     IF (nks > 1) THEN
!        READ (iunigk) igk
!     endif
!    
!    ! read wfcs from file and compute becp
!    evc = (0.d0, 0.d0)
!    if (present(istepp)) then 
!        call get_buffer (evc, nwordwfc, iunwfc, ik)
!    else
!        call get_buffer(evc, nwordwfc, iunevcn, ik)
!    endif
!    call init_us_2(npw, igk, xk(1,ik), vkb)

!*******************************************************************
! OLD PART, CHECK IF OK: IT IS NOT
!    ! initialize at k-point k 
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin = g2kin * tpiba2
    call init_us_2(npw, igk, xk(1,ik), vkb)
!
    ! read wfcs from file and compute becp
    evc_temp = (0.d0, 0.d0)
    if (present(istepp)) then 
        call get_buffer (evc_temp, nwordwfc, iunwfc, ik)
    else
        call get_buffer(evc_temp, nwordwfc, iunevcn, ik)
    endif
!
!*******************************************************************
    call calbec(npw, vkb, evc_temp, becp)

 
    ! calculate H |psi_current>
!    call h_psi(npwx, npw, nbnd_occ(ik), evc_temp, hpsi)
    call h_psi(npwx, npw, nbnd, evc_temp, hpsi)
    
!    ! add to eband
!    do ibnd = 1, nbnd_occ(ik)
  do ibnd = 1, nbnd
        et(ibnd,ik) = real(zdotc(npw,hpsi(1,ibnd),1,evc_temp(1,ibnd),1),dp)
    enddo
    do ibnd = 1, nbnd_occ(ik)
        eband = eband + wg(ibnd,ik)*real(zdotc(npw,hpsi(1,ibnd),1,evc_temp(1,ibnd),1),dp)
    enddo
  enddo
  call mp_sum(et, intra_pool_comm)
  call mp_sum(et, inter_pool_comm)
  call mp_sum(eband, intra_pool_comm)
  call mp_sum(eband, inter_pool_comm)
  
  eband = eband + delta_e()
  
  deallocate(hpsi)
  deallocate(evc_temp)
  call deallocate_bec_type(becp)
  
  CONTAINS
  
   !-----------------------------------------------------------------------
   FUNCTION delta_e()
     !-----------------------------------------------------------------------
     ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
     !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
     !               - \sum rho%ns       v%ns       [for LDA+U]
     !               - \sum becsum       D1_Hxc     [for PAW]
     USE scf,              ONLY : scf_type, rho, v
     USE funct,            ONLY : dft_is_meta
     USE fft_base,         ONLY : dfftp
     USE noncollin_module, ONLY : noncolin
     USE mp_bands,         ONLY : intra_bgrp_comm
     USE paw_variables,    ONLY : okpaw, ddd_paw
     IMPLICIT NONE
     REAL(DP) :: delta_e, delta_e_hub
     !
     delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
     !
     IF ( dft_is_meta() ) &
        delta_e = delta_e - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )
     !
     delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
     !
     CALL mp_sum( delta_e, intra_bgrp_comm )
     !
     if (lda_plus_u) then
       if (noncolin) then
         delta_e_hub = - real(SUM (rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:)),kind=dp)
         delta_e = delta_e + delta_e_hub
       else
         delta_e_hub = - SUM (rho%ns(:,:,:,:)*v%ns(:,:,:,:))
         if (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
         delta_e = delta_e + delta_e_hub
       endif
     end if
     !
     IF (okpaw) delta_e = delta_e - SUM(ddd_paw(:,:,:)*rho%bec(:,:,:))
     !
     RETURN
     !
   END FUNCTION delta_e

 END SUBROUTINE calculate_eband

   
   
 !-----------------------------------------------------------------------
 SUBROUTINE calculate_eband_deltat(rho_deltat,v_deltat)
  !-----------------------------------------------------------------------
  !
  USE kinds,                       ONLY : dp
!  USE wavefunctions_module,        ONLY : evc
  USE io_global,                   ONLY : stdout
  USE io_files,                    ONLY : nwordwfc, iunwfc, iunigk
  USE klist,                       ONLY : nks
  USE wvfct,                       ONLY : nbnd, npwx, npw, wg, current_k, et, igk, g2kin
  USE lsda_mod,                    ONLY : current_spin, isk, nspin
  USE becmod,                      ONLY : becp, calbec, allocate_bec_type, &
                                          is_allocated_bec_type
  USE mp,                          ONLY : mp_sum, mp_barrier
  USE mp_pools,                    ONLY : intra_pool_comm, inter_pool_comm
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE uspp,                        ONLY : nkb, vkb
  USE pwcom
  USE tddft_module
  use scf, only : scf_type
  ! ... Calculate band energy
  !
  implicit none
  complex(dp), allocatable :: hpsi(:,:)
  complex(dp), allocatable :: tddft_psi_deltat(:,:,:)
  integer :: ik, ibnd
  complex(dp), external :: zdotc
  type (scf_type) :: rho_deltat, v_deltat
 
  if (.not. is_allocated_bec_type(becp)) call allocate_bec_type(nkb, nbnd, becp)
  allocate(hpsi(npwx,nbnd))
   allocate (tddft_psi_deltat(npwx,nbnd,2))
  eband = 0.d0
  ! loop over k-points     
  if (nks > 1) rewind (iunigk)
  do ik = 1, nks
    current_k = ik
    current_spin = isk(ik)
   
!*******************************************************************
!! NEW PART, taken from PW/src
!    npw = ngk(ik)
!     IF (nks > 1) THEN
!        READ (iunigk) igk
!     endif
!    
!    ! read wfcs from file and compute becp
!    evc = (0.d0, 0.d0)
!    if (present(istepp)) then 
!        call get_buffer (evc, nwordwfc, iunwfc, ik)
!    else
!        call get_buffer(evc, nwordwfc, iunevcn, ik)
!    endif
!    call init_us_2(npw, igk, xk(1,ik), vkb)

!*******************************************************************
! OLD PART, CHECK IF OK: IT IS NOT
!    ! initialize at k-point k 
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin = g2kin * tpiba2
    call init_us_2(npw, igk, xk(1,ik), vkb)
!
    ! read wfcs from file and compute becp
        call get_buffer(tddft_psi_deltat, 2*nwordwfc, iuntdwfc, ik)
!
!*******************************************************************
    call calbec(npw, vkb, tddft_psi_deltat(:,:,2), becp)

 
    ! calculate H |psi_current>
!    call h_psi(npwx, npw, nbnd_occ(ik), evc, hpsi)
    call h_psi(npwx, npw, nbnd, tddft_psi_deltat(:,:,2), hpsi)
    
    ! add to eband
!    do ibnd = 1, nbnd_occ(ik)
    do ibnd = 1, nbnd
        et(ibnd,ik) = real(zdotc(npw,hpsi(1,ibnd),1,tddft_psi_deltat(1,ibnd,2),1),dp)
    enddo
    do ibnd = 1, nbnd_occ(ik)
        eband = eband + wg(ibnd,ik)*real(zdotc(npw,hpsi(1,ibnd),1,tddft_psi_deltat(1,ibnd,2),1),dp)
    enddo
  enddo
  call mp_sum(et, intra_pool_comm)
  call mp_sum(et, inter_pool_comm)
  call mp_sum(eband, intra_pool_comm)
  call mp_sum(eband, inter_pool_comm)
  
  eband = eband + delta_e_deltat(rho_deltat,v_deltat)
  
  deallocate(hpsi)
  deallocate(tddft_psi_deltat)
  
  CONTAINS
  
   !-----------------------------------------------------------------------
   FUNCTION delta_e_deltat(rho,v)
     !-----------------------------------------------------------------------
     ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
     !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
     !               - \sum rho%ns       v%ns       [for LDA+U]
     !               - \sum becsum       D1_Hxc     [for PAW]
     USE scf,              ONLY : scf_type
     USE funct,            ONLY : dft_is_meta
     USE fft_base,         ONLY : dfftp
     USE noncollin_module, ONLY : noncolin
     USE mp_bands,         ONLY : intra_bgrp_comm
     USE paw_variables,    ONLY : okpaw, ddd_paw
     IMPLICIT NONE
     REAL(DP) :: delta_e_deltat, delta_e_hub
     type (scf_type), intent(inout) :: rho, v
     !
     delta_e_deltat = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
     !
     IF ( dft_is_meta() ) &
        delta_e_deltat = delta_e_deltat - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )
     !
     delta_e_deltat = omega * delta_e_deltat / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
     !
     CALL mp_sum( delta_e_deltat, intra_bgrp_comm )
     !
     if (lda_plus_u) then
       if (noncolin) then
         delta_e_hub = - real(SUM (rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:)),kind=dp)
         delta_e_deltat = delta_e_deltat + delta_e_hub
       else
         delta_e_hub = - SUM (rho%ns(:,:,:,:)*v%ns(:,:,:,:))
         if (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
         delta_e_deltat = delta_e_deltat + delta_e_hub
       endif
     end if
     !
     IF (okpaw) delta_e_deltat = delta_e_deltat - SUM(ddd_paw(:,:,:)*rho%bec(:,:,:))
     !
     RETURN
     !
   END FUNCTION delta_e_deltat

 END SUBROUTINE calculate_eband_deltat
