subroutine crank_nicolson

  USE tddft_module
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, g2kin, current_k, ecutwfc
  USE wavefunctions_module,        ONLY : evc
  USE becmod,                      ONLY : becp, allocate_bec_type, is_allocated_bec_type, deallocate_bec_type
  USE uspp,                        ONLY : nkb, vkb
  USE lsda_mod,                    ONLY : current_spin, isk, nspin
  USE gvect,                       ONLY : ngm, g
  USE io_files,                    ONLY : nwordwfc, iunwfc, iunigk
  USE klist,                       ONLY : nks, xk
  USE io_global,                   ONLY : stdout, ionode, ionode_id
  USE cell_base,                   ONLY : tpiba2
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE scf,                         only : rho_core, rhog_core, kedtau, scf_type_copy, vltot, rho
  USE fft_base,   ONLY : dfftp
  USE ener,          ONLY : etxc, ehart, vtxc, ehart
  USE ldaU,          ONLY : eth
  use fde, only : do_fde, rho_old, nfragments
  USE noncollin_module, ONLY : nspin_lsda
  USE gvecs,         ONLY : doublegrid
  USE mp,                          ONLY : mp_sum, mp_bcast
  USE mp_pools,  ONLY : intra_pool_comm
  use cnmod
  USE mp_images,            ONLY : inter_fragment_comm, intra_image_comm

  implicit none
  include 'mpif.h'

  integer     :: ik, ibnd, lter, flag_global, ip, is, fde_cycle
  external tddft_ch_psi_all
  real(dp) :: charge

! For the predictor corrector scheme
  
call start_clock('crank_nicolson')
iploop: do ip = 1, pmax
    if (ip==1) then
!       if ((do_fde).and.(.not.update).and.(.not.coupled)) call scf_type_copy(rho, rho_old)
     continue
    else
        if ((do_fde).and.(.not.update)) call scf_type_copy(rho_deltat, rho_old)
    endif

  if (.not. is_allocated_bec_type(becp)) call allocate_bec_type(nkb, nbnd, becp)
!*******************************************************************
    ! loop over k-points     
    if (nks > 1) rewind (iunigk)
    do ik = 1, nks
!*******************************************************************
      current_k = ik
      current_spin = isk(ik)
!*******************************************************************
!      OLD PART
!      ! initialize at k-point k 
      call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      g2kin = g2kin * tpiba2
      call init_us_2(npw, igk, xk(1,ik), vkb)

      ! read wfcs from file and compute becp
      if ((istep==1).and.(itefield==1).and.(.not.l_tddft_restart).and.(ip==1)) then
        call get_buffer (evc, nwordwfc, iunwfc, ik)
      else
        call get_buffer (evc, nwordwfc, iunevcn, ik) 
      endif
      call get_buffer (tddft_psi, nwordtdwfc*2, iuntdwfc, ik) 
      
!*******************************************************************
        
     ! guess the wavefunction at the next timestep
     if (ip==1.and.istep==1) then ! fix later
!        do ibnd = 1, nbnd_occ(ik) 
        do ibnd = 1, nbnd
          tddft_psi(:,ibnd,2) = 2.d0*evc(:,ibnd) - tddft_psi(:,ibnd,2)
        enddo
     endif

     ! calculate H |psi_current>, S |psi_current>
     if (ip ==1 ) then
!         call h_psi(npwx, npw, nbnd_occ(ik), evc, tddft_hpsi)
         call h_psi(npwx, npw, nbnd, evc, tddft_hpsi)
     else
!         call h_psi       (npwx, npw, nbnd_occ(ik),                       evc, tddft_hpsi)
!         call h_psi_deltat(npwx, npw, nbnd_occ(ik), vrs_deltat, v_deltat, evc, tddft_hpsi_deltat)
         call h_psi       (npwx, npw, nbnd,                       evc, tddft_hpsi)
         call h_psi_deltat(npwx, npw, nbnd, vrs_deltat, v_deltat, evc, tddft_hpsi_deltat)
         tddft_hpsi=(tddft_hpsi+tddft_hpsi_deltat)*0.5d0
     endif
         
     if (ehrenfest .and. istep > 1) then
!!        call apply_capital_P(npwx, npw, nbnd_occ(ik), evc, tddft_Ppsi)
        call apply_capital_P(npwx, npw, nbnd, evc, tddft_Ppsi)
        tddft_hpsi = tddft_hpsi + tddft_Ppsi
     endif
!     call s_psi(npwx, npw, nbnd_occ(ik), evc, tddft_spsi)
     call s_psi(npwx, npw, nbnd, evc, tddft_spsi)
  
     ! calculate (S - H*dt*i/2) |\psi_current>
     b = (0.d0, 0.d0)
!!     b(1:npw, 1:nbnd_occ(ik)) = tddft_spsi(1:npw,1:nbnd_occ(ik)) - ee * tddft_hpsi(1:npw,1:nbnd_occ(ik))
     b(1:npw, 1:nbnd) = tddft_spsi(1:npw,1:nbnd) - ee * tddft_hpsi(1:npw,1:nbnd)
  
     ! solve A * x = b
      call flush_unit(stdout)
!!     call tddft_cgsolver(tddft_ch_psi_all, b, tddft_psi(:,:,2), npwx, npw, &
 !!                   conv_threshold, lter, flag_global, nbnd_occ(ik), ee, ip, vrs_deltat, v_deltat)
     call tddft_cgsolver(tddft_ch_psi_all, b, tddft_psi(:,:,2), npwx, npw, &
                    conv_threshold, lter, flag_global, nbnd, ee, ip, vrs_deltat, v_deltat)
     if (iverbosity > 10.and.ionode) write(stdout,'(5x,a,i4)')'Number of CGS iterations=',lter
      call flush_unit(stdout)
     tddft_psi(:,:,1)=tddft_psi(:,:,2)
     call save_buffer (tddft_psi, nwordtdwfc*2, iuntdwfc, ik)
      ! update the wavefunctions
     if (pmax == 1) then
!  So tddft_psi(:,:,1) has the psi(t), tddft_psi(:,:,2) the tddft_psi(:,:,2) and evc the psi(t+deltat)
!!      evc(:,1:nbnd_occ(ik)) = tddft_psi(:,1:nbnd_occ(ik),2)
      evc(:,1:nbnd) = tddft_psi(:,1:nbnd,2)

     ! save wavefunctions to disk
      call save_buffer (evc, nwordwfc, iunevcn, ik)
     endif
        
!*******************************************************************
  enddo ! ik
!*******************************************************************
     call deallocate_bec_type(becp)
!*******************************************************************
  ! set up the density and potential at the t+deltat for the next predictor-corrector iteration
!----------------------
 if (pmax > 1) then
!----------------------

   if (ip>1) call scf_type_copy(rho_deltat, rho)

   call sum_band_tddft(iuntdwfc,rho_deltat)
   if (do_fde) then
      if (update) then
        call update_rho_fde(rho_deltat, .true.) 
     else 
        call update_rho_fde(rho_deltat, .false.) 
     endif
   endif
! test whether we are converged
   if (ip > 1) then
     rho%of_r = rho_deltat%of_r-rho%of_r
     norm_rho = dnrm2(dfftp%nnr,rho%of_r,1)
     call mp_sum(norm_rho,intra_pool_comm)
     if (.not. do_fde) then 
        if ((norm_rho <= pc_conv_threshold).or.(ip==pmax)) then
          do ik=1,nks
             call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
             g2kin = g2kin * tpiba2
             call init_us_2(npw, igk, xk(1,ik), vkb)
             call get_buffer(tddft_psi,nwordtdwfc*2,iuntdwfc,ik)
!             evc(:,1:nbnd_occ(ik)) = tddft_psi(:,1:nbnd_occ(ik),2)
             evc(:,1:nbnd) = tddft_psi(:,1:nbnd,2)
     !     save wavefunctions to disk
             call save_buffer (evc, nwordwfc, iunevcn, ik)
          enddo
        if (iverbosity > 10.and.ionode) write(stdout,'(5x,a,i4)')'Number of Predictor-Corrector steps=',ip
        exit iploop
       endif ! norm_rho
     else ! do_fde
! A check if all images converged or we need to do another round
     fde_cycle = 0
     if (ionode) then
        if ((norm_rho <= pc_conv_threshold).or.(ip==pmax)) fde_cycle = 1
        call mp_sum(fde_cycle, inter_fragment_comm)
     endif
     call mp_bcast( fde_cycle, ionode_id, intra_image_comm )
     if (fde_cycle==nfragments) then
        do ik=1,nks
           call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
           g2kin = g2kin * tpiba2
           call init_us_2(npw, igk, xk(1,ik), vkb)
           call get_buffer(tddft_psi,nwordtdwfc*2,iuntdwfc,ik)
!           evc(:,1:nbnd_occ(ik)) = tddft_psi(:,1:nbnd_occ(ik),2)
           evc(:,1:nbnd) = tddft_psi(:,1:nbnd,2)
     !   save wavefunctions to disk
           call save_buffer (evc, nwordwfc, iunevcn, ik)
        enddo
        if (iverbosity > 10.and.ionode) write(stdout,'(5x,a,i4)')'Number of Predictor-Corrector steps=',ip
        exit iploop
     endif ! fde_cycle
    endif ! do_fde
   endif ! ip>1
   
! set potential going with the new psi(t+deltat)
  call v_of_rho( rho_deltat, rho_core, rhog_core, ehart, etxc, vtxc, eth, .false., charge, v_deltat )
   if (itdpulse>1.and.(dt*(istep-1)<=ttend)) then
     DO is = 1, nspin_lsda
        CALL add_efield_t(v_deltat%of_r(1,is),  .false.)
     END DO
   endif 
  call set_vrs(vrs_deltat, vltot, v_deltat%of_r, kedtau, v_deltat%kin_r, dfftp%nnr, nspin, doublegrid)    
!----------------------
 endif ! pmax > 1
!----------------------
end do iploop 
    
!*******************************************************************

call stop_clock('crank_nicolson')
return
end subroutine crank_nicolson
