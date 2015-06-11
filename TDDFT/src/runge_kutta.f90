subroutine runge_kutta()

  USE tddft_module
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, g2kin, current_k, ecutwfc
  USE wavefunctions_module,        ONLY : evc
  USE becmod,                      ONLY : becp, allocate_bec_type, is_allocated_bec_type, deallocate_bec_type
  USE uspp,                        ONLY : nkb, vkb
  USE lsda_mod,                    ONLY : current_spin, isk, nspin
  USE gvect,                       ONLY : ngm, g
  USE io_files,                    ONLY : nwordwfc, iunigk
  USE klist,                       ONLY : nks, xk
!  USE io_global,                   ONLY : stdout, ionode
  USE cell_base,                   ONLY : tpiba2
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE scf,           ONLY : rho, rho_core, rhog_core, vltot, v, kedtau, vrs, scf_type_copy
  USE fft_base,      ONLY : dfftp
  use ruku
  use fde, only : do_fde, rho_old 
  USE gvecs,         ONLY : doublegrid
  USE ener,          ONLY : etxc, ehart, vtxc
  USE ldaU,          ONLY : eth
  USE noncollin_module, ONLY : nspin_lsda

  implicit none
  include 'mpif.h'

  integer     :: ik, is
  real(dp) :: charge
  call start_clock('runge_kutta')

!*******************************************************************
irkloop: do irk4 = 1,4
!*******************************************************************
! for each RK step, have a loop over k-points, calculating the RK vectors
!*******************************************************************
if ((do_fde).and.(.not.update)) call scf_type_copy(rho, rho_old)
!   ! loop over k-points     
    if (nks > 1) rewind (iunigk)
    do ik = 1, nks
!*******************************************************************
      current_k = ik
      current_spin = isk(ik)


!      OLD PART
!      ! initialize at k-point k 
      call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      g2kin = g2kin * tpiba2
      call init_us_2(npw, igk, xk(1,ik), vkb)

      call get_buffer (evc, nwordwfc, iunevcn, ik)
      if (irk4==1) then
      call get_buffer (tddft_psi, nwordtdwfc*2, iuntdwfc, ik)
      tddft_psi(:,:,2)=evc(:,:) ! save the original psi
      call save_buffer (tddft_psi, nwordtdwfc*2, iuntdwfc, ik)
      endif
!*******************************************************************
        
     ! calculate H |psi_current>, S^-1 |Hpsi_current>
     if (.not. is_allocated_bec_type(becp)) call allocate_bec_type(nkb, nbnd, becp)
     call h_psi(npwx, npw, nbnd_occ(ik), evc, tddft_hpsi)
     if (ehrenfest .and. istep > 1) then
        call apply_capital_P(npwx, npw, nbnd_occ(ik), evc, tddft_Ppsi)
        tddft_hpsi = tddft_hpsi + tddft_Ppsi
     endif
     call sinv_hpsi(npwx, nbnd_occ(ik), tddft_hpsi, tddft_spsi)
    
     if (irk4==1) then
          k1(1:npw,1:nbnd_occ(ik))= - 2.0d0 * ee * tddft_spsi(1:npw,1:nbnd_occ(ik)) 
         call save_buffer(k1,nwordtdwfc,iunruku1,ik)
     endif
     if (irk4==2) then
          k1(1:npw,1:nbnd_occ(ik))= - 2.0d0 * ee * tddft_spsi(1:npw,1:nbnd_occ(ik)) 
         call save_buffer(k1,nwordtdwfc,iunruku2,ik)
     endif
     if (irk4==3) then
          k1(1:npw,1:nbnd_occ(ik))= - 2.0d0 * ee * tddft_spsi(1:npw,1:nbnd_occ(ik)) 
         call save_buffer(k1,nwordtdwfc,iunruku3,ik)
     endif
     if (irk4==4) then
          k1(1:npw,1:nbnd_occ(ik))= - 2.0d0 * ee * tddft_spsi(1:npw,1:nbnd_occ(ik)) 
         call save_buffer(k1,nwordtdwfc,iunruku4,ik)
     endif
     
     
!*******************************************************************
  enddo ! ik
    
!*******************************************************************
! update psi with the k1 vector
if (irk4==1) then
  do ik = 1, nks
      call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      g2kin = g2kin * tpiba2
      call init_us_2(npw, igk, xk(1,ik), vkb)
      CALL start_clock('stop_k')
    call get_buffer(tddft_psi,nwordtdwfc*2,iuntdwfc,ik)
    call get_buffer(k1,nwordtdwfc,iunruku1,ik)
    tddft_psi(1:npw,1:nbnd_occ(ik),1)=tddft_psi(1:npw,1:nbnd_occ(ik),2)+0.5d0*k1(1:npw,1:nbnd_occ(ik))
    call save_buffer(tddft_psi(:,:,1),nwordwfc, iunevcn, ik)
  enddo
!  call update_hamiltonian(istep)
  call sum_band_tddft(iunevcn,rho)
  if (do_fde) then
     if (update) then
       call update_rho_fde(rho, .true.) 
    else 
       call update_rho_fde(rho, .false.) 
    endif
  endif
  call v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, .false., charge, v )
   if (itdpulse>1.and.(dt*(istep-1)<=ttend)) then
     DO is = 1, nspin_lsda
        CALL add_efield_t(v%of_r(1,is),  .false.)
     END DO
   endif 
  call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    
elseif (irk4==2) then
  do ik = 1, nks
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin = g2kin * tpiba2
    call init_us_2(npw, igk, xk(1,ik), vkb)
    call get_buffer(tddft_psi,nwordtdwfc*2,iuntdwfc,ik)
    call get_buffer(k1,nwordtdwfc,iunruku2,ik)
    if (rkmax==2) then 
       tddft_psi(1:npw,1:nbnd_occ(ik),1)=tddft_psi(1:npw,1:nbnd_occ(ik),2)+k1(1:npw,1:nbnd_occ(ik))
      tddft_psi(:,:,2) = tddft_psi(:,:,1)
      evc(:,:) = tddft_psi(:,:,1)
      call save_buffer(evc,nwordwfc, iunevcn, ik)
      call save_buffer (tddft_psi, nwordtdwfc*2, iuntdwfc, ik)
    else
      tddft_psi(1:npw,1:nbnd_occ(ik),1)=tddft_psi(1:npw,1:nbnd_occ(ik),2)+0.5d0*k1(1:npw,1:nbnd_occ(ik))
      call save_buffer(tddft_psi(:,:,1),nwordwfc, iunevcn, ik)
     endif
  enddo
  if (rkmax==2) exit irkloop
  call sum_band_tddft(iunevcn,rho)
  if (do_fde) then
     if (update) then
       call update_rho_fde(rho, .true.) 
    else 
       call update_rho_fde(rho, .false.) 
    endif
  endif
  call v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, .false., charge, v )
   if (itdpulse>1.and.(dt*(istep-1)<=ttend)) then
     DO is = 1, nspin_lsda
        CALL add_efield_t(v%of_r(1,is),   .false.)
     END DO
   endif 
  call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    
elseif (irk4==3) then
  do ik = 1, nks
      call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      g2kin = g2kin * tpiba2
      call init_us_2(npw, igk, xk(1,ik), vkb)
    call get_buffer(tddft_psi,nwordtdwfc*2,iuntdwfc,ik)
    call get_buffer(k1,nwordtdwfc,iunruku3,ik)
    tddft_psi(1:npw,1:nbnd_occ(ik),1)=tddft_psi(1:npw,1:nbnd_occ(ik),2)+      k1(1:npw,1:nbnd_occ(ik))
    call save_buffer(tddft_psi(:,:,1),nwordwfc, iunevcn, ik)
  enddo
  call sum_band_tddft(iunevcn,rho)
  if (do_fde) then
     if (update) then
       call update_rho_fde(rho, .true.) 
    else 
       call update_rho_fde(rho, .false.) 
    endif
  endif
  call v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, .false., charge, v )
   if (itdpulse>1.and.(dt*(istep-1)<=ttend)) then
     DO is = 1, nspin_lsda
        CALL add_efield_t(v%of_r(1,is),   .false.)
     END DO
   endif 
  call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    
elseif (irk4==4) then
  do ik = 1, nks
      call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      g2kin = g2kin * tpiba2
      call init_us_2(npw, igk, xk(1,ik), vkb)
    call get_buffer(tddft_psi,nwordtdwfc*2,iuntdwfc,ik)
    call get_buffer(k1,nwordtdwfc,iunruku1,ik)
    tddft_psi(1:npw,1:nbnd_occ(ik),1)=tddft_psi(1:npw,1:nbnd_occ(ik),2)+1.0d0/6.0d0*k1(1:npw,1:nbnd_occ(ik)) 
    call get_buffer(k1,nwordtdwfc,iunruku2,ik)
    tddft_psi(1:npw,1:nbnd_occ(ik),1)=tddft_psi(1:npw,1:nbnd_occ(ik),1) +1.0d0/3.0d0*k1(1:npw,1:nbnd_occ(ik)) 
    call get_buffer(k1,nwordtdwfc,iunruku3,ik)
    tddft_psi(1:npw,1:nbnd_occ(ik),1)=tddft_psi(1:npw,1:nbnd_occ(ik),1) +1.0d0/3.0d0*k1(1:npw,1:nbnd_occ(ik)) 
    call get_buffer(k1,nwordtdwfc,iunruku4,ik)
    tddft_psi(1:npw,1:nbnd_occ(ik),1)=tddft_psi(1:npw,1:nbnd_occ(ik),1)+1.0d0/6.0d0*k1(1:npw,1:nbnd_occ(ik)) 
    ! update the wavefunctions
    tddft_psi(:,:,2) = tddft_psi(:,:,1)
    evc(:,:) = tddft_psi(:,:,1)
    call save_buffer(evc,nwordwfc, iunevcn, ik)
   call save_buffer (tddft_psi, nwordtdwfc*2, iuntdwfc, ik)
   enddo
 endif
enddo irkloop


  call stop_clock('runge_kutta')
return
end subroutine runge_kutta

