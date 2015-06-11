!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE apply_electric_field_t(tddftpsi)
  !-----------------------------------------------------------------------
  !
  ! ... Apply an electric field impuse at t = 0, a homogeneous phase-shift
  ! ... to each band
  !
  USE kinds,        ONLY : dp
  USE mp,           ONLY : mp_sum
  USE fft_base,     ONLY : dffts
  USE cell_base,    ONLY : alat
  USE gvecs,        ONLY : nls
  USE wvfct,        ONLY : current_k, igk, npw, npwx, nbnd
  USE io_files,     ONLY : nwordwfc
  USE buffers,      ONLY : save_buffer
  USE wavefunctions_module, ONLY : evc
  USE fft_base,       ONLY : dffts
  USE fft_interfaces, ONLY : invfft, fwfft
  USE tddft_module
  use constants,      only : au_sec, pi
  USE io_global,                   ONLY : stdout, ionode
  implicit none

  complex(dp), intent(out) :: tddftpsi(npwx,nbnd)
  integer :: ik, ibnd, ir
  complex(dp) :: psic(dffts%nnr)
  real(dp) :: phase 
  double precision ::  es, es2, t
  
  ik = current_k

  t = dt*(istep-1)
  if (itdpulse==2) then
    es = e_strength * cos(pi*(t-2.0d0*tau0-t0)/(2.0d0*tau0))
    es2 = real(es*exp(cmplx(0,omega0*t,kind=dp)),kind=dp)
    if (ionode) write(stdout,*)'pulse info',t,es,es2
  elseif (itdpulse==3) then
    es = e_strength * exp(-omega0*(t-t0)**2)
    if (ionode) write(stdout,*)'pulse info',t,es
  endif
  do ibnd = 1, nbnd_occ(ik)
    ! transform wavefunction from reciprocal space into real space
    psic = (0.d0, 0.d0)
    psic(nls(igk(1:npw))) = evc(1:npw, ibnd)
    call invfft ('Wave', psic, dffts)  

    do ir = 1, dffts%nnr
      phase = es2 * r_pos_s(e_direction,ir) * alat  !! CAREFUL WITH ULTRASOFT
      psic(ir) = psic(ir) * cmplx(cos(phase),sin(phase),dp)
    enddo

    call fwfft ('Wave', psic, dffts)  
       
    evc(1:npw,ibnd) = psic(nls(igk(1:npw)))
  enddo
    
  call save_buffer (evc, nwordwfc, iunevcn, ik)
    
  tddftpsi(1:npwx,1:nbnd_occ(ik)) = evc(1:npwx,1:nbnd_occ(ik))
    
END SUBROUTINE apply_electric_field_t
  
