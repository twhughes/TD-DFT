!
! Copyright (C) 2001-2015 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
!#ifdef __BANDS
!subroutine tddft_ch_psi_all (n, h, ah, ee, m, ibnd_start, ibnd_end)
!#else
subroutine tddft_ch_psi_all (n, h, ah, ee, m, ip, vrs_deltat, v_deltat)
!#endif
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( S + ee * H), where, ee = i * dt/2
  ! to a vector h. The result is given in Ah.
  !
  USE kinds,        ONLY : dp
  USE wvfct,        ONLY : npwx
  USE becmod,       ONLY : calbec
  USE mp,           ONLY : mp_sum
  use scf,          only : scf_type
  use io_global,    only : stdout
  USE fft_base,   ONLY : dfftp
  USE lsda_mod,   ONLY : nspin
#ifdef __BANDS
  USE mp_bands,     ONLY : intra_bgrp_comm
#endif
  implicit none

  integer :: n, m
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point
#ifdef __BANDS
  integer, intent(in) :: ibnd_start, ibnd_end
#endif
  integer, intent(in) :: ip
  type (scf_type), intent(inout) :: v_deltat
  real(dp), intent(inout) :: vrs_deltat(dfftp%nnr,nspin)
  complex(DP) :: ee
  ! input: i*dt/2

  complex(DP) :: h (npwx, m), ah (npwx, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ig

  complex(dp), allocatable :: hpsi (:,:), hpsi_deltat(:,:)
  
  call start_clock ('ch_psi')
  allocate (hpsi( npwx , m))
  if (ip>1) allocate(hpsi_deltat(npwx,m))
  hpsi (:,:) = (0.d0, 0.d0)
  if (ip>1) hpsi_deltat (:,:) = (0.d0, 0.d0)

  !   compute the product of the hamiltonian H and ultrasoft projection operator S with the h vector
  !   TODO: band version
  if (ip==1) then
     call h_psi (npwx, n, m, h, hpsi)
  else
     call h_psi (npwx, n, m, h, hpsi)
     call h_psi_deltat(npwx, n, m, vrs_deltat, v_deltat, h, hpsi_deltat)
     hpsi=(hpsi+hpsi_deltat)*0.5d0
  endif
  call s_psi (npwx, n, m, h, ah)

  call start_clock ('last')
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) +  ee * hpsi (ig, ibnd)
     enddo
  enddo

  deallocate (hpsi)
  if (ip>1) deallocate(hpsi_deltat)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return

end subroutine tddft_ch_psi_all
