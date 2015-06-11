!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
!SUBROUTINE sinv_hpsi( lda, n, m, hpsi, spsi )
SUBROUTINE sinv_hpsi( lda, m, hpsi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the S matrix to m wavefunctions psi
  ! ... and puts the results in spsi.
  ! ... Requires the products of psi with all beta functions
  ! ... in array becp(nkb,m) (calculated in h_psi or by calbec)
  !
  ! ... input:
  !
  ! ...    lda   leading dimension of arrays psi, spsi
  ! ...    n     true dimension of psi, spsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  !
  ! ...    spsi  S*psi
  !
  USE kinds,      ONLY : DP
!  USE uspp,       ONLY : vkb, nkb, qq, okvan
  USE uspp,       ONLY : nkb, okvan
!  USE uspp_param, ONLY : upf, nh 
!  USE ions_base,  ONLY : nat, nsp, ityp
!  USE control_flags,    ONLY: gamma_only 
  USE noncollin_module, ONLY: npol
!  USE noncollin_module, ONLY: npol, noncolin
!  USE realus,     ONLY :  real_space, fft_orbital_gamma, initialisation_level,&
!                          bfft_orbital_gamma, calbec_rs_gamma, s_psir_gamma
  !
  IMPLICIT NONE
  !
!  INTEGER, INTENT(IN) :: lda, n, m
  INTEGER, INTENT(IN) :: lda, m
  COMPLEX(DP), INTENT(IN) :: hpsi(lda*npol,m)
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !
!  INTEGER :: ibnd
  !
  ! ... initialize  spsi
  !
  spsi = hpsi
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) then
      RETURN
  else
      call errore ('sinv_hpsi', 'Not implemented yet for ultrasoft', 1)
  endif
  !
END SUBROUTINE sinv_hpsi
