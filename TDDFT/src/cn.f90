!-----------------------------------------------------------------------
MODULE cnmod
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for the Runge-Kutta integration scheme
  !
  USE kinds,            ONLY : dp
  USE wvfct,            ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k, ecutwfc, et
  USE klist,            ONLY : nks, nkstot, wk, xk, nelec, ngk
  USE buffers,          ONLY : close_buffer
  use tddft_module,     only : nbnd_occ_max
  use scf,              only : scf_type, create_scf_type, destroy_scf_type
  USE fft_base,         ONLY : dfftp
  use lsda_mod,         only : nspin
  IMPLICIT NONE
  SAVE

  complex(dp), allocatable :: tddft_hpsi_deltat(:,:)
  type (scf_type) :: v_deltat,  rho_deltat
  real(dp), allocatable :: vrs_deltat(:,:)
  real(dp) :: norm_rho
  real(dp), external :: dnrm2
 
  
Contains
  
  subroutine allocate_cn
!   allocate (tddft_hpsi_deltat(npwx,nbnd_occ_max))
   allocate (tddft_hpsi_deltat(npwx,nbnd))
   call create_scf_type(v_deltat, do_not_allocate_becsum = .true.)
!   call create_scf_type(v_deltat, do_not_allocate_becsum = .false.)
   call create_scf_type(rho_deltat)
   allocate(vrs_deltat(dfftp%nnr, nspin))
  end subroutine allocate_cn

  subroutine deallocate_cn
!   if (pmax > 1) then
     deallocate (vrs_deltat)
     call destroy_scf_type(rho_deltat)
     call destroy_scf_type(v_deltat)
     deallocate (tddft_hpsi_deltat)
!   endif 
  end subroutine deallocate_cn
 
END MODULE cnmod

