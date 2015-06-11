!
! Copyright (C) 2001-2010 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! TODO: nsave, restart_mode

!-----------------------------------------------------------------------
MODULE tddft_module
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for TDDFT calculations
  !
  USE kinds,               ONLY : dp
  USE dynamics_module,     ONLY : dt           ! timestep
  USE control_flags,       ONLY : istep, nstep ! current step and number of steps
  USE fixed_occ,           ONLY : tfixed_occ, f_inp ! occupations from input
  IMPLICIT NONE
  SAVE
  
  character(80) :: job             ! 'optical'
  integer  :: e_direction          ! impulse electric field direction: 1-x 2-y 3-z
  real(dp) :: e_strength           ! impulse electron field strength
  real(dp) :: conv_threshold       ! cg convergence threshold
  real(dp) :: pc_conv_threshold       ! predictor-corrector convergence threshold
  integer  :: nupdate_Dnm          ! update USPP Dnm matrix every n steps
  logical  :: l_circular_dichroism ! calculate circular dichroism
  logical  :: l_tddft_restart      ! restart propagation from the last step
  integer  :: iverbosity           ! verbosity level (default = 1)
  logical  :: molecule             ! use molecular routuines
  logical  :: ehrenfest            ! .true. if moving the ions
  integer  :: itdpulse              ! .true. if time dependent electric field
  real(dp) :: t0, tau0, omega0, ttend, sigma0 ! parameters for tdpulse
  real(dp) :: alpha
  integer  :: itefield
  logical  :: cn, rk, em, etrs, magnus, update, coupled, print_bands
  integer  :: pmax, rkmax
  REAL(DP) :: &
  ! both emaxpos and eopreg need to be in the vacuum, close to each other 
      emaxpos,       &! position of the maximum of the field (0<emaxpos<1)
      eopreg,        &! amplitude of the inverse region (0<eopreg<1)
      eamp,          &
      etotefield      
  
  complex(dp), parameter :: i_complex = (0.0_dp,1.0_dp)

  real(dp), allocatable :: r_pos(:,:)     ! position operator in real space
  real(dp), allocatable :: r_pos_s(:,:)   ! position operator in real space (smooth grid)
  integer, allocatable :: nbnd_occ(:)     ! occupied bands for each k-point
  integer :: nbnd_occ_max                 ! max number of occupied bands

  integer, parameter :: iuntdwfc = 51     ! to save TDDFT intermediate wfcs
  integer :: nwordtdwfc 
  integer, parameter :: iunevcn = 52      ! evc for restart
  real(dp) :: alpha_pv                    ! shift of conduction levels

  integer :: tddft_exit_code = 0

  complex(dp) :: ee                     ! i*dt/2

  complex(dp), allocatable :: tddft_psi(:,:,:), b(:,:)
  complex(dp), allocatable :: tddft_hpsi(:,:), tddft_spsi(:,:), tddft_Ppsi(:,:)

END MODULE tddft_module

