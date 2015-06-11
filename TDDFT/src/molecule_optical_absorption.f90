! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine molecule_optical_absorption
  !----------------------------------------------------------------------
  !  ... Compute optical absorption spectrum by real-time TDDFT 
  !  ... References:
  !      (1) Phys. Rev. B 73, 035408 (2006)
  !      (2) http://www.netlib.org/linalg/html_templates/Templates.html
  !                                             Xiaofeng Qian, MIT (2008)
  !----------------------------------------------------------------------
  USE kinds,                       ONLY : dp
  USE io_global,                   ONLY : stdout, ionode
  USE io_files,                    ONLY : nwordwfc, iunwfc, iunigk
  USE cell_base,                   ONLY : at, tpiba, tpiba2, alat
  USE wavefunctions_module,        ONLY : evc
  USE klist,                       ONLY : nks, xk
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, g2kin, current_k, ecutwfc, et
  USE lsda_mod,                    ONLY : current_spin, isk, nspin
  USE becmod,                      ONLY : allocate_bec_type, is_allocated_bec_type, deallocate_bec_type, becp
  USE mp_pools,                    ONLY : inter_pool_comm
  USE mp,                          ONLY : mp_sum, mp_barrier
  USE gvect,                       ONLY : ngm, g
  USE fft_base,                    ONLY : dfftp
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE uspp,                        ONLY : nkb, vkb
  USE scf,                         ONLY : create_scf_type, rho, destroy_scf_type, scf_type_copy
  USE ener,                        ONLY : etot, ehart, eband, ewld, etxc, deband, etxcc
  USE tddft_module
  USE dynamics_module              
  USE fde
  use ruku
  use cnmod
  use constants,                   ONLY : rytoev

  IMPLICIT NONE
  include 'mpif.h'

  !-- tddft variables ----------------------------------------------------
  real(dp), allocatable :: charge(:), dipole(:,:), quadrupole(:,:,:)
  complex(dp), allocatable :: circular(:,:), circular_local(:)
  real(dp), parameter :: Ha2Ryd = 2.d0
!  real(dp), allocatable :: precondition(:)

  integer :: istep_
  integer :: ik, is, ibnd, i
  logical :: file_exists
  integer, external :: find_free_unit
!

!  ! TODO: gk_sort
!     if (nks > 1) rewind (iunigk)
!     do ik = 1, nks
!        current_k = ik
!        current_spin = isk(ik)
!        call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
!        g2kin = g2kin * tpiba2
!        call init_us_2(npw, igk, xk(1,ik), vkb)
!     enddo
!  ! TODO: restart

!!!!
  ! allocate memory
  call allocate_optical()
  if (rk) then
     call allocate_ruku()
  elseif (cn) then
     call allocate_cn()
  endif
  if (.not.coupled) then
     call create_scf_type ( rho_old )
     call scf_type_copy( rho, rho_old ) 
  elseif (.not.update) then
     call create_scf_type(rho_old)
  endif
!  allocate (precondition(npwx))
!  precondition(:) = 1.0d0

  ee = i_complex * dt / 2.d0  ! i*dt/2: do not change
  
  evc = (0.d0,0.d0)
  call tddft_cgsolver_initialize(npwx)
  if (iverbosity > 0 .and. ionode) then
    write(stdout,'(5X,''Done with tddft_cgsolver_initialize'')')
    call flush_unit(stdout)
  endif
 
  ! print the legend
  if (ionode) call print_legend
  
  ! check if we are restarting
  if (l_tddft_restart) then

     if (nks > 1) rewind (iunigk)
     do ik = 1, nks
        current_k = ik
        current_spin = isk(ik)
!*******************************************************************
!!       NEW PART, taken from PW/src
!         npw = ngk(ik)
!         IF (nks > 1) THEN
!            READ (iunigk) igk
!         endif
!         call init_us_2(npw, igk, xk(1,ik), vkb)
!!    
!!    ! read wfcs from file and compute becp
!         evc = (0.d0, 0.d0)
!         call get_buffer(evc, nwordwfc, iunevcn, ik)
        
!*******************************************************************
!        OLD PART
        ! initialize at k-point k 
        CALL start_clock('init_k')
        call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
        g2kin = g2kin * tpiba2
        call init_us_2(npw, igk, xk(1,ik), vkb)
        CALL stop_clock('init_k')
       
!        ! read wfcs from file and compute becp
        evc = (0.d0, 0.d0)
        call get_buffer (evc, nwordwfc, iunevcn, ik)
!*******************************************************************
     end do
     call update_hamiltonian(-1)
     
  end if
 
  if (iverbosity > 0) then
    write(stdout,'(5X,''Done with restart'')')
    call flush_unit(stdout)
  endif

  ! for the time being, kill the MD file
  call seqopn(4, 'md', 'formatted', file_exists)
  close(unit=4, status='delete')


  call calculate_eband(-1)

 ! apply electric field to wavefunction, unless it is done through v_of_rho
 if (itefield==2.and.(.not.l_tddft_restart)) then
   if (itdpulse==1) then
     do ik = 1,nks
        current_k = ik
!       OLD PART
!       ! initialize at k-point k 
        CALL start_clock('init_k')
        call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
        g2kin = g2kin * tpiba2
        call init_us_2(npw, igk, xk(1,ik), vkb)
        CALL stop_clock('init_k')

        if (iverbosity > 0) write(stdout,*)'shifting wavefunction at istep==',istep, 'e_strength=',e_strength
    call flush_unit(stdout)
        call get_buffer (evc, nwordwfc, iunwfc, ik)

!!        tddft_psi(:,:,2) = evc(:,1:nbnd_occ(ik))  ! not shifted
        tddft_psi(:,:,2) = evc(:,1:nbnd)  ! not shifted
        call apply_electric_field(tddft_psi(:,:,1))
        ! save wavefunctions to disk
        call save_buffer (tddft_psi(:,:,1), nwordwfc, iunevcn, ik)
        call save_buffer (tddft_psi, nwordtdwfc*2, iuntdwfc, ik)
     enddo
   endif
 endif

!  if (.not.l_tddft_restart) then
!             ! update the density to be as specified with occupations in input
!             ! file. Done for next steps in update_hamiltonian
!          if (is_allocated_bec_type(becp)) call deallocate_bec_type(becp)
!          call sum_band_tddft(iunwfc)
!          if (do_fde) call update_rho_fde(rho, .true.)
!! add electric field to the potential
!! done for the next steps in update_hamiltonian 
!          if (itefield==1) then
!               DO is = 1, nspin_lsda
!                  CALL add_efield_t(v%of_r(1,is),  rho%of_r)
!               END DO
!               call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    
!          endif
! endif
! in case of field applied through potential, transfer wave function to iunevcn which is used in update_ham.f90
  if ((.not.l_tddft_restart).and.(itefield==1)) then
     do ik = 1,nks
      CALL start_clock('init_k')
      call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      g2kin = g2kin * tpiba2
      call init_us_2(npw, igk, xk(1,ik), vkb)
      CALL stop_clock('init_k')
      call get_buffer(evc,nwordwfc,iunwfc,ik)
      call save_buffer(evc,nwordwfc,iunevcn,ik)
!!      tddft_psi(:,:,1)=evc(:,1:nbnd_occ(ik))
!!      tddft_psi(:,:,2)=evc(:,1:nbnd_occ(ik))
      tddft_psi(:,:,1)=evc(:,1:nbnd)
      tddft_psi(:,:,2)=evc(:,1:nbnd)
      call save_buffer (tddft_psi, nwordtdwfc*2, iuntdwfc, ik)
     enddo
  endif
! PRINT STARTING VALUES
   if (ehrenfest) call molecule_setup_r
   call update_hamiltonian(0)
   if (.not.coupled) call scf_type_copy(rho,rho_old)
   call molecule_compute_dipole( charge, dipole )
#ifdef __PARA
    ! reduce over k-points
    if (l_circular_dichroism) call mp_sum(circular, inter_pool_comm)
    call mp_sum(charge, inter_pool_comm)
    call mp_sum(dipole, inter_pool_comm)
#endif
    if (ionode) then
       if (.not.coupled) write(stdout,*)'UNCOUPLED CALCULATION: ONLY DIPOLE MOMENT PRINTED OUT'
       if (do_fde) then
          if (coupled) write(stdout,'(''ENERGY '',2X,I6,6F16.10)') istep, etot_fde, etot, eband, ehart, etxc, ewld
       else
          write(stdout,'(''ENERGY '',2X,I6,5F16.10)') 0, etot, eband, ehart, etxc, ewld
       endif
       do is = 1, nspin
          if  (iverbosity>10) write(stdout,'(''CHARGE '',I1,1X,I6,3E16.6)') is, 0, charge(is)
          write(stdout,'(''DIP    '',I1,1X,I6,3E16.6)') is, 0, dipole(:,is)
       enddo
    endif
    call flush_unit(stdout)

!*******************************************************************
  ! enter the main TDDFT loop 
  do istep_ = 1, nstep
!*******************************************************************
    istep = istep_
!*******************************************************************
      if (cn) then
        call crank_nicolson
      elseif (rk) then
        call runge_kutta(istep)
      else
         call errore ('molecule_optical_absorption', 'No other integrator schemes implemented yet', 1)
      endif
!*******************************************************************

        
       
      
        
    call update_hamiltonian(istep)
    if (.not.coupled) call scf_type_copy(rho,rho_old)


    call flush_unit(stdout)
  
    if (ehrenfest) then   
!ARK Calculate the forces
     call forces()

!ARK Calculate new positions and velocities with Verlet
     call verlet()
    endif


!    call flush_unit(70+currfrag)
    if (ehrenfest) call trajectory()

!ARK: adding hinit1 because the positions of the atoms were moved
!  terms of the hamiltonian that depend on the nuclear positions are 
!  updated here
    if (ehrenfest)  call hinit1()

!  moved dipole calculation to here, makes more sense at the END of a step...
    if (ehrenfest) call molecule_setup_r
    call molecule_compute_dipole( charge, dipole )
   ! calculate circular dichroism along x, y, and z direction
    if (l_circular_dichroism)  then
      circular_local = (0.d0, 0.d0)
      circular = (0.d0, 0.d0)
      call compute_circular_dichroism(circular_local)
      circular(1:3, current_spin) = circular_local(1:3)
    end if
#ifdef __PARA
    ! reduce over k-points
    if (l_circular_dichroism) call mp_sum(circular, inter_pool_comm)
    call mp_sum(charge, inter_pool_comm)
    call mp_sum(dipole, inter_pool_comm)
#endif
    ! print observables
    if (ionode) then
      do is = 1, nspin
        if  (iverbosity>10) write(stdout,'(''CHARGE '',I1,1X,I6,3E16.6)') is, istep, charge(is)
        write(stdout,'(''DIP    '',I1,1X,I6,3E16.6)') is, istep, dipole(:,is)
      enddo
       if (do_fde) then
          if (coupled) write(stdout,'(''ENERGY '',2X,I6,6F16.10)') istep, etot_fde, etot, eband, ehart, etxc, ewld
!          WRITE( 80, 9060 ) &
!               ( eband + deband ), ehart, ( etxc - etxcc ), ewld
!                write(80,*)
!                write(80,'(5X,''FDE: sum of frag energies ='',F17.8,'' Ry'')') etot_sum
!                write(80,'(5X,''FDE: ewald double counting='',F17.8,'' Ry'')') edc_fde
!                write(80,'(5X,''FDE: non additive kinetic ='',F17.8,'' Ry'')') ekin_nadd
!                write(80,'(5X,''FDE: non additive E_XC    ='',F17.8,'' Ry'')') etxc_nadd
!                write(80,'(5X,''FDE: total energy         ='',F17.8,'' Ry'')') etot_fde
!9060 FORMAT(/'     The total energy is the sum of the following terms:',/,&
!            /'     one-electron contribution =',F17.8,' Ry' &
!            /'     hartree contribution      =',F17.8,' Ry' &
!            /'     xc contribution           =',F17.8,' Ry' &
!            /'     ewald contribution        =',F17.8,' Ry' )
       else
          write(stdout,'(''ENERGY '',2X,I6,5F16.10)') istep, etot, eband, ehart, etxc, ewld
       endif
    endif
    if (print_bands) then
      if (nks>1) then
           DO ik = 1, nks
               WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
               WRITE( stdout, 9030 ) ( et(ibnd,ik) * rytoev, ibnd = 1, nbnd )
           enddo
      else
               WRITE( stdout, 9021 ) 
               WRITE( stdout, 9030 ) ( et(ibnd,1) , ibnd = 1, nbnd )
      endif
    endif
9020 FORMAT('          k =',3F7.4,'     band energies (ev):')
9021 FORMAT('          band energies (ev):')
9030 FORMAT( '  ',8F9.4 )
     
    call flush_unit(stdout)
     
!*******************************************************************
  enddo      ! end of TDDFT loop
!*******************************************************************

  ! finish  
  call tddft_cgsolver_finalize()
  call deallocate_optical()
  call deallocate_dyn_vars()
  if (rk) then 
        call deallocate_ruku() 
        call ruku_closefil()
  elseif (cn) then
        call deallocate_cn()
  endif
  if ((.not.coupled).or.(.not.update)) call destroy_scf_type(rho_old)
   
    
CONTAINS

  !====================================================================
  ! Print the legend key
  !====================================================================    
  SUBROUTINE print_legend
    write(stdout,'(5X,''Output quantities:'')')
    write(stdout,'(5X,''  CHARGE spin  istep  charge'')')
    write(stdout,'(5X,''  DIP    spin  istep  dipole(1:3)'')')
    if (do_fde) then
       if (coupled) write(stdout,'(5X,''  ENERGY istep etot_fde etot eband ehart etxc ewld'')')
    else
       write(stdout,'(5X,''  ENERGY istep etot eband ehart etxc ewld'')')
    endif
    write(stdout,*)
    call flush_unit(stdout)
  END SUBROUTINE print_legend

  
  !====================================================================
  ! Initialize and allocate memory
  !====================================================================    
  SUBROUTINE allocate_optical()
    USE becmod, ONLY : becp, allocate_bec_type
    IMPLICIT NONE
    integer :: ik
    
    nbnd_occ_max = 0
    do ik = 1, nks
      if (nbnd_occ(ik) > nbnd_occ_max) nbnd_occ_max = nbnd_occ(ik)
    enddo

    call allocate_bec_type(nkb, nbnd, becp)
   
    allocate (tddft_psi (npwx,nbnd,2))
!!    allocate (tddft_hpsi(npwx,nbnd_occ_max))
!!    allocate (tddft_spsi(npwx,nbnd_occ_max))
!!    if (ehrenfest) allocate (tddft_Ppsi(npwx,nbnd_occ_max))
!!    allocate (b(npwx,nbnd_occ_max))
    allocate (tddft_hpsi(npwx,nbnd))
    allocate (tddft_spsi(npwx,nbnd))
    if (ehrenfest) allocate (tddft_Ppsi(npwx,nbnd))
    allocate (b(npwx,nbnd))
!    tddft_psi = (0.d0,0.
    tddft_psi = 0.0d0
    tddft_hpsi = (0.d0,0.d0)
    tddft_spsi = (0.d0,0.d0)
    b = (0.d0,0.d0)

    allocate (charge(nspin), dipole(3,nspin), quadrupole(3,3,nspin))
    allocate (circular(3,nspin), circular_local(3))
    charge = 0.d0
    dipole = 0.d0
    quadrupole = 0.d0
    circular = (0.d0, 0.d0)
    circular_local = (0.d0, 0.d0)

    allocate (r_pos(3,dfftp%nnr), r_pos_s(3,dfftp%nnr))
    call molecule_setup_r
    
  END SUBROUTINE allocate_optical
  
  
  !====================================================================
  ! Deallocate memory
  !====================================================================    
  SUBROUTINE deallocate_optical()
    USE becmod, ONLY : becp, deallocate_bec_type
    IMPLICIT NONE

    call deallocate_bec_type(becp)
    deallocate (tddft_psi, tddft_hpsi, tddft_spsi, b)
    if (ehrenfest) deallocate(tddft_Ppsi)
    deallocate (charge, dipole, quadrupole, circular, circular_local)
    deallocate (r_pos, r_pos_s)
    
  END SUBROUTINE deallocate_optical
   
   
  !====================================================================
  ! compute circular dichroism (EXPERIMENTAL, NORM-CONSERVING ONLY)
  !====================================================================      
  subroutine compute_circular_dichroism(cd)
    USE fft_base,               ONLY : dfftp
    USE fft_interfaces,         ONLY : invfft
    USE mp_global,              ONLY : me_pool
    USE gvecs,                  ONLY : nls
    IMPLICIT NONE
    REAL(DP) :: xx(dfftp%nnr), yy(dfftp%nnr), zz(dfftp%nnr), gk
    INTEGER  :: ik, ibnd, i, ii, jj, kk, index0, index, ir, ipol, ind, i_current_spin, ig
    complex(dp) :: p_psi(npwx), p_psi_r(dfftp%nnr, 3), cd(3, nspin), psic1(dfftp%nnr)
    
    xx(:) = 0.d0
    yy(:) = 0.d0
    zz(:) = 0.d0
    
    index0 = 0

#ifdef __PARA
  do i = 1, me_pool
    index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  enddo
#endif

  ! loop over real space grid
  do ir = 1, dfftp%nnr
    index = index0 + ir - 1
    kk     = index / (dfftp%nr1x*dfftp%nr2x)
    index = index - (dfftp%nr1x*dfftp%nr2x)*kk
    jj     = index / dfftp%nr1x
    index = index - dfftp%nr1x*jj
    ii     = index

             xx(ir) = &
                  dble( ii-1 )/dble(dfftp%nr1) * at(1,1) * alat + &
                  dble( jj-1 )/dble(dfftp%nr2) * at(1,2) * alat + &
                  dble( kk-1 )/dble(dfftp%nr3) * at(1,3) * alat
             
             yy(ir) = &
                  dble( ii-1 )/dble(dfftp%nr1) * at(2,1) * alat + &
                  dble( jj-1 )/dble(dfftp%nr2) * at(2,2) * alat + &
                  dble( kk-1 )/dble(dfftp%nr3) * at(2,3) * alat
             
             zz(ir) = &
                  dble( ii-1 )/dble(dfftp%nr1) * at(3,1) * alat + &
                  dble( jj-1 )/dble(dfftp%nr2) * at(3,2) * alat + &
                  dble( kk-1 )/dble(dfftp%nr3) * at(3,3) * alat
             
    end do
    
    cd(:,:) = (0.d0, 0.d0)
    
    do ik = 1, nks
       
!!       if (nbnd_occ(ik) > 0) then
       if (nbnd > 0) then
          
! AK: Not sure if calling gk_sort here is correct...
       i_current_spin = isk(ik)
       CALL start_clock('init_k')
       call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
       g2kin(:) = g2kin(:) * tpiba2
       CALL stop_clock('init_k')
       
!!       do ibnd = 1, nbnd_occ(ik)
       do ibnd = 1, nbnd
          
          p_psi_r(:, :) = (0.d0, 0.d0)
          do ipol = 1, 3
             p_psi(:) = (0.d0, 0.d0)
             do ig = 1, npw
                gk = xk(ipol,ik) + g(ipol,igk(ig))
                p_psi(ig) = gk * tpiba * tddft_psi(ig, ibnd, ik)
             end do
             psic1(:) = (0.d0, 0.d0)
             psic1(nls(igk(1:npw))) = p_psi(:)
             call invfft('Wave', psic1, dfftp)
             p_psi_r(:,ipol) = psic1(:)
          end do
          
          ! transform wavefunction from reciprocal space into real space
          psic1(:) = (0.d0, 0.d0)
          psic1(nls(igk(1:npw))) = tddft_psi(1:npw, ibnd, ik)
          call invfft('Wave', psic1, dfftp)
          
          do ind = 1, dfftp%nnr
             cd(1, i_current_spin) = cd(1, i_current_spin) + &
                  conjg(psic1(ind)) * ( yy(ind) * p_psi_r(ind,3) - zz(ind) * p_psi_r(ind,2) )
             cd(2, i_current_spin) = cd(2, i_current_spin) + &
                  conjg(psic1(ind)) * ( zz(ind) * p_psi_r(ind,1) - xx(ind) * p_psi_r(ind,3) )
             cd(3, i_current_spin) = cd(3, i_current_spin) + &
                  conjg(psic1(ind)) * ( xx(ind) * p_psi_r(ind,2) - yy(ind) * p_psi_r(ind,1) )
          end do
          
          
       end do

       end if
       
    end do
    cd = cd  / dble(dfftp%nnr)
    
    
    RETURN
  end subroutine compute_circular_dichroism

END SUBROUTINE molecule_optical_absorption
 

