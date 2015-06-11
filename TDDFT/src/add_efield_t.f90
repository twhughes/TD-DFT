!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... written by J. Tobik
!
! Changes 30/06/2003 (ADC) : 
!               Calculation of corrections to energy and forces due
!               to the field.
!               Added possibility to subtract the dipole field 
!               for slab or molecule calculation.
!               (See Bengtsson PRB 59, 12 301 (1999) and
!                    Meyer and Vanderbilt, PRB 63, 205426 (2001).)
!
!          25/06/2009 (Riccardo Sabatini)
!               reformulation using a unique saw(x) function (included in 
!               cell_base) in all e-field related routines and inclusion of 
!               a macroscopic electronic dipole contribution in the mixing 
!               scheme. 
!

!
!--------------------------------------------------------------------------
SUBROUTINE add_efield_t(vpoten,write_field)
  !--------------------------------------------------------------------------
  !
  !   This routine adds an electric field to the local potential. The
  !   field is made artificially periodic by introducing a saw-tooth
  !   potential. The field is parallel to a reciprocal lattice vector bg, 
  !   according to the index e_direction.
  !
  !   if dipfield is false the electric field correction is added to the
  !   potential given as input (the bare local potential) only
  !   at the first call to this routine. In the following calls
  !   the routine exit.
  !
  !   if dipfield is true the dipole moment per unit surface is calculated
  !   and used to cancel the electric field due to periodic boundary
  !   conditions. This potential is added to the Hartree and xc potential
  !   in v_of_rho. NB: in this case the electric field contribution to the 
  !   band energy is subtracted by deband.
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  use tddft_module
  USE io_global,     ONLY : stdout
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dfftp
  USE mp,            ONLY : mp_bcast, mp_sum
  use constants,      only : au_sec, pi
  
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr)! ef is added to this potential
! REAL(DP),INTENT(INOUT) :: etotefield       ! contribution to etot due to ef
  !
  ! local variables
  !
  INTEGER :: index0, i, j, k
  INTEGER :: ir
  REAL(DP) :: length, vamp, value, sawarg, e_dipole, ion_dipole, t
  REAL(DP) :: tot_dipole, bmod

  LOGICAL :: first=.TRUE.
  logical :: write_field
  SAVE first
  
  !---------------------
  !  Execution control
  !---------------------

  IF (itefield>1) RETURN ! itefield==1 is external field
  ! efield only needs to be added on the first iteration, if dipfield
  ! is not used. note that for relax calculations it has to be added
  ! again on subsequent relax steps.

   IF ((e_direction.lt.1).or.(e_direction.gt.3)) THEN
      CALL errore('add_efield',' wrong e_direction',1)
   ENDIF

  !---------------------
  !  Variable initialization
  !---------------------

  bmod=SQRT(bg(1,e_direction)**2+bg(2,e_direction)**2+bg(3,e_direction)**2)

  tot_dipole=0._dp
  e_dipole  =0._dp
  ion_dipole=0._dp
  
  !---------------------
  !  Calculate dipole
  !---------------------
  
!  if (dipfield) then
!  !
!  ! dipole correction is active 
!  !
!     CALL compute_el_dip(emaxpos, eopreg, e_direction, rho, e_dipole)
!     CALL compute_ion_dip(emaxpos, eopreg, e_direction, ion_dipole)
!    
!     tot_dipole  = -e_dipole + ion_dipole
!     CALL mp_bcast(tot_dipole, 0, intra_image_comm)
!  !  
!  !  E_{TOT} = -e^{2} \left( eamp - dip \right) dip \frac{\Omega}{4\pi} 
!  !
!     etotefield=-e2*(eamp-tot_dipole/2.d0)*tot_dipole*omega/fpi 
!
!  !---------------------
!  !  Define forcefield
!  !  
!  !  F_{s} = e^{2} \left( eamp - dip \right) z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
!  !---------------------
!    
!     IF (lforce) THEN
!        DO na=1,nat
!           DO ipol=1,3
!              forcefield(ipol,na)= e2 *(eamp - tot_dipole) &
!                               *zv(ityp(na))*bg(ipol,e_direction)/bmod
!           ENDDO
!        ENDDO
!     ENDIF
!
!  else
!  !
!  ! dipole correction is not active
  !

     CALL compute_ion_dip(emaxpos, eopreg, e_direction, ion_dipole)

  !  
  !  E_{TOT} = -e^{2} eamp * iondip \frac{\Omega}{4\pi} 
  !
     t = dt*istep
     if (itdpulse==2) eamp = real(e_strength * cos(pi*(t-2.0d0*tau0-t0)/(2.0d0*tau0))*exp(cmplx(0,omega0*t,kind=dp)),kind=dp)
     if (itdpulse==3) eamp = e_strength * exp(-omega0*(t-t0)**2)
     etotefield=-e2*eamp*ion_dipole*omega/fpi 
     if (itdpulse==4) eamp = e_strength * cos(omega0*t)
     if (itdpulse==5) eamp = e_strength * cos(omega0*t) * exp(-0.05*t)
     if (itdpulse==6) eamp = e_strength * cos(omega0*t) * exp(-(t-t0)**2/(2*sigma0**2)) ! make half width be a bit less than place of maximum
     if (istep < nstep.and.write_field) write(stdout,*)'efield info',istep,eamp,etotefield
!      write(stdout,*)'e_strength',e_strength
!      write(stdout,*)'t0,tau0',t0,tau0
!      write(stdout,*)'omega0',omega0

  !---------------------
  !  Define forcefield
  !  
  !  F_{s} = e^{2}  eamp z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
  !---------------------
    
!!     IF (lforce) THEN
!        DO na=1,nat
!           DO ipol=1,3
!              forcefield(ipol,na)= e2 *eamp &
!                               *zv(ityp(na))*bg(ipol,e_direction)/bmod
!           ENDDO
!        ENDDO
!     ENDIF

!  end if

  !
  !  Calculate potential and print values 
  !   
  
  length=(1._dp-eopreg)*(alat*SQRT(at(1,e_direction)**2+at(2,e_direction)**2+at(3,e_direction)**2))
  
  vamp=e2*(eamp-tot_dipole)*length

! IF (ionode) THEN
       !
       ! Output data
       !
!      WRITE( stdout,*)
!      WRITE( stdout,'(5x,"Adding external electric field":)')

!       IF (dipfield) then
!          WRITE( stdout,'(/5x,"Computed dipole along e_direction(",i1,") : ")' ) e_direction
!
!          !
!          !  If verbose prints also the different components
!          !
!          IF ( iverbosity > 0 ) THEN
!              WRITE( stdout, '(8X,"Elec. dipole ",1F15.4," Ry au,  ", 1F15.4," Debye")' ) &
!                                            e_dipole, (e_dipole*au_debye)
!              WRITE( stdout, '(8X,"Ion. dipole  ",1F15.4," Ry au,", 1F15.4," Debye")' ) &
!                                          ion_dipole, (ion_dipole*au_debye)
!          ENDIF
!
!          WRITE( stdout, '(8X,"Dipole       ",1F15.4," Ry au, ", 1F15.4," Debye")' ) &
!                                            (tot_dipole* (omega/fpi)),   &
!                                            ((tot_dipole* (omega/fpi))*au_debye)  
!
!          WRITE( stdout, '(8x,"Dipole field     ", f11.4," Ry au")') tot_dipole
!          WRITE( stdout,*)
!
!       ENDIF

!       IF (abs(eamp)>0._dp) WRITE( stdout, &
!          '(8x,"E field amplitude [Ha a.u.]: ", es11.4)') eamp 
!        
!       WRITE( stdout,'(8x,"Potential amp.   ", f11.4," Ry")') vamp 
!       WRITE( stdout,'(8x,"Total length     ", f11.4," bohr")') length
!       WRITE( stdout,*)     
!  ENDIF
!
!
  !
  !------------------------------
  !  Add potential
  !  
  !  V\left(ijk\right) = e^{2} \left( eamp - dip \right) z_{v} 
  !          Saw\left( \frac{k}{nr3} \right) \frac{alat}{bmod} 
  !          
  !---------------------

  ! Index for parallel summation
  !
  index0 = 0
#if defined (__MPI)
  !
  DO i = 1, me_bgrp
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
  !
#endif
  !
  ! Loop in the charge array
  !

  DO ir = 1, dfftp%nnr
     !
     ! ... three dimensional indexes
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
     
     if (e_direction.eq.1) sawarg = DBLE(i)/DBLE(dfftp%nr1)
     if (e_direction.eq.2) sawarg = DBLE(j)/DBLE(dfftp%nr2)
     if (e_direction.eq.3) sawarg = DBLE(k)/DBLE(dfftp%nr3)
     
     value = e2*(eamp - tot_dipole)*saw(emaxpos,eopreg,sawarg) * (alat/bmod)

     vpoten(ir) = vpoten(ir) + value

  END DO
  
  
  RETURN

END SUBROUTINE add_efield_t
