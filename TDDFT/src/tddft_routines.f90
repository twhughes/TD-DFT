!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! TODO: eventually, separate input files for FDE, and bcast to images,
! not to comm_world

!-----------------------------------------------------------------------
SUBROUTINE tddft_readin(filin)
  !-----------------------------------------------------------------------
  !
  ! ... Read in the tddft input file. The input file consists of a
  ! ... single namelist &inputtddft. See doc/user-manual.pdf for the
  ! ... list of input keywords.
  !
  USE tddft_module
  USE io_files,         ONLY : prefix, tmp_dir  
  USE io_global,        ONLY : stdout, stdin
  USE constants,        ONLY : bohr_radius_angs, au_sec, pi
  USE input_parameters, ONLY : max_seconds, wf_collect
  USE control_flags,          ONLY : twfcollect
  USE fde
  ! -- parameters --------------------------------------------------------
  implicit none
  character(*), intent(in) :: filin
  ! -- local variables ---------------------------------------------------
  integer :: ios
  character(len=256), external :: trimcheck
  character(len=80) :: verbosity, tdpulse, tefield, integrator, predictor
  namelist /inputtddft/ job, prefix, tmp_dir, conv_threshold, verbosity, &
                        dt, e_strength, e_direction, nstep, nupdate_Dnm, &
                        l_circular_dichroism, l_tddft_restart, max_seconds, &
                        molecule, ehrenfest, tdpulse, t0, tau0, omega0, sigma0, alpha, ttend, &
                        tefield, emaxpos, eopreg, wf_collect, integrator, predictor, update, &
                        pc_conv_threshold, coupled, print_bands

                        

!  if (do_fde) then
     open(unit=stdin, file=filin, form='formatted', status='old')
!  else
!     call input_from_file()
!  endif
!  if (.not. ionode) goto 400
  
  ! define input defult values
  call get_env( 'ESPRESSO_TMPDIR', tmp_dir )
  if (trim(tmp_dir) == ' ') tmp_dir = './scratch/'
  tmp_dir = trimcheck(tmp_dir)
  job          = ''
  prefix       = 'pwscf'
  tmp_dir      = './scratch/'    
  verbosity    = 'high'
  dt           = 2.d0                      ! time step (default: 2 attosecond)
  e_strength   = 0.01d0                    ! impulse electric field strength (default: 0.01/Ang)
  e_direction  = 1                         ! impulse electric field direction: 1-x 2-y 3-z
  conv_threshold = 1.0d-12                 ! convergence threshold    
  nstep        = 1000                      ! total time steps
  nupdate_Dnm  = 1                         ! update USPP Dnm every step
  l_circular_dichroism = .false.
  l_tddft_restart      = .false.
  max_seconds  =  1.d7
  molecule     = .true.
  ehrenfest    = .false.
  tdpulse      = 'cosine'
  t0           = 200.0d0                    ! 0.2 femtoseconds
  tau0         = 200.0d0                    ! 0.2 femtoseconds
  omega0       = 18.0d0                     ! in femtoseconds^-1
  alpha       = 1.0d0                     ! in femtoseconds^-1
  ttend         = 4000.0d0                    ! end time of pulse, 4 femtoseconds
  tefield       = 'wfshift'
  emaxpos       = 0.9d0
  eopreg        = 0.1d0
  eamp          = e_strength
  wf_collect    = .false.
! for the integrators
  integrator = 'cn'
  cn = .false.
  rk = .false.
  em = .false.
  etrs = .false.
  magnus = .false.
  update = .true.
  coupled = .true.
  pmax = 1
  rkmax = 2
  pc_conv_threshold = 1.0d-8
  print_bands = .false.

  
  ! read input    
  read(stdin, inputtddft, err = 200, iostat = ios )

  ! check input
  if (max_seconds < 0.1d0) call errore ('tddft_readin', ' wrong max_seconds', 1)
200 call errore('tddft_readin', 'reading inputtddft namelist', abs(ios))

  select case (trim(verbosity))
     case('low')
       iverbosity = 1
     case('medium')
       iverbosity = 11
     case('high')
       iverbosity = 21
     case default
       call errore('tdddft_readin', 'verbosity can be ''low'', ''medium'' or ''high''', 1)
  end select
   
   select case(trim(tdpulse))
      case('uniform')
          itdpulse=1
          write(stdout,*)'Applying uniform electric field'
      case('cosine')
          itdpulse=2
          write(stdout,*)'Applying damped td electric field'
      case('gaussian')
          itdpulse=3
          write(stdout,*)'Applying damped td electric field'
      case('cosine2')
          itdpulse=4
          ttend=nstep*dt
      case('cosine3')
          itdpulse=5
      case('cosine4')
          itdpulse=6
          write(stdout,*)'Applying undamped td electric field'
      case default
       call errore('tdddft_readin', 'tdpulse can be ''uniform'' or ''cosine''', 1)
  end select
   select case(trim(tefield))
      case('external')
          write(stdout,*)'Applying field through external potential'
          itefield=1
      case('wfshift')
          write(stdout,*)'Applying field through a shift of wavefunction'
          itefield=2
      case default
       call errore('tdddft_readin', 'itefield can be ''external'' or ''wfshift''', 1)
  end select
   select case(trim(integrator))
      case('cnfull')
          write(stdout,*)'Integrating using Self-Consistent Cranck-Nicolson'
          cn = .true.
          pmax = 100
      case('cn4')
          write(stdout,*)'Integrating using Cranck-Nicolson 4th order'
          cn = .true.
          pmax = 4
      case('cn3')
          write(stdout,*)'Integrating using Cranck-Nicolson 3rd order'
          cn = .true.
          pmax = 3
      case('cn2')
          write(stdout,*)'Integrating using Cranck-Nicolson 2nd order'
          cn = .true.
          pmax = 2
      case('cn ')
          write(stdout,*)'Integrating using Cranck-Nicolson first order'
          cn = .true.
          pmax = 1
      case('rk4')
          write(stdout,*)'Integrating using 4th order Runge-Kutta'
          rk = .true.
          rkmax = 4
      case('rk2')
          write(stdout,*)'Integrating using 2nd order Runge-Kutta'
          rk = .true.
          rkmax = 2
      case('em')
          write(stdout,*)'Integrating using Exponentional Middle Point'
          em = .true.
      case('etrs')
          write(stdout,*)'Integrating using Exponentional Time Reversal Symmetry'
          etrs = .true.
      case('magnus')
          write(stdout,*)'Integrating using 4th order Magnus Expansion'
          magnus = .true.
      case default
          write(stdout,*)'Default: Integrating using 2nd order Cranck-Nicolson'
          cn = .true.
          pmax = 2
  end select
  
  twfcollect = wf_collect
  
  if (.not.coupled.and.update) then
     write(stdout,*)'update and uncoupled are not combinable, update set to .false.'
     update = .false.
  endif
          

  ! convert to atomic units
  e_strength = e_strength * bohr_radius_angs ! change from Ryd/Ang to Ryd/Bohr
  dt = dt * 1.d-18 / (2.d0*au_sec)           ! change from femto-second to a.u. (Rydberg unit)
  if (itdpulse==2) then
     t0 = t0 * 1.d-18 / (2.d0*au_sec)           ! change from 0.2 femto-second to a.u. (Rydberg unit)
     tau0 = tau0 * 1.d-18 / (2.d0*au_sec)           ! change from 0.2 femto-second to a.u. (Rydberg unit)
     omega0 = omega0 / 1.d-15 * (2.d0*au_sec)           ! change 18 fs-1 to a.u.
      ttend = ttend * 1.d-18 / (2.d0*au_sec)                   ! end time of pulse, 0.4 femtoseconds
   elseif (itdpulse==3) then
      t0 = t0* 1.d-18 / (2.d0*au_sec) 
   elseif (itdpulse>=4) then
      t0 = t0* 1.d-18 / (2.d0*au_sec) 
      sigma0= sigma0* 1.d-18 / (2.d0*au_sec) 
       omega0 = omega0 * 0.036749309d0 * 2.0d0! omega0 given in eV and changed to au
  endif
 
  eamp = e_strength
!

#ifdef __MPI
  ! broadcast input variables  
  call tddft_bcast_input
#endif

END SUBROUTINE tddft_readin


#ifdef __MPI
!-----------------------------------------------------------------------
SUBROUTINE tddft_bcast_input
  !-----------------------------------------------------------------------
  !
  ! ... Broadcast input data to all processors 
  !
  USE mp_world,         ONLY : world_comm
  USE mp_images,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE io_files,         ONLY : prefix, tmp_dir
  USE input_parameters, ONLY : max_seconds
  USE tddft_module
  USE fde,               ONLY : do_fde

  implicit none
  integer, parameter :: root = 0    

  if (do_fde) then
  call mp_bcast(job, root, intra_image_comm)
  call mp_bcast(prefix, root, intra_image_comm)
  call mp_bcast(tmp_dir, root, intra_image_comm)
  call mp_bcast(dt, root, intra_image_comm)
  call mp_bcast(e_strength, root, intra_image_comm)
  call mp_bcast(e_direction, root, intra_image_comm)
  call mp_bcast(conv_threshold, root, intra_image_comm)
  call mp_bcast(pc_conv_threshold, root, intra_image_comm)
  call mp_bcast(nstep, root, intra_image_comm)
  call mp_bcast(nupdate_Dnm , root, intra_image_comm)
  call mp_bcast(l_circular_dichroism, root, intra_image_comm)
  call mp_bcast(l_tddft_restart, root, intra_image_comm)
  call mp_bcast(iverbosity, root, intra_image_comm)
  call mp_bcast(max_seconds, root, intra_image_comm)
  call mp_bcast(molecule, root, intra_image_comm)
  call mp_bcast(ehrenfest, root, intra_image_comm)
  else
  call mp_bcast(job, root, world_comm)
  call mp_bcast(prefix, root, world_comm)
  call mp_bcast(tmp_dir, root, world_comm)
  call mp_bcast(dt, root, world_comm)
  call mp_bcast(e_strength, root, world_comm)
  call mp_bcast(e_direction, root, world_comm)
  call mp_bcast(conv_threshold, root, world_comm)
  call mp_bcast(pc_conv_threshold, root, world_comm)
  call mp_bcast(nstep, root, world_comm)
  call mp_bcast(nupdate_Dnm , root, world_comm)
  call mp_bcast(l_circular_dichroism, root, world_comm)
  call mp_bcast(l_tddft_restart, root, world_comm)
  call mp_bcast(iverbosity, root, world_comm)
  call mp_bcast(max_seconds, root, world_comm)
  call mp_bcast(molecule, root, world_comm)
  call mp_bcast(ehrenfest, root, world_comm)
  endif

END SUBROUTINE tddft_bcast_input
#endif
  

!-----------------------------------------------------------------------
SUBROUTINE tddft_allocate
  !-----------------------------------------------------------------------
  !
  ! ... Allocate memory for TDDFT
  !
  USE tddft_module
  USE klist,         ONLY : nkstot
  USE wvfct,         ONLY : btype, nbndx
  USE tddft_module

  implicit none

  ! needed by sum_band
  allocate(btype(nbndx,nkstot))
  btype = 1
    
END SUBROUTINE tddft_allocate


!-----------------------------------------------------------------------
SUBROUTINE tddft_summary
  !-----------------------------------------------------------------------
  !
  ! ... Print a short summary of the calculation
  !
  USE io_global,        ONLY : stdout
  USE lsda_mod,         ONLY : nspin
  USE wvfct,            ONLY : nbnd
  USE tddft_module
  implicit none
  integer :: is
    
  write(stdout,*)

  write(stdout,'(5X,''Calculation type      : '',A12)') job
  if (molecule) then
     write(stdout,'(5X,''System is             : molecule'')')
  else
     write(stdout,'(5X,''System is             : crystal'')')
  endif
  if (ehrenfest) write(stdout,'(5X,''Ehrenest dynamics'')')
  write(stdout,'(5X,''Number or steps       : '',I12)') nstep
  write(stdout,'(5X,''Time step             : '',F12.4,'' rydberg_atomic_time'')') dt
  write(stdout,'(5X,''Electric field dir.   : '',I12,'' (1=x,2=y,3=z)'')') e_direction
  write(stdout,'(5X,''Electric field impulse: '',F12.4,'' bohrradius^-1'')') e_strength

  write(stdout,*)

  if (tfixed_occ) then
     write(stdout,'(5X,''Occupations from input:'')')
     do is = 1, nspin
       write(stdout,'(5X,''ispin='',I1,'': '')',advance='no') is
       write(stdout,'(10(F4.2,2X))') f_inp(1:nbnd,is)
     enddo
    write(stdout,*)
  endif
     
  call flush_unit( stdout )

END SUBROUTINE tddft_summary
  
  

!-----------------------------------------------------------------------
SUBROUTINE tddft_openfil
  !
  ! ... Open files needed for TDDFT
  !
  USE tddft_module   
  USE wvfct,            ONLY : nbnd, npwx
  USE ldaU,             ONLY : lda_plus_U, nwfcU
  USE io_files,         ONLY : iunhub, iunwfc, &
                               nwordwfcU, nwordwfc,  seqopn, iunigk
  USE noncollin_module, ONLY : npol
  USE buffers,          ONLY : open_buffer
  USE control_flags,    ONLY : io_level    
  use ruku,             ONLY : iunruku1, iunruku2, iunruku3, iunruku4
  IMPLICIT NONE  
  character*1, parameter :: dir(3) = (/'x', 'y', 'z'/)
  logical :: exst, opnd

  !
  ! ... nwordwfc is the record length (IN REAL WORDS)
  ! ... for the direct-access file containing wavefunctions
  ! ... io_level > 0 : open a file; io_level <= 0 : open a buffer
  !
  nwordwfc = nbnd*npwx*npol
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst )

  ! do not overwrite wfc
  nwordwfc = nbnd*npwx*npol
  CALL open_buffer( iunevcn, 'wfc'//dir(e_direction), nwordwfc, io_level, exst )

  if (rk) then
  ! for ruku
    nwordtdwfc = nbnd*npwx*npol
    CALL open_buffer( iunruku1, 'ruku1', nwordtdwfc, io_level, exst )
    CALL open_buffer( iunruku2, 'ruku2', nwordtdwfc, io_level, exst )
    CALL open_buffer( iunruku3, 'ruku3', nwordtdwfc, io_level, exst )
    CALL open_buffer( iunruku4, 'ruku4', nwordtdwfc, io_level, exst )
  endif

  ! for restart
  nwordtdwfc = nbnd*npwx*npol
  CALL open_buffer( iuntdwfc, 'tmp'//dir(e_direction), nwordtdwfc*2, io_level, exst )

  ! For atomic wavefunctions
  INQUIRE( UNIT = iunigk, OPENED = opnd )
  IF(.NOT. opnd) CALL seqopn( iunigk, 'igk', 'UNFORMATTED', exst )

  ! ... Needed for LDA+U
  ! ... iunhub contains the (orthogonalized) atomic wfcs * S
  nwordwfcU = npwx*nwfcU*npol
  IF ( lda_plus_u ) &
     CALL open_buffer( iunhub, 'hub', nwordwfcU, io_level, exst )

END SUBROUTINE tddft_openfil


!-----------------------------------------------------------------------
SUBROUTINE tddft_closefil
  !-----------------------------------------------------------------------
  !
  ! ... Close files opened by TDDFT
  !
  USE ldaU,             ONLY : lda_plus_U  
  USE io_files,         ONLY : iunhub, iunwfc
  USE buffers,          ONLY : close_buffer

  call close_buffer( iunwfc, 'keep' )
  call close_buffer( iunevcn, 'keep' )
  if ( lda_plus_u ) call close_buffer ( iunhub, status = 'keep' )

END SUBROUTINE tddft_closefil



!-----------------------------------------------------------------------
SUBROUTINE print_clock_tddft
  !-----------------------------------------------------------------------
  !
  ! ... Print clocks
  !
  USE io_global,  ONLY : stdout
  IMPLICIT NONE

  write(stdout,*) '    Initialization:'
  call print_clock ('tddft_setup')
  call print_clock ('gksort')
  write(stdout,*)
  write(stdout,*) '    SCF routines'
  call print_clock ('greenf')
  call print_clock ('ch_psi')
  call print_clock ('h_psi')
  call print_clock ('s_psi')
  call print_clock ('v_of_rho')
  call print_clock ('v_xc')
  call print_clock ('v_h')
  call print_clock ('fde_nonadd')
  call print_clock ('set_vrs')
  call print_clock ('init_us_2')
  call print_clock ('update_rho')
  write(stdout,*)
  write(stdout,*) '    Real time evolution'
  call print_clock ('updateH')
  call print_clock ('crank_nicolson')
  call print_clock ('runge_kutta')
  call print_clock ('dipole')
  call print_clock ('quadrupole')
  call print_clock ('circular')
  write(stdout,*)
  write(stdout,*) '    General routines'
  call print_clock ('calbec')
  call print_clock ('fft')
  call print_clock ('ffts')
  call print_clock ('fftw')
  call print_clock ('cinterpolate')
  call print_clock ('davcio')
  call print_clock ('write_rec')
  write(stdout,*)

#ifdef __MPI
  write(stdout,*) '    Parallel routines'
  call print_clock ('reduce')  
  call print_clock( 'fft_scatter' )
  call print_clock( 'ALLTOALL' )
  write(stdout,*)
#endif
  call print_clock ('TDFFT') 

END SUBROUTINE print_clock_tddft



!-----------------------------------------------------------------------
SUBROUTINE tddft_memory_report
  !-----------------------------------------------------------------------
  !
  ! ... Print estimated memory usage
  !
  USE io_global,                 ONLY : stdout
  USE noncollin_module,          ONLY : npol
  USE uspp,                      ONLY : nkb
  USE fft_base,                  ONLY : dffts
  USE pwcom
  IMPLICIT NONE
  integer, parameter :: Mb=1024*1024, complex_size=16, real_size=8

  ! the conversions to double prevent integer overflow in very large run
  write(stdout,'(5x,"Largest allocated arrays",5x,"est. size (Mb)",5x,"dimensions")')

  write(stdout,'(8x,"KS wavefunctions at k     ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     complex_size*nbnd*npol*DBLE(npwx)/Mb, npwx*npol,nbnd

  write(stdout,'(8x,"First-order wavefunctions ",f10.2," Mb",5x,"(",i8,",",i5,",",i3")")') &
     complex_size*nbnd*npol*DBLE(npwx)*10/Mb, npwx*npol,nbnd,10

  write(stdout,'(8x,"Charge/spin density       ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     real_size*dble(dffts%nnr)*nspin/Mb, dffts%nnr, nspin
  
  write(stdout,'(8x,"NL pseudopotentials       ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     complex_size*nkb*DBLE(npwx)/Mb, npwx, nkb
  write(stdout,*)

END SUBROUTINE tddft_memory_report


!-----------------------------------------------------------------------
SUBROUTINE tddft_read_cards
  !-----------------------------------------------------------------------
  !
  ! ... Read in extra cards
  !
  USE mp_world,         ONLY : world_comm
  USE mp_images,     ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE io_global,        ONLY : ionode, stdout, stdin
  USE fixed_occ,        ONLY : f_inp, tfixed_occ
  USE lsda_mod,         ONLY : nspin
  USE wvfct,            ONLY : nbnd
  USE parser,           ONLY : read_line
  USE tddft_module
  USE fde,               ONLY : do_fde
  implicit none
  character(len=256)         :: input_line
  character(len=80)          :: card
  character(len=1), external :: capital
  logical                    :: tend
  integer                    :: i, is
  integer, parameter         :: root = 0
  
  allocate(f_inp(nbnd,nspin))
  tfixed_occ = .false.
 
  if (.not. ionode) goto 400
 
  read_cards: do
    call read_line(input_line, end_of_file=tend, ionode_only=.true.)
    if (tend) exit read_cards
    if (input_line == ' ' .or. input_line(1:1) == '#' .or. input_line(1:1) == '!' ) &
      cycle read_cards
    read(input_line,*) card
    do i = 1, len_trim(input_line)
       input_line(i: ) = capital(input_line(i:i))
    enddo
  
    if (trim(card) == 'OCCUPATIONS') then
       tfixed_occ = .true.
       do is = 1, nspin
         read(stdin,*) f_inp(1:nbnd,is)
       enddo
   
    else
       write(stdout,*) 'Warning: card '//trim(card)//' ignored'
    endif
    
  enddo read_cards
400 continue
   
  if (do_fde) then
  call mp_bcast(tfixed_occ, root, intra_image_comm)
  call mp_bcast(f_inp, root, intra_image_comm)
  else
  call mp_bcast(tfixed_occ, root, world_comm)
  call mp_bcast(f_inp, root, world_comm)
  endif

  if (tfixed_occ) call weights
!  call weights
  
  return
END SUBROUTINE tddft_read_cards
