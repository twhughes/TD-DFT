!
! Copyright (C) 2001-2014 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Fortran implementation of the Conjugate Gradient Square solver
! author: Xiaofeng Qian, MIT (2008)
!

!-----------------------------------------------------------------------
MODULE tddft_cgsolver_module
  !-----------------------------------------------------------------------
  !
  ! ... this module contains variables and temporaries for the CGS solver
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  SAVE

  integer     :: flag
  real(dp)    :: tolb, normr, relres, n2b, normrmin, realnormr
  complex(dp) :: rho, rho1, alpha, beta, rtvh, aconst
  complex(dp), allocatable :: r(:),  rt(:), vh(:), u(:), &
                              q(:), qh(:), p(:), ax(:)

END MODULE tddft_cgsolver_module


!-----------------------------------------------------------------------
SUBROUTINE tddft_cgsolver_initialize(ndmx)
  !-----------------------------------------------------------------------
  !
  ! ... allocate memory for the solver
  !
  USE tddft_cgsolver_module
  implicit none
  integer, intent(in) :: ndmx
  integer :: ierr
    
  allocate (r(ndmx), rt(ndmx),  u(ndmx), p(ndmx), q(ndmx), &
            qh(ndmx),  vh(ndmx), ax(ndmx), stat=ierr )
  if (ierr/=0) call errore('tddft_cgsolver_initialize','error allocating',ierr)
    
END SUBROUTINE tddft_cgsolver_initialize
  
  
!-----------------------------------------------------------------------
SUBROUTINE tddft_cgsolver_finalize()
  !-----------------------------------------------------------------------
  !
  ! ... deallocate memory
  !
  USE tddft_cgsolver_module
  implicit none
  deallocate( r, rt,  u, p, q, qh,  vh , ax)
    
END SUBROUTINE tddft_cgsolver_finalize
  

!----------------------------------------------------------------------
SUBROUTINE tddft_cgsolver (A, b, x, ndmx, ndim, tol, iter, flag_global,  &
                           nbnd, ee, ip, vrs_deltat, v_deltat)
  !----------------------------------------------------------------------
  !
  ! ... Conjugate-Gradient Square method for solving:   A * x = b
  ! ... where: A*x is evaluated by subroutine 'A', and 'A' is implicit
  ! ... general square-matrix.
  !                                            Xiaofeng Qian, MIT (2008)
  USE kinds,     ONLY : dp
  USE mp_pools,  ONLY : intra_pool_comm
  USE mp,        ONLY : mp_sum
  USE tddft_cgsolver_module
  use io_global, only : stdout
  use scf,       only : scf_type
  USE fft_base,   ONLY : dfftp
  USE lsda_mod,   ONLY : nspin
  !----------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  integer, intent(in) :: ndmx               ! the maximum dimension of the vectors
  integer, intent(in) :: ndim               ! the actual dimension of the vectors
  integer, intent(in) :: nbnd               ! the number of bands
  integer, intent(out) :: iter              ! counter on iterations
  integer, intent(out) :: flag_global
  integer, intent(in)  :: ip
  real(dp), intent(in) :: tol               ! threshold for convergence
  real(dp), intent(in) :: vrs_deltat(dfftp%nnr,nspin)
  complex(dp), intent(in) :: ee             ! i*dt/2
  complex(dp), intent(in) :: b(ndmx,nbnd)   ! input: the known term
  complex(dp), intent(inout) :: x(ndmx,nbnd)  ! the solution of the linear system
  type (scf_type), intent(inout) :: v_deltat
  external A                                ! the subroutine computing A*x
  !----------------------------------------------------------------------
  integer, parameter :: maxit = 200          ! the maximum number of iterations
  complex(dp), external :: zdotc
  real(dp), external    :: ddot
  integer :: imin, stag, i, ibnd
! for the timing
  double precision :: step_start, step_end
  

  if (.not. allocated(r)) call errore('tddft_cgsolver', 'cgsolver not initialized', 1)
  
  call start_clock ('cgsolver')

  ! initialize module variables
  tolb = 0.d0
  n2b  = 0.d0
  relres = 0.d0
  normr= 0.d0
  rho  = 0.d0
  rho1 = 0.d0
  r    = (0.d0, 0.d0)
  aconst = (1.0d0, 0.0d0)
  iter = 0
  
  flag           = 1
  imin           = 0
  
  !----------------------------------------------------------------------
  ! loop over bands
  !----------------------------------------------------------------------
  do ibnd = 1, nbnd

     n2b = dble(zdotc(ndim, b(1,ibnd), 1, b(1,ibnd), 1))
#ifdef __MPI
     call mp_sum(n2b, intra_pool_comm)
#endif
     n2b = dsqrt(n2b)
     tolb = tol * n2b

     call A(ndim, x(1,ibnd), r, ee, 1, ip, vrs_deltat, v_deltat)
!      r = (1.0d0,1.0d0)
     
!     r = b(:,ibnd) - Ax
     call zaxpy(ndim, -1.0d0, b(1,ibnd), 1, r, 1)
     call zscal(ndim, -1.0d0, r, 1)
     normr = dble(zdotc( ndim, r, 1, r, 1))
#ifdef __MPI
     call mp_sum(normr, intra_pool_comm)
#endif  
     normr = dsqrt(normr)
     
     if (normr < tolb) then
        flag = 0
        cycle
     endif
     
     rt = r
     normrmin = normr
     stag = 0
     rho  = cmplx(1.d0, 0.d0, kind=dp)
     
     ! CG iteration
     do i = 1, maxit
        step_start = mpi_wtime() 
        rho1 = rho
        rho = zdotc(ndim, rt, 1, r, 1)
#ifdef __MPI
        call mp_sum(rho, intra_pool_comm)
#endif
        
        if (rho == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag  = 2
           flag_global = flag
           return
        endif

        if (i == 1) then
!           u = r
           call zcopy(ndim,r,1,u,1)
!           p = u
           call zcopy(ndim,u,1,p,1)
        else
           beta = rho / rho1
           if ( beta == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
              flag  = 3
              flag_global = flag
              return
           endif
!           u = r + beta * q
! translating to blas in 2 steps
!          u = r
           call zcopy(ndim,r,1,u,1)
!          u = beta*q + u
           call zaxpy(ndim,beta,q,1,u,1)
!           p = u + beta * (q + beta * p)
! translating to blas in three steps
!           q = beta*p + q
           call zaxpy(ndim, beta, p, 1, q, 1)
!           p = u
           call zcopy(ndim,u,1,p,1)
!           p = beta*q + p
           call zaxpy(ndim, beta, q, 1, p, 1)

        endif
        
        call A(ndim, p, vh, ee, 1, ip, vrs_deltat, v_deltat)
        
        rtvh = zdotc(ndim, rt, 1, vh, 1)
#ifdef __MPI
        call mp_sum(rtvh, intra_pool_comm)
#endif
        
        if (rtvh == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag = 4
           flag_global = flag
           return
        else
           alpha = rho / rtvh
        endif
        
        if (alpha == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag = 5
           stag = 1
           flag_global = flag
           return
        endif
        
!        q  = u - alpha * vh
! translate to blas in two steps:
!       q = u
           call zcopy(ndim,u,1,q,1)
!       q = -alpha*vh + q
           call zaxpy(ndim, -1.d0*alpha, vh, 1, q, 1)
        
!        uh = u + q
!         uh = u
!           call zcopy(ndim, u, 1, uh, 1)
!         uh = aconst*q + uh
           call zaxpy(ndim,aconst,  q, 1, u, 1)
        

        
!        x(:,ibnd)  = x(:,ibnd) + alpha * uh
! chagned to:        x(:,ibnd)  = x(:,ibnd) + alpha * u
           call zaxpy(ndim, alpha, u, 1, x(1,ibnd), 1)

! instead of actually calculating B-Ax each iteration with call A, how about
! looking at the norm of the residue r? If r is small, we can still do the check
! just once.
        
        call A(ndim, u, qh, ee, 1, ip, vrs_deltat, v_deltat)
        
!        r  = r - alpha * qh
           call zaxpy(ndim, -alpha, qh, 1, r, 1)
! new part for estimation of convergence
         normr = dble(zdotc( ndim, r, 1, r, 1))
#ifdef __MPI
         call mp_sum(normr, intra_pool_comm)
#endif  
    
        if (normr <= tolb) then
           flag = 0
           flag_global = flag
           iter = i
           exit
        endif
        
        if (stag == 1) then
           flag  = 5
           flag_global = flag
           return
        endif
        
        if (normr < normrmin)  then
           normrmin = normr
           imin = i
        endif
     
        step_end = mpi_wtime() 
        
     enddo ! i

  end do ! ibnd
  !----------------------------------------------------------------------
  ! end of the loop over bands
  !----------------------------------------------------------------------
  
  if (flag > 0) then
     call errore('tddft_cgsolver', 'cgsolver cannot achieve convergence', flag)
     stop
  end if
  call stop_clock ('cgsolver')
  
  return
  
END SUBROUTINE tddft_cgsolver


