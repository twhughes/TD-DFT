!-----------------------------------------------------------------------
MODULE ruku
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for the Runge-Kutta integration scheme
  !
  USE kinds,               ONLY : dp
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k, ecutwfc, et
  USE klist,                       ONLY : nks, nkstot, wk, xk, nelec, ngk
  USE buffers,          ONLY : close_buffer
  IMPLICIT NONE
  SAVE

!  complex(dp), allocatable :: k1(:,:),k2(:,:),k3(:,:),k4(:,:)
  complex(dp), allocatable :: k1(:,:)
  integer :: irk4
 
  integer, parameter :: iunruku1 = 53 ! to save ruku intermediate vectors
  integer, parameter :: iunruku2 = 54 ! to save ruku intermediate vectors
  integer, parameter :: iunruku3 = 55 ! to save ruku intermediate vectors
  integer, parameter :: iunruku4 = 56 ! to save ruku intermediate vectors
  
Contains
  
  subroutine allocate_ruku
     allocate(k1(npwx,nbnd))
!     allocate(k2(npwx,nbnd))
!     allocate(k3(npwx,nbnd))
!     allocate(k4(npwx,nbnd))
  end subroutine allocate_ruku

  subroutine deallocate_ruku
!     deallocate(k4,k3,k2,k1)
     deallocate(k1)
  end subroutine deallocate_ruku
 
  subroutine ruku_closefil
  call close_buffer( iunruku1, 'delete' )
  call close_buffer( iunruku2, 'delete' )
  call close_buffer( iunruku3, 'delete' )
  call close_buffer( iunruku4, 'delete' )
  end subroutine ruku_closefil
END MODULE ruku

