! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - An example to test MINPACK routine hybrd
!     - The example is from IMSL
! Author:
!     Xin Tang @ IMF, Summer 2019
!  
! Record of Revisions:
!         Date:                 Description of Changes
!     ===========        =================================
!      07/12/2019:                 Original Code
!
! Compiling Environment:
!   GNU gfortran on WSL1
!
! Library Used:
!   MINPACK with source code
! =============================================================================

program minpack_ex

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: np = 3
    integer, parameter :: lwa = (np*(3*np+13))/2
    integer :: iflag, info
    real(dp), dimension(np) :: x, fvec
    real(dp), dimension(lwa) :: wa
    real(dp) :: fnorm, tol

    external :: fcn

    x = (/ 4.0_dp, 4.0_dp, 4.0_dp /)
    fvec = 0.0_dp
    fnorm = 10000.0_dp
    tol = 1.0d-7

    ! call fcn(np,x,fvec,iflag)
    call hybrd1(fcn,np,x,fvec,tol,info,wa,lwa)

    fnorm = sqrt(sum(fvec**2,1))

    write (*,*) 'fnorm = ', fnorm
    write (*,*) 'info = ', info
    write (*,*) 'x = ', x

end program minpack_ex

subroutine fcn(np,x,fvec,iflag)

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, intent(in) :: np
    integer, intent(out) :: iflag
    real(dp), dimension(np), intent(in) :: x
    real(dp), dimension(np), intent(out) :: fvec

    fvec(1) = x(1) + exp(x(1)-1.0_dp) + (x(2)+x(3))**2 - 27.0_dp
    fvec(2) = exp(x(2)-2.0_dp)/x(1) + x(3)**2 - 10.0_dp 
    fvec(3) = x(3) + sin(x(2)-2.0_dp) + x(2)**2 - 7.0_dp
    iflag = 0

end subroutine fcn