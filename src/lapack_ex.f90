! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - An example to test LAPACK routine dgels
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
!   LAPACK, BLAS and ATLAS in binary form
! =============================================================================

program lapack_ex

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: m = 9, n = 4, nrhs = 1
    integer, parameter :: lda = 9, ldb = 9
    integer :: lwork = -1, info = 0
    real(dp), dimension(m,n) :: A
    real(dp), dimension(m) :: B
    real(dp), dimension(1) :: workquery
    real(dp), allocatable, dimension(:) :: work

    A = 1.0_dp
    A(:,2) = (/ 7.0_dp, 2.0_dp, 7.0_dp, -3.0_dp, 2.0_dp, 2.0_dp, -3.0_dp, 2.0_dp, 2.0_dp /)
    A(:,3) = (/ 5.0_dp, -1.0_dp, 3.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp /)       
    A(:,4) = (/ 6.0_dp, 6.0_dp, 5.0_dp, 4.0_dp, 0.0_dp, 7.0_dp, 3.0_dp, 1.0_dp, 4.0_dp /)        

    B = (/ 7.0_dp, -5.0_dp, 6.0_dp, 5.0_dp, 5.0_dp, -2.0_dp, 0.0_dp, 8.0_dp, 3.0_dp /)

    ! Query the optimal workspace
    call dgels('N',m,n,nrhs,A,lda,B,ldb,workquery,lwork,info)

    lwork = int(workquery(1))
    allocate(work(lwork))

    ! Do the OLS
    call dgels('N',m,n,nrhs,A,lda,B,ldb,work,lwork,info)

    ! Print the results
    write (*,*) 'OLS Estimates = ', B(1:n)

end program lapack_ex