!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
!
! Copyright 2013
! All rights reserved
!=================================================================

!-----------------------------------------------------------------
! Driver to test IRKSolve on a simple problem with analytical 
! solution.
!-----------------------------------------------------------------
program driver

  ! Declarations
  implicit none

  ! general problem variables
  integer, parameter :: neq=1
  integer, parameter :: order=9
  integer, parameter :: ntimes=11
  real*8,  parameter :: t0 = 0.d0
  real*8,  parameter :: tf = 10.d0
  real*8  :: y0 = 0.d0
  real*8  :: rtol = 1.d-4
  real*8  :: atol = 1.d-10
  real*8  :: times(ntimes) 
  real*8  :: yout(neq,ntimes)
  real*8  :: data(1)
  real*8  :: ytrue
  integer :: i, ierr


  !------- Internals ---------

  ! set up output time vector
  do i=1,ntimes
     times(i) = t0 + (i-1)*(tf-t0)/(ntimes-1)
  end do

  ! call IRKSolve to compute the solution
  call IRKSolve(y0, neq, times, ntimes, rtol, atol, order, yout, data)

  ! check solution error
  do i=1,ntimes
     ytrue = log(((times(i)+1.d0)**2 + 1.d0)/2.d0)
     print '(A,es8.1,2(A,es12.5),A,es8.1)', &
          '  t = ',times(i), &
          ',  ytrue =',ytrue, &
          ',  y =',yout(1,i), &
          ',  error =',abs(yout(1,i)-ytrue)
  enddo

end program driver
!-----------------------------------------------------------------


subroutine f(t, y, ydot, data)
  !-----------------------------------------------------------------
  ! Description: rhs function
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  real*8, intent(in)  :: t, y
  real*8, intent(out) :: ydot
  real*8 :: data(1)

  !======= Internals =========
  ydot = (t + 1.d0)*exp(-y)

  return

end subroutine f
!-----------------------------------------------------------------


subroutine J(t, y, Jac, data)
  !-----------------------------------------------------------------
  ! Description: Jacobian function
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  real*8, intent(in)  :: t, y
  real*8, intent(out) :: Jac
  real*8  :: data(1)

  !======= Internals =========
  Jac = -(t + 1.d0)*exp(-y)

  return

end subroutine J
!-----------------------------------------------------------------


subroutine update_data(data)
  !-----------------------------------------------------------------
  ! Description: data update function
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  real*8  :: data(1)

  !======= Internals =========

  return

end subroutine update_data
!-----------------------------------------------------------------
