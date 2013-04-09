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
  integer, parameter :: neq=3
  integer, parameter :: order=9
  integer, parameter :: ntimes=11
  real*8,  parameter :: t0 = 0.d0
  real*8,  parameter :: tf = 0.05d0
  real*8  :: rtol = 1.d-6
  real*8  :: lam = -100.d0
  real*8  :: y0(neq), atol(neq), ytrue(neq)
  real*8  :: times(ntimes) 
  real*8  :: yout(neq,ntimes)
  real*8  :: data(1), expDt(3,3)
  real*8, dimension(3,3) :: V, Vinv, D
  integer :: i, ierr


  !------- Internals ---------

  ! set initial conditions, tolerances, data
  y0 = 1.d0
  atol = 1.d-10
  data(1) = lam

  ! set problem-defining matrices
  V = reshape( (/1.d0,-1.d0,0.d0,-1.d0,2.d0,-1.d0,1.d0,1.d0,2.d0/), (/3,3/) )
  Vinv = 0.25d0 * reshape( (/5.d0,2.d0,1.d0,1.d0,2.d0,1.d0,-3.d0,-2.d0,1.d0/), (/3,3/) )
  D = reshape( (/-0.5d0,0.d0,0.d0,0.d0,-0.1d0,0.d0,0.d0,0.d0,lam/), (/3,3/) )

  ! set up output time vector
  do i=1,ntimes
     times(i) = t0 + (i-1)*(tf-t0)/(ntimes-1)
  end do

  ! call IRKSolve to compute the solution
  call IRKSolve(y0, neq, times, ntimes, rtol, atol, order, yout, data)

  ! check solution error
  do i=1,ntimes
     expDt = 0.d0
     expDt(1,1) = exp(times(i) * D(1,1))
     expDt(2,2) = exp(times(i) * D(2,2))
     expDt(3,3) = exp(times(i) * D(3,3))
     ytrue = matmul(V, matmul(expDt, matmul(Vinv, y0)))
     print '(A,es8.1,A,3(es12.4),A,3(es12.4),A,es8.1)', &
          '  t = ',times(i), &
          ',  ytrue =',ytrue, &
          ',  y =',yout(:,i), &
          ',  error =',maxval(abs(yout(:,i)-ytrue))
  enddo

end program driver
!-----------------------------------------------------------------


subroutine f(t, y, ydot, data)
  !-----------------------------------------------------------------
  ! Description: rhs function
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  real*8, intent(in)  :: t, y(3)
  real*8, intent(out) :: ydot(3)
  real*8 :: data(1)

  ! local variables
  real*8, dimension(3,3) :: V, Vinv, D
  real*8 :: lam

  !======= Internals =========
  lam = data(1)
  V = reshape( (/1.d0,-1.d0,0.d0,-1.d0,2.d0,-1.d0,1.d0,1.d0,2.d0/), (/3,3/) )
  Vinv = 0.25d0 * reshape( (/5.d0,2.d0,1.d0,1.d0,2.d0,1.d0,-3.d0,-2.d0,1.d0/), (/3,3/) )
  D = reshape( (/-0.5d0,0.d0,0.d0,0.d0,-0.1d0,0.d0,0.d0,0.d0,lam/), (/3,3/) )

  ydot = matmul(V, matmul(D, matmul(Vinv, y)))

  return

end subroutine f
!-----------------------------------------------------------------


subroutine J(t, y, Jac, data)
  !-----------------------------------------------------------------
  ! Description: Jacobian function
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  real*8, intent(in)  :: t, y(3)
  real*8, intent(out) :: Jac(3,3)
  real*8  :: data(1)

  ! local variables
  real*8, dimension(3,3) :: V, Vinv, D
  real*8 :: lam

  !======= Internals =========
  lam = data(1)
  V = reshape( (/1.d0,-1.d0,0.d0,-1.d0,2.d0,-1.d0,1.d0,1.d0,2.d0/), (/3,3/) )
  Vinv = 0.25d0 * reshape( (/5.d0,2.d0,1.d0,1.d0,2.d0,1.d0,-3.d0,-2.d0,1.d0/), (/3,3/) )
  D = reshape( (/-0.5d0,0.d0,0.d0,0.d0,-0.1d0,0.d0,0.d0,0.d0,lam/), (/3,3/) )

  Jac = matmul(V, matmul(D, Vinv))

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
