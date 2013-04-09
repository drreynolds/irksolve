!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
!
! Copyright 2013
! All rights reserved
!=================================================================

!-----------------------------------------------------------------
! Driver to test utility functions from IRKSolve.F90
!-----------------------------------------------------------------
program driver

  ! Declarations
  implicit none

  ! general problem variables
  integer, parameter :: n=10
  integer, parameter :: neq=2
  integer, parameter :: s=3
  real*8  :: A(n,n), x(n), xtrue(n), b(n)
  real*8  :: y(neq), y2(neq), z(neq,s), ff(neq,s), Amat(neq*s,neq*s), data(5)
  real*8  :: h, t
  integer :: i, ierr


  !------- Internals ---------

  ! test 1: GE_BS with nonsingular matrix
  print *, '  '
  call random_number(A)
  call random_number(xtrue)
  do i=1,10
     A(i,i) = A(i,i) + 10.d0
  end do
  b = matmul(A,xtrue)
  call GE_BS(A, x, b, 10, ierr)
  print *, 'GE_BS on nonsingular problem returned ierr =',ierr
  print *, '   ||x-xtrue||_inf =',maxval(abs(x-xtrue))

  ! test 2: GE_BS with nonsingular matrix requiring pivoting
  print *, '  '
  call random_number(A)
  call random_number(xtrue)
  do i=1,10
     A(i,i) = A(i,i) + 10.d0
  end do
  A(1,2) = 10.d0
  A(2,1) = 10.d0
  A(1,1) = 0.d0
  A(2,2) = 0.d0
  b = matmul(A,xtrue)
  call GE_BS(A, x, b, 10, ierr)
  print *, 'GE_BS on nonsingular problem requiring pivoting returned ierr =',ierr
  print *, '   ||x-xtrue||_inf =',maxval(abs(x-xtrue))

  ! test 3: GE_BS with singular
  print *, '  '
  call random_number(A)
  call random_number(xtrue)
  do i=1,10
     A(i,i) = A(i,i) + 10.d0
  end do
  A(1,2) = 10.d0
  A(:,9) = A(:,3)
  b = matmul(A,xtrue)
  call GE_BS(A, x, b, 10, ierr)
  print *, 'GE_BS on singular problem returned ierr =',ierr
  print *, '   ||x-xtrue||_inf =',maxval(abs(x-xtrue))

  ! test 4: F_IRK test
  print *, '  '
  z = 1.d0
  t = 0.d0
  h = 1.d0
  y = 1.d0
  call F_IRK(z, t, h, y, ff, 2, 3, data, ierr)
  print *, 'F_IRK with f(t,y)=-y, yold = 1, h=1, z=1 gives ierr =',ierr
  print *, '   ff =',ff

  ! test 5: A_IRK test
  print *, '  '
  z = 1.d0
  t = 0.d0
  h = 1.d0
  call A_IRK(z, t, h, Amat, 2, 3, data, ierr)
  print *, 'A_IRK with f(t,y)=-y, h=1, z=1 gives ierr =',ierr
  print *, '   Amat ='
  do i=1,6
     print '(6(f8.5,1x))', Amat(i,:)
  end do

  ! test 6: Y_IRK test
  print *, '  '
  z = 1.d0
  t = 0.d0
  h = 1.d0
  y = 2.d0
  call Y_IRK(z, t, h, y, y2, 2, 3, data, ierr)
  print *, 'Y_IRK with f(t,y)=-y, yold = 2, h=1, z=1 gives ierr =',ierr
  print *, '   y2 =',y2
  
  ! test 7: Newton test
  print *, '  '
  z = 1.d0
  t = 0.d0
  h = 1.d0
  y = 1.d0
  call Newton(z, t, h, y, 2, 3, 1.d-5, 20, data, ierr)
  print *, 'Newton with f(t,y)=-y, yold = 1, h=1, z=1 gives ierr =',ierr
  print *, '   z =',z
  print *, '  '

end program driver
!-----------------------------------------------------------------


subroutine f(t, y, ydot, data)
  !-----------------------------------------------------------------
  ! Description: rhs function
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  real*8, intent(in)  :: t, y(2)
  real*8, intent(out) :: ydot(2)
  real*8 :: data(5)

  !======= Internals =========
  ydot = -y

  return

end subroutine f
!-----------------------------------------------------------------


subroutine J(t, y, Jac, data)
  !-----------------------------------------------------------------
  ! Description: Jacobian function
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  real*8, intent(in)  :: t, y(2)
  real*8, intent(out) :: Jac(2,2)
  real*8  :: data(5)
  integer :: i

  !======= Internals =========
  Jac = 0.d0
  do i=1,2
     Jac(i,i) = -1.d0
  enddo

  return

end subroutine J
!-----------------------------------------------------------------


subroutine update_data(data)
  !-----------------------------------------------------------------
  ! Description: data update function
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  real*8  :: data(5)

  !======= Internals =========
  call random_number(data)

  return

end subroutine update_data
!-----------------------------------------------------------------
