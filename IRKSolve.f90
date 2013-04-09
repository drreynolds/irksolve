!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
!
! Copyright 2013
! All rights reserved
!=================================================================

subroutine IRKSolve(y0, neq, times, ntimes, rtol, atol, q, y, data)
  !-----------------------------------------------------------------
  ! Description: time-stepping driver/solver for Radau IIA IRK 
  ! methods.
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  integer, intent(in)  :: neq, ntimes, q
  real*8,  intent(in)  :: y0(neq), times(ntimes), rtol, atol(neq)
  real*8,  intent(out) :: y(neq,ntimes)
  real*8               :: data(*)

  !======= Local variables =========
  real*8, allocatable :: z(:,:), c(:)
  real*8, dimension(neq) :: ynew, yold, y1, y2, yerr, ewt
  integer :: ierr, s, i, tstep, st_fail, newt_maxit
  real*8  :: t, t_this, c1, c2, hmin, h, h_this, newt_tol, tstage
  real*8  :: eps, err_step, safety, dt_growth, alpha
  real*8  :: ONEMSM, ONEPSM, ERRTOL, h_reduce


  !======= Internals =========

  ! check inputs
  if (mod(q,2) == 0) then
     print *, 'IRKSolve error: order must be one of {3,5,9}'
     return
  end if

  if (neq < 1) then
     print *, 'IRKSolve error: problem must have at least one equation'
     return
  end if

  if (ntimes < 2) then
     print *, 'IRKSolve error: problem must a multiple time outputs'
     return
  end if


  ! set solver parameters
  c1 = -1.d0/(2**q - 1.d0)     ! Richardson extrapolation factor for h step
  c2 = 1.d0 - c1               ! Richardson extrapolation factor for h/2 step
  newt_maxit = 20              ! max number of Newton iterations
  newt_tol   = 1.d-10          ! Newton solver tolerance
  h_reduce   = 0.1d0           ! failed step reduction factor 
  ONEMSM     = 1.d0-sqrt(eps)  ! coefficients to account for
  ONEPSM     = 1.d0+sqrt(eps)  !     floating-point roundoff
  ERRTOL     = 1.1d0           ! upper bound on allowed step error (WRMS norm)
  eps        = 2.d-16          ! approximation of roudoff
  safety     = 0.9d0
  dt_growth  = 10.d0
  alpha      = 1.d0/q

  ! set up stage structure
  s = (q+1)/2
  allocate(z(neq,s),c(s))

  ! set up Butcher table times
  if (s == 2) then
     c = (/ 1.d0/3.d0, 1.d0 /)
  elseif (s == 3) then
     c = (/ (4.d0-sqrt(6.d0))/10.d0, (4.d0+sqrt(6.d0))/10.d0, 1.d0 /)
  else 
     c = (/ 0.05710419611451768d0, 0.2768430136381238d0, &
          0.5835904323689168d0, 0.8602401356562195d0, 1.d0 /)
  end if

  ! set initial condition into solution, output
  ynew = y0
  y(:,1) = y0
  t = times(1)

  ! set time-stepping parameters
  hmin = 1.d-5 * minval(times(2:ntimes)-times(1:ntimes-1))
  h = hmin

  ! initialize error weight vector
  ewt = 1.d0 / (rtol*ynew + atol)

  ! outer time-stepping loop (between output times)
  do tstep = 2,ntimes

     do while (t < times(tstep)*ONEMSM)

        ! bound internal time step 
        h = max(h, hmin)            ! enforce minimum time step size
        h = min(h, times(tstep)-t)  ! stop at output time

        ! update "old" state
        yold = ynew

        ! reset solve failure flag
        st_fail = 0

        ! update user data structure for time step
        call update_data(data)


        !----- solve with time step h -----
        h_this = h

        ! set Newton initial guesses as previous step solution
        z = 0.d0
        do i = 1,s
           z(:,i) = ynew
        end do

        ! call Newton solver to update solution in time
        call Newton(z, t, h_this, yold, neq, s, newt_tol, newt_maxit, data, ierr)
        if (ierr /= 0)  st_fail = 1

        ! compute solution with this h
        call Y_IRK(z, t, h_this, yold, y1, neq, s, data, ierr)

        !----- solve with steps of size h/2 -----
        if (st_fail == 0) then

           h_this = 0.5*h
           yold = ynew

           ! set Newton initial guesses as linear interpolants of full step solutions
           z = 0.d0
           do i = 1,s
              tstage = c(i)*0.5
              z(:,i) = (1.d0-tstage)*ynew + tstage*y1
           end do

           ! call Newton solver to update solution in time
           call Newton(z, t, h_this, yold, neq, s, newt_tol, newt_maxit, data, ierr)
           if (ierr /= 0)  st_fail = 1

           ! compute half-step solution
           call Y_IRK(z, t, h_this, yold, y2, neq, s, data, ierr)

           ! if first half-step succeeded, take second half-step
           if (st_fail == 0) then

              yold = y2
              t_this = t + h_this

              ! set Newton initial guesses as linear interpolants of full step solutions
              z = 0.d0
              do i = 1,s
                 tstage = 0.5d0 + c(i)*0.5
                 z(:,i) = (1.d0-tstage)*ynew + tstage*y1
              end do

              ! call Newton solver to update solution in time
              call Newton(z, t_this, h_this, yold, neq, s, newt_tol, newt_maxit, data, ierr)
              if (ierr /= 0)  st_fail = 1

              ! compute full-step solution (store back in y2)
              call Y_IRK(z, t_this, h_this, yold, y2, neq, s, data, ierr)

           end if  ! end second half-step

        end if ! half-step solutions


        ! if solves succeeded, check step accuracy
        if (st_fail == 0) then

           ! estimate error in current step
           yerr = y1-y2
           err_step = max( maxval(abs(Yerr*ewt)), eps)

           ! if error too high, flag step as a failure (will be be recomputed)
           if (err_step > ERRTOL*ONEPSM) then
              st_fail = 1
           end if

        end if


        ! if step was successful (solves succeeded, and error acceptable)
        if (st_fail == 0) then

           ! update time for last successful step
           t = t + h

           ! update solution using Richardson extrapolation
           ynew = c1*y1 + c2*y2

           ! update error weight vector
           ewt = 1.d0 / (rtol*ynew + atol)

           ! estimate error and update time step
           h_this = safety * h * err_step**(-alpha)
           h = min(dt_growth*h, h_this)

        ! if step solves or error test failed
        else

           ! if already at minimum step, just return with failure
           if (h <= hmin) then
              print *, 'Cannot achieve desired accuracy.'
              print *, 'Consider reducing hmin or increasing rtol.'
              ierr = 1
              return
           end if

           ! otherwise, reset guess, reduce time step, retry solve
           h = h * h_reduce

        end if  ! end logic tests for step success/failure

     end do  ! end while loop attempting to solve steps to next output time

     ! store updated solution in output array
     y(:,tstep) = ynew

  end do

  ! clean up
  deallocate(z,c)

end subroutine IRKSolve
!=================================================================



subroutine Newton(z, t, h, yold, neq, s, tol, maxit, data, ierr)
  !-----------------------------------------------------------------
  ! Description: Newton solver
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  integer, intent(in)  :: neq, s, maxit
  real*8,  intent(in)  :: t, h, yold(neq), tol
  integer, intent(out) :: ierr
  real*8               :: z(neq,s), data(*)

  !======= Local variables =========
  real*8, dimension(s*neq)       :: y, delt, F
  real*8, dimension(s*neq,s*neq) :: A
  integer :: inewt


  !======= Internals =========

  ! reshape input to form y vector, initialize update to ones
  y = reshape( source=z, shape=(/s*neq/) )
  delt = 1.d0

  ! perform Newton iterations
  do inewt=1,maxit

     ! check for success
     if (maxval(abs(delt)) < tol) then
        ierr = 0
        return
     end if

     ! compute the residual and Jacobian
     call F_IRK(z, t, h, yold, F, neq, s, data, ierr)
     if (ierr /= 0) return
     call A_IRK(z, t, h, A, neq, s, data, ierr)
     if (ierr /= 0) return
     
     ! solve for Newton update
     call GESPP_BS(A, delt, F, s*neq, ierr)
     if (ierr /= 0) return

     ! perform update, and copy results into z
     y = y - delt
     z = reshape( source=y, shape=(/neq,s/) )

  end do

  ! if we've made it this far, then Newton did not converge
  ierr = 1
  return

end subroutine Newton
!=================================================================



subroutine GESPP_BS(A, x, b, n, ierr)
  !-----------------------------------------------------------------
  ! Description: Gaussian Elimination with scaled partial pivoting 
  ! routine; destroys both A and b
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  integer, intent(in)  :: n
  integer, intent(out) :: ierr
  real*8  :: A(n,n), x(n), b(n)

  !======= Local variables =========
  integer :: i, j, k, p, parr(1)
  real*8, dimension(n) :: s, tmp
  real*8 :: tmp2, m

  !======= Internals =========

  ! set up scale factor array
  do k=1,n
     s(k) = maxval(abs(A(i,:)))
  end do

  ! outer Gaussian elimination loop over pivot columns
  do k=1,n-1

     ! find the pivot row p
     parr = maxloc( abs(A(k:n,k))/s(k:n) )
     p = parr(1)+k-1

     ! swap rows in A
     tmp = A(p,:)
     A(p,:) = A(k,:)
     A(k,:) = tmp

     ! swap rows in b
     tmp2 = b(p)
     b(p) = b(k)
     b(k) = tmp2
   
     ! swap rows in s
     tmp2 = s(p)
     s(p) = s(k)
     s(k) = tmp2

     ! check for failure
     if (abs(A(k,k)) < 1.d-13) then
        print *, 'GESPP_BS error: singular matrix!'
        ierr = 1
        return
     end if

     ! perform elimination
     m = 1.d0/A(k,k)
     do j = k+1,n
        A(k+1:n,j) = A(k+1:n,j) - m*A(k,j)*A(k+1:n,k)
     end do
     b(k+1:n) = b(k+1:n) - m*b(k)*A(k+1:n,k)

  end do

  ! check for singularity in last row/column
  if (abs(A(n,n)) < 1.d-13) then
     print *, 'GESPP_BS error: singular matrix!'
     ierr = 1
     return
  end if

  ! perform column-oriented back-substitution 
  do j=n,1,-1
     x(j) = b(j) / A(j,j)
     do i=1,j-1
        b(i) = b(i) - A(i,j)*x(j)
     end do
  end do

  ierr = 0
  return

end subroutine GESPP_BS
!=================================================================



subroutine F_IRK(z, t, h, yold, ff, neq, s, data, ierr)
  !-----------------------------------------------------------------
  ! Description: nonlinear residual vector
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  integer, intent(in)  :: neq, s
  real*8,  intent(in)  :: z(neq,s), t, h, yold(neq)
  real*8,  intent(out) :: ff(neq,s)
  integer, intent(out) :: ierr
  real*8 :: data(*)

  !======= Local variables =========
  real*8, allocatable :: y(:), frhs(:), A(:,:), c(:), allf(:,:)
  real*8 :: tstage
  integer :: is, js


  !======= Internals =========

  allocate(y(neq), frhs(neq), A(s,s), c(s), allf(neq,s))

  ! set up Butcher table
  A = 0.d0
  c = 0.d0
  if (s == 2) then
     A(:,1) = (/ 5.d0/12.d0, 9.d0/12.d0 /)
     A(:,2) = (/ -1.d0/12.d0, 3.d0/12.d0 /)
     c = (/ 1.d0/3.d0, 1.d0 /)
  elseif (s == 3) then
     A(:,1) = (/ (88.d0-7.d0*sqrt(6.d0))/360.d0, &
          (296.d0+169.d0*sqrt(6.d0))/1800.d0, &
          (16.d0-sqrt(6.d0))/36.d0 /)
     A(:,2) = (/ (296.d0-169.d0*sqrt(6.d0))/1800.d0, &
          (88.d0+7.d0*sqrt(6.d0))/360.d0, &
          (16.d0+sqrt(6.d0))/36.d0 /)
     A(:,3) = (/ (-2.d0+3.d0*sqrt(6.d0))/225.d0, &
          (-2.d0-3.d0*sqrt(6.d0))/225.d0, &
          1.d0/9.d0 /)
     c = (/ (4.d0-sqrt(6.d0))/10.d0, (4.d0+sqrt(6.d0))/10.d0, 1.d0 /)
  else 
     A(:,1) = (/ 0.07299886431790337d0, &
          0.1537752314791824d0, &
          0.1400630456848099d0, &
          0.1448943081095342d0, &
          0.1437135607912259d0 /)
     A(:,2) = (/ -0.02673533110794565d0, &
          0.1462148678474935d0, &
          0.2989671294912833d0, &
          0.2765000687601608d0, &
          0.2813560151494621d0 /)
     A(:,3) = (/ 0.01867692976398445d0, &
          -0.03644456890512816d0, &
          0.1675850701352492d0, &
          0.3257979229104191d0, &
          0.3118265229757413d0 /)
     A(:,4) = (/ -0.01287910609330652d0, &
          0.02123306311930480d0, &
          -0.03396910168661794d0, &
          0.1287567532549115d0, &
          0.2231039010835707d0 /)
     A(:,5) = (/ 0.005042839233882052d0, &
          -0.007935579902728813d0, &
          0.01094428874419233d0, &
          -0.01570891737880607d0, &
          0.04d0 /)
     c = (/ 0.05710419611451768d0, 0.2768430136381238d0, &
          0.5835904323689168d0, 0.8602401356562195d0, 1.d0 /)
  end if

  ! call f at our guesses
  do is = 1,s
     tstage = t + h*c(is)
     y = z(:,is)
     call f(tstage, y, frhs, data)
     allf(:,is) = frhs
  end do

  ! form the IRK residuals
  !    Fs = zs - y_n - h*sum(a(s,j)*fj)
  ff = 0.d0
  do is = 1,s
     ff(:,is) = z(:,is) - yold
     do js = 1,s
        ff(:,is) = ff(:,is) - h*A(is,js)*allf(:,js)
     end do
  end do

  deallocate(y, frhs, A, c, allf)

  ierr = 0
  return

end subroutine F_IRK
!=================================================================



subroutine A_IRK(z, t, h, Amat, neq, s, data, ierr)
  !-----------------------------------------------------------------
  ! Description: nonlinear residual Jacobian
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  integer, intent(in)  :: neq, s
  real*8,  intent(in)  :: z(neq,s), t, h
  real*8,  intent(out) :: Amat(neq*s,neq*s)
  integer, intent(out) :: ierr
  real*8 :: data(*)

  !======= Local variables =========
  real*8, allocatable :: y(:), Jac(:,:), A(:,:), c(:), allJ(:,:,:)
  real*8 :: tstage
  integer :: is, js, irow, jcol


  !======= Internals =========
  allocate(y(neq), Jac(neq,neq), A(s,s), c(s), allJ(neq,neq,s))

  ! set up Butcher table
  A = 0.d0
  c = 0.d0
  if (s == 2) then
     A(:,1) = (/ 5.d0/12.d0, 9.d0/12.d0 /)
     A(:,2) = (/ -1.d0/12.d0, 3.d0/12.d0 /)
     c = (/ 1.d0/3.d0, 1.d0 /)
  elseif (s == 3) then
     A(:,1) = (/ (88.d0-7.d0*sqrt(6.d0))/360.d0, &
          (296.d0+169.d0*sqrt(6.d0))/1800.d0, &
          (16.d0-sqrt(6.d0))/36.d0 /)
     A(:,2) = (/ (296.d0-169.d0*sqrt(6.d0))/1800.d0, &
          (88.d0+7.d0*sqrt(6.d0))/360.d0, &
          (16.d0+sqrt(6.d0))/36.d0 /)
     A(:,3) = (/ (-2.d0+3.d0*sqrt(6.d0))/225.d0, &
          (-2.d0-3.d0*sqrt(6.d0))/225.d0, &
          1.d0/9.d0 /)
     c = (/ (4.d0-sqrt(6.d0))/10.d0, (4.d0+sqrt(6.d0))/10.d0, 1.d0 /)
  else 
     A(:,1) = (/ 0.07299886431790337d0, &
          0.1537752314791824d0, &
          0.1400630456848099d0, &
          0.1448943081095342d0, &
          0.1437135607912259d0 /)
     A(:,2) = (/ -0.02673533110794565d0, &
          0.1462148678474935d0, &
          0.2989671294912833d0, &
          0.2765000687601608d0, &
          0.2813560151494621d0 /)
     A(:,3) = (/ 0.01867692976398445d0, &
          -0.03644456890512816d0, &
          0.1675850701352492d0, &
          0.3257979229104191d0, &
          0.3118265229757413d0 /)
     A(:,4) = (/ -0.01287910609330652d0, &
          0.02123306311930480d0, &
          -0.03396910168661794d0, &
          0.1287567532549115d0, &
          0.2231039010835707d0 /)
     A(:,5) = (/ 0.005042839233882052d0, &
          -0.007935579902728813d0, &
          0.01094428874419233d0, &
          -0.01570891737880607d0, &
          0.04d0 /)
     c = (/ 0.05710419611451768d0, 0.2768430136381238d0, &
          0.5835904323689168d0, 0.8602401356562195d0, 1.d0 /)
  end if

  ! call J at our guesses
  do is = 1,s
     tstage = t + h*c(is)
     y = z(:,is)
     call J(tstage, y, Jac, data)
     allJ(:,:,is) = Jac
  end do

  ! form the IRK Jacobian
  Amat = 0.d0
  do js = 1,s
     do is = 1,s
        Amat(neq*(is-1)+1:neq*is,neq*(js-1)+1:neq*js) = A(is,js)*allJ(:,:,js)
     end do
  end do
  Amat = -h*Amat
  do jcol = 1,neq*s
     Amat(jcol,jcol) = Amat(jcol,jcol) + 1.d0
  end do

  deallocate(y, Jac, A, c, allJ)

  ierr = 0
  return

end subroutine A_IRK
!=================================================================



subroutine Y_IRK(z, t, h, yold, ynew, neq, s, data, ierr)
  !-----------------------------------------------------------------
  ! Description: solution construction routine
  !-----------------------------------------------------------------
  implicit none

  !======= Input variables =========
  integer, intent(in)  :: neq, s
  real*8,  intent(in)  :: z(neq,s), t, h, yold(neq)
  real*8,  intent(out) :: ynew(neq)
  integer, intent(out) :: ierr
  real*8               :: data(*)

  !======= Local variables =========
  real*8, dimension(neq) :: y, frhs
  real*8, allocatable :: b(:), c(:), allf(:,:)
  real*8  :: tstage
  integer :: is


  !======= Internals =========

  ! set up Butcher table
  allocate(c(s), b(s), allf(neq,s))
  if (s == 2) then
     b = (/ 3.d0/4.d0, 1.d0/4.d0 /)
     c = (/ 1.d0/3.d0, 1.d0 /)
  elseif (s == 3) then
     b = (/ (16.d0-sqrt(6.d0))/36.d0, (16.d0+sqrt(6.d0))/36.d0, 1.d0/9.d0 /)
     c = (/ (4.d0-sqrt(6.d0))/10.d0, (4.d0+sqrt(6.d0))/10.d0, 1.d0 /)
  else 
     b = (/ 0.1437135607912259d0, 0.2813560151494621d0, &
          0.3118265229757413d0, 0.2231039010835707d0, 0.04d0 /)
     c = (/ 0.05710419611451768d0, 0.2768430136381238d0, &
          0.5835904323689168d0, 0.8602401356562195d0, 1.d0 /)
  end if

  ! call f at our guesses
  do is = 1,s
     tstage = t + h*c(is)
     y = z(:,is)
     call f(tstage, y, frhs, data)
     allf(:,is) = frhs
  end do

  ! form the solution
  !    ynew = yold + h*sum(b(j)*fj)
  ynew = yold
  do is = 1,s
     ynew = ynew + h*b(is)*allf(:,is)
  end do

  deallocate(c, b, allf)

  ierr = 0
  return

end subroutine Y_IRK
!=================================================================
