!The numerical solution to the harmonic oscillator
!with Hamiltonian H=(1/2)(p^2+w^2x^2) where w=1
!yields a coupled system of differential equations
!solved via 4th order Runge-Kutta Method
!~~Joe Wimmergren~~

!-------------------------------------------------------------
!to personalize, edit the functions, h, and initial conditions
!-------------------------------------------------------------




!-------------------------------------------------------------
!function for dx/dt
function functionx(t_,x_,p_) result(out)
    real, intent(in) :: p_,t_,x_    ! input
    real             :: out         ! output
    out=p_                          ! dx/dt=out
 end function functionx

!function for dp/dt
 function functionp(t_,x_,p_) result(out)
     real, intent(in) :: p_,t_,x_  ! input
     real            :: out        ! output
     out=-x_                      ! dp/dt=out
  end function functionp
!-------------------------------------------------------------
program rk
implicit none

real, allocatable :: x(:),t(:), p(:)
real :: k1p, k2p, k3p, k4p, h, t0, x0, tf, k1x, k2x, k3x, k4x, p0, functionx, functionp
integer :: i, k


parameter (h=0.001)             !step size
parameter (t0=0)                !starting time
parameter (tf=1)               !ending time
parameter (k= ((tf - t0)/h)+1)  !# of outter loops (rows)
parameter (x0=0)                !initial condition for x
parameter (p0=1)                !initial condition for p

open (unit=14, file='rk.txt', status='unknown')

!start arrays at 0
allocate (x(0:k))
allocate (p(0:k))
allocate (t(0:k))


!initialize starting conditions
x(0)=x0
p(0)=p0
t(0)=t0



do i=0, k-1 !loop for k's
  k1x=h*(functionx(t(i),x(i),p(i)))
  k1p=h*(functionp(t(i),x(i),p(i)))
  k2x=h*(functionx(t(i)+h/2,x(i)+k1x/2,p(i)+k1p/2))
  k2p=h*(functionp(t(i)+h/2,x(i)+k1x/2,p(i)+k1p/2))
  k3x=h*(functionx(t(i)+h/2,x(i)+k2x/2,p(i)+k2p/2))
  k3p=h*(functionp(t(i)+h/2,x(i)+k2x/2,p(i)+k2p/2))
  k4x=h*(functionx(t(i)+h,x(i)+k3x,p(i)+k3p))
  k4p=h*(functionp(t(i)+h,x(i)+k3x,p(i)+k3p))

  x(i+1)=x(i)+(k1x+2*k2x+2*k3x+k4x)/6
  p(i+1)=p(i)+(k1p+2*k2p+2*k3p+k4p)/6
  t(i+1)=t(i)+h
  write (14,*) x(i), p(i)
end do
  write (14,*) x(k), p(k)

deallocate (x)
deallocate (p)
deallocate (t)

close (14)

end program rk
