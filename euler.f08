!The numerical solution to the harmonic oscillator
!with Hamiltonian H=(1/2)(p^2+w^2x^2) where w=1
!yields a coupled system of differential equations
!solved via Euler's Method
!~~Joe Wimmergren~~

!function for dx/dt
function functionx(t_,x_,p_) result(out)
    real, intent(in) :: p_,t_,x_ ! input
    real             :: out ! output
    out=p_
 end function functionx

!function for dp/dt
 function functionp(t_,x_,p_) result(out)
     real, intent(in) :: p_,t_,x_ ! input
     real            :: out ! output
     out=-x_
  end function functionp



program EulerOscillator
implicit none

real, allocatable :: x(:), p(:), t(:)
real :: h, t0, tf, x0, p0, functionx, functionp
integer :: k, i


parameter (h=1)               !step size
parameter (t0=0)                !starting time
parameter (tf=1)                !ending time
parameter (k= ((tf - t0)/h)+1)  !# of outter loops (rows)
parameter (x0=0)                !initial condition for x
parameter (p0=1)                !initial condition for p


open (unit=13, file='euler.txt', status='unknown')

allocate (x(0:k+1))
allocate (p(0:k+1))
allocate (t(0:k+1))

!initialize t at t0
t(0)=t0

!initialize x at x0
x(0)=x0

!initialize y at y0
p(0)=p0
write(13,*)
do i=0, k-1
  x(i+1)=x(i)+h*(functionx(t(i),x(i),p(i)))
  p(i+1)=p(i)+h*(functionp(t(i),x(i),p(i)))
  write(13,*)  t(i), x(i),p(i)
  t(i+1)=t(i)+h
end do


deallocate (x)
deallocate (p)
deallocate (t)

close (13)

end program EulerOscillator
