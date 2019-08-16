!function
!-----------------------------------------------------
!x'=b
function ad(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	out=b
end function ad


function bd(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	out=-muk*(b-f*sin(e)*sin(g)+h*cos(e)*cos(g))
end function bd


function cd(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	out=d
end function cd


function dd(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	out=-muk*(d+f*sin(e)*cos(g)+h*cos(e)*sin(g))
end function dd


function ed(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	out=f
end function ed


function fd(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	double precision             :: newfk1
	newfk1=muk*(b*sin(g)-d*cos(g)-f*sin(e))
	out=(((4*f*f-h*h)*cos(e)-2*h*j)*sin(e)-4*9.8*cos(e)+4*newfk1*sin(e))/(3+2*cos(2*e))
end function fd



function gd(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	out=h
end function gd


function hd(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	out=2*f*j/sin(e)
end function hd


function id(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	out=j
end function id


function jd(a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta) result(out)
	double precision, intent(in) :: a,b,c,d,e,f,g,h,i,j,muk,fk1,fk2,beta
	double precision 		 :: out
	double precision 		 :: newfk2
	newfk2=muk*(b*cos(g)+d*sin(g)+h*cos(e))
	out=h*f*sin(e)-(2*f*j*cos(e)/sin(e))-2*(newfk2+beta*j)
end function jd




!-----------------------------------------------------

program rung
implicit none

double precision, allocatable :: a(:),b(:),c(:),d(:),e(:),f(:),g(:),h(:),i(:),j(:),t(:)      !arrays
double precision :: muk,fk1,fk2,beta                             !constants
double precision :: ad,bd,cd,dd,ed,fd,gd,hd,id,jd                !functions
double precision :: a0,b0,c0,d0,e0,f0,g0,h0,i0,j0,tf,t0,hs       !initial conditions
double precision :: k1a,k1b,k1c,k1d,k1e,k1f,k1g,k1h,k1i,k1j      !k1's
double precision :: k2a,k2b,k2c,k2d,k2e,k2f,k2g,k2h,k2i,k2j      !k2's
double precision :: k3a,k3b,k3c,k3d,k3e,k3f,k3g,k3h,k3i,k3j      !k3's
double precision :: k4a,k4b,k4c,k4d,k4e,k4f,k4g,k4h,k4i,k4j      !k4's
integer :: n,k


!initial conditions
parameter (a0=7.6552)
parameter (b0=0)
parameter (c0=6.551)
parameter (d0=0)
parameter (e0=0.174533) !1.48353 radians
parameter (f0=0)
parameter (g0=0.79788)
parameter (h0=6.0141)
parameter (i0=5.3352)
parameter (j0=1.74)
parameter (beta=100)
parameter (muk=1000)
parameter (hs=0.00001)
parameter (t0=0)
parameter (tf=2)
parameter (k=((tf-t0)/hs)+1)


!solve for fk2

open (unit=12, file='./row1/85_100.txt', status='unknown')  !open output file


!allocate arrays
allocate (a(0:k))
allocate (b(0:k))
allocate (c(0:k))
allocate (d(0:k))
allocate (e(0:k))
allocate (f(0:k))
allocate (g(0:k))
allocate (h(0:k))
allocate (i(0:k))
allocate (j(0:k))
allocate (t(0:k))

a(0)=a0
b(0)=b0
c(0)=c0
d(0)=d0
e(0)=e0
f(0)=f0
g(0)=g0
h(0)=h0
i(0)=i0
j(0)=j0
t(0)=t0

do n=0, k-1 !loop for k's

	!k1

	k1a=hs*ad(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1b=hs*bd(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1c=hs*cd(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1d=hs*dd(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1e=hs*ed(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1f=hs*fd(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1g=hs*gd(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1h=hs*hd(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1i=hs*id(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)
	k1j=hs*jd(a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n),muk,fk1,fk2,beta)

	!k2
	k2a=hs*ad(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	k2b=hs*bd(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	k2c=hs*cd(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	k2d=hs*dd(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)
	
	k2e=hs*ed(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	k2f=hs*fd(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	k2g=hs*gd(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	k2h=hs*hd(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	k2i=hs*id(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	k2j=hs*jd(a(n)+k1a/2,b(n)+k1b/2,c(n)+k1c/2,d(n)+k1d/2,e(n)+k1e/2,f(n)+k1f/2,&
	g(n)+k1g/2,h(n)+k1h/2,i(n)+k1i/2,j(n)+k1j/2,muk,fk1,fk2,beta)

	!k3
	k3a=hs*ad(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3b=hs*bd(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3c=hs*cd(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3d=hs*dd(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3e=hs*ed(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3f=hs*fd(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3g=hs*gd(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3h=hs*hd(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3i=hs*id(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	k3j=hs*jd(a(n)+k2a/2,b(n)+k2b/2,c(n)+k2c/2,d(n)+k2d/2,e(n)+k2e/2,f(n)+k2f/2,&
	g(n)+k2g/2,h(n)+k2h/2,i(n)+k2i/2,j(n)+k2j/2,muk,fk1,fk2,beta)

	!k4
	k4a=hs*ad(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4b=hs*bd(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4c=hs*cd(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4d=hs*dd(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4e=hs*ed(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4f=hs*fd(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4g=hs*gd(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4h=hs*hd(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4i=hs*id(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	k4j=hs*jd(a(n)+k3a,b(n)+k3b,c(n)+k3c,d(n)+k3d,e(n)+k3e,f(n)+k3f,&
	g(n)+k3g,h(n)+k3h,i(n)+k3i,j(n)+k3j,muk,fk1,fk2,beta)

	!statements

	a(n+1)=a(n)+(k1a+2*k2a+2*k3a+k4a)/6
	b(n+1)=b(n)+(k1b+2*k2b+2*k3b+k4b)/6
	c(n+1)=c(n)+(k1c+2*k2c+2*k3c+k4c)/6
	d(n+1)=d(n)+(k1d+2*k2d+2*k3d+k4d)/6
	e(n+1)=e(n)+(k1e+2*k2e+2*k3e+k4e)/6
	f(n+1)=f(n)+(k1f+2*k2f+2*k3f+k4f)/6	
	g(n+1)=g(n)+(k1g+2*k2g+2*k3g+k4g)/6	
	h(n+1)=h(n)+(k1h+2*k2h+2*k3h+k4h)/6
	i(n+1)=i(n)+(k1i+2*k2i+2*k3i+k4i)/6
	j(n+1)=j(n)+(k1j+2*k2j+2*k3j+k4j)/6

	t(n+1)=t(n)+hs

!	write (12,*) t(n),a(n),b(n),c(n),d(n),e(n),f(n),g(n),h(n),i(n),j(n) 
	write (12,*) cos(e(n))*sin(g(n)),-1*cos(e(n))*cos(g(n))
!	if (e(n).lt.0) then
!		write (13,*) '0 point at ',t(n-1)
!		exit
!	end if

end do
	
!	write (12,*) t(k),a(k),b(k),c(k),d(k),e(k),f(k),g(k),h(k),i(k),j(k)
!	write (12,*) t(k),e(k)


close (12)
!close (13)

deallocate (a)
deallocate (b)
deallocate (c)
deallocate (d)
deallocate (e)
deallocate (f)
deallocate (g)
deallocate (h)
deallocate (i)
deallocate (j)

end program rung



