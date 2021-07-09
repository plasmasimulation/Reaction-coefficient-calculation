module modulesimpson
contains
    
function simpson_formula(a,b,s,T)
implicit none
Real(8),parameter :: PI = 3.141592653589793238D0
Real(8),parameter :: exp = 2.71828d0
Real(8),parameter :: me = 9.109d-31
Real(8),parameter :: e = 1.602d-19
Real(8),parameter :: kb = 1.3806488d-23
integer(4) :: n,m,p,q
Real(8),intent(in):: a,b,s,T
real(8) :: h,f,simpson_formula
real :: x(1000001),y(1000001)

h=1.d1
n=(b-a)/h
m=n/2
f=0.d0

do p=1,(n+1)
    x(p)=a+(p-1)*h
    y(p)=s*x(p)*exp**(-me*x(p)*x(p)/(2*T*11605*kb))*4*PI*x(p)*x(p)
end do

do q=1,m
    f=f+2*y(2*q-1)+4*y(2*q)
end do
simpson_formula=(f-y(1)+y(n+1))*h/3
return
end function simpson_formula
end module modulesimpson