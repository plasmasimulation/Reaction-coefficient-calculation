program reaction_rate_solution
    use modulesimpson
implicit none

Real(8),parameter :: PI = 3.141592653589793238D0
Real(8),parameter :: exp = 2.71828d0
Real(8),parameter :: me = 9.109d-31
Real(8),parameter :: e = 1.602d-19
Real(8),parameter :: kb = 1.3806488d-23
integer(4):: i=1
!real(8),demention(10001) t,y1,y2
real :: t(10001)
real(8) :: y1(10001),y2(10001),Te(10001),sec1(10001),sec2(10001),sec3(10001)
real(8) :: fun1(10001),fun2(10001),I1(10001),I2(10001),k1(10001),k2(10001)
!real,external :: simpson_formula

open(10,file="reaction rate.dat")
open(20,file="SN           1 .dat")


do i=1,10001
    t(i)=0.1*i
    y1(i)=2.336d-14*(t(i)**1.609)*exp**(0.0168*log(t(i))**2-0.1171*log(t(i))**3)
    y2(i)=2.34d-14*(t(i)**0.59)*exp**(-17.44/t(i))
    
    read(20,*) Te(i),sec1(i),sec2(i),sec3(i)
    fun1(i)=simpson_formula(0.d0,1.d7,sec1(i),Te(i))
    fun2(i)=simpson_formula(0.d0,1.d7,sec3(i),Te(i))
    
    k1(i)=(me/(2*PI*kb*Te(i)*11605))**1.5*fun1(i)
    k2(i)=(me/(2*PI*kb*Te(i)*11605))**1.5*fun2(i)


write(10,FMt="(*(es21.14,1x))") t(i),y1(i),y2(i),Te(i),k1(i),k2(i)

end do

close(10)
close(20)



stop
end program 
