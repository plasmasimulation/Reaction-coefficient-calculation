clear all
close all
clc
me=9.109*10^(-31);
e=1.602*10^(-19);
kb=1.3806488*10.^(-23);
s=load('C:\Users\Bruce Xie\Desktop\matlab\Ar.txt');
Te=[];
n=1;
k1=[];
k2=[];
t=0:0.1:1000;
y1=2.336*10^(-14).*t.^(1.609).*exp(0.0618.*((log(t)).^2)-0.1171.*((log(t)).^3));
y2=2.34*10.^(-14).*t.^(0.59).*exp(-17.44./t);
while(n<=10001)
Te(n)=s(n,1);
fun1=@(v) s(n,2).*v.*exp(-me.*(v.^2)./(2.*Te(n).*11605.*kb)).*4.*pi.*v.^2;
fun2=@(v) s(n,4).*v.*exp(-me.*(v.^2)./(2.*Te(n).*11605.*kb)).*4.*pi.*v.^2;
I1=integral(fun1,0,inf);
I2=integral(fun2,0,inf);
k1(n)=(me./(2*pi*kb.*Te(n).*11605)).^1.5.*I1;
k2(n)=(me./(2*pi*kb.*Te(n).*11605)).^1.5.*I2;
n=n+1;
end
figure(1)
loglog(Te,k1,'c')
hold on
loglog(t,y1)
legend('Elastic') 
title('Reaction rates') 
xlabel('T(eV)') 
ylabel('k(m^3sec^{-1})')
figure(2)
loglog(Te,k2,'c')
hold on
loglog(t,y2)
legend('Ionization') 
title('Reaction rates') 
xlabel('T(eV)') 
ylabel('k(m^3sec^{-1})')

