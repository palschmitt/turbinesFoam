clear all
close all
%Plot analytical solution against FEA AL results

%AL Data
L=1.5
AoA=10
c=0.1
v=1
nu=1E-6
rho=1000
Re=v*c/nu
%FEA Data
E=11E9
Iz=7.0E-11
Ix=7.0E-11
A=0.003
Ip=0.05E-6

Polar= [            -90 -1.45 0.21;
                     -18 -1.45 0.21; 
                     -16 -1.3 0.165 ;
                     -14 -1.1 0.125 ;
                     -12 -0.95 0.092 ;
                     -10 -0.8 0.07 ;
                     -8 -0.64 0.05 ;
                     -6 -0.5 0.04 ;
                     -4 -0.32 0.028; 
                     -2 -0.18 0.022 ;
                     0 0.0 0.02 ;
                     2 0.18 0.022; 
                     4 0.32 0.028 ;
                     6 0.5 0.04 ;
                     8 0.64 0.05 ;
                     10 0.8 0.07 ;
                     12 0.95 0.092; 
                     14 1.1 0.125 ;
                     16 1.3 0.165 ;
                     18 1.45 0.21 ;
                     90 1.45 0.21 ];

Cl=interp1(Polar(:,1),Polar(:,2),AoA)



F=v^2/2*Cl*c
F2=0.0365630107212
x=(0:L/100:L);

%2 supports so dv=ivide results by 2
defl=F*L^3/(3*E*Iz)/2
defl2=F2*L^3/(3*E*Iz)/2

x=(0:L/100:L);
defly=F*L^3/(6*E*Iz)*(2-3*x./L+x.^3/L^3)/2;
defly2=F2*L^3/(6*E*Iz)*(2-3*x./L+x.^3/L^3)/2;



phi=rad2deg(atan(F*L^2/(2*E*Iz)))/2

AoAnew=AoA-phi
Clnew=interp1(Polar(:,1),Polar(:,2),AoAnew);
Fnew=v^2/2*Clnew*c;
deflnew=Fnew*L^3/(6*E*Iz)*(2-3*x./L+x.^3/L^3)/2;

phi2=rad2deg(atan(F2*L^2/(2*E*Iz)))/2

%Extract simulation data 
cmd="head -n 65 postProcessing/actuatorFlexibleLines/0/VTK/leftblade_000000000003.vtk | tail -n 60 > nodepos.mat"
system(cmd);
res=load("nodepos.mat");

figure
%plot(res(:,3),res(:,1))

plot(res(1:10,1),res(1:10,2),'+','linewidth',1.5)
hold on
plot(x-1.5,flip(defly),'linewidth',1.5)
plot(x-1.5,flip(defly2),'linewidth',1.5)
plot(x-1.5,flip(deflnew),'linewidth',1.5)
%plot(x,flip(defly))
legend('AL','Eq','Eq_{Sim}','Eq_{Iter}','location','northwest')
ylabel('Deflection [m]')
xlabel('[m]')
filename='AoAdeflextest.png'
print(filename)
max(res(:,2))