clear all
close all
%Plot analytical solution against FEA AL results

%AL Data
L=1
Cl=0.8
c=0.1
v=1
nu=1E-6
rho=1000
Re=v*c/nu
%FEA Data
E=11E9
Poi=0.2
Iz=7.E-6
Ix=7.E-6
A=0.003
Ip=0.05E-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bending

q=v^2/2*Cl*c*rho
x=(0:L/100:L);
defly=q*L^4/(24*E*Iz)*(x.^4/L^4-4*x/L+3);
deflz=q*L^4/(24*E*Ix)*(x.^4/L^4-4*x/L+3);

%Extract simulation data 
cmd="head -n 55 postProcessing/actuatorFlexibleLines/0/VTK/leftblade_000000000019.vtk | tail -n 50 > nodepos.mat"
system(cmd);
res=load("nodepos.mat");


figure
%plot(res(:,3),res(:,1))

plot(res(:,3),res(:,2))
hold on
plot(x,flip(deflz))
%plot(x,flip(defly))
legend('AL','Analytical')
ylabel('Deflection [m]')
xlabel('[m]')
filename='Bendingtest.png'
print(filename)
max(res(:,2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Torsion
Mt=q/2
G=E/(2*(1+Poi))
deltaphi=Mt*1/(Ip*G)
deltaphideg=rad2deg(deltaphi)
torquedisp=asin(deltaphi)