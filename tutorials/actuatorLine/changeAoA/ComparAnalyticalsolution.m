clear all
close all
%Plot analytical solution against FEA AL results

%AL Data
L=1.5
Cl=0.8
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



F=v^2/2*Cl*c
F2=0.05
x=(0:L/100:L);

%2 supports so dv=ivide results by 2
defl=F*L^3/(3*E*Iz)/2
defl2=F2*L^3/(3*E*Iz)/2

phi=rad2deg(atan(F*L^2/(2*E*Iz)))/2
phi2=rad2deg(atan(F2*L^2/(2*E*Iz)))/2

%Extract simulation data 
##cmd="head -n 55 postProcessing/actuatorFlexibleLines/0/VTK/leftblade_000000000019.vtk | tail -n 50 > nodepos.mat"
##system(cmd);
##res=load("nodepos.mat");


