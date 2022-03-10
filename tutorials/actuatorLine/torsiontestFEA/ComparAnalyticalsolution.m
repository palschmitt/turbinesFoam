clear all
close all
%Plot analytical solution against FEA AL results

%AL Data
L=1
Ls=0.9
Cl=0.8
c=0.1
v=1
nu=1E-6
rho=1
Re=v*c/nu
%FEA Data
E=11E9
Iz=7.E-6
Ix=7.E-6
A=0.003
Ip=0.05E-6
Poi=0.2

q=v^2/2*Cl*c*rho

%Torsion
F=q*L
Mt=q*L^2/2
G=E/(2*(1+Poi))
deltaphi=Mt*Ls/(Ip*G)
deltaphideg=rad2deg(deltaphi)
torquedisp=asin(deltaphi)

%Equivalent load
YForcefromLog=0.0424614199887
qequi=YForcefromLog*rho/L
Mtequi=qequi/2
deltaphiequi=Mtequi*Ls/(Ip*G)
deltaphiequideg=rad2deg(Mtequi*Ls/(Ip*G))

torquedispequi=asin(deltaphiequi)


%Read simulation results (element positions)
element0=csvread('postProcessing/actuatorBernoulliLineElements/0/leftblade.element0.csv');
element59=csvread('postProcessing/actuatorBernoulliLineElements/0/leftblade.element59.csv');

P1=element0(end,3:5);
P2=element59(end,3:5);
%Torsion angle is sin of last to first element position
deltaphisim=asin((P2(2)-P1(2))/(P2(3)-P1(3)))
deltaphisimdeg=rad2deg(deltaphisim)




figure
bar([deltaphideg, deltaphiequideg, deltaphisimdeg])
ylabel('Angle [deg]')
labels = ['Eq'; 'Eq_{Sim}'; 'Sim'];
set(gca, 'XTickLabel', labels);  
filename='Torsiontest.png'
print(filename)


figure
phis=[[deltaphideg; deltaphiequideg; deltaphisimdeg] zeros(3,1) ];
errors=[zeros(3,1) [(deltaphideg-deltaphisimdeg)/deltaphisimdeg; (deltaphiequideg-deltaphisimdeg)/deltaphisimdeg; 0 ]*100];
[AX,H1,H2] =plotyy([1:3],phis, [1:3],errors, 'bar', 'bar');
#set(H1,'FaceColor','r') % a
#set(H2,'FaceColor','b') % b
labels = [ 'Eq'; 'Eq_{Sim}'; 'Sim'];
set(AX, 'XTickLabel', labels);  
set(AX(1), 'xlim', [0 3.5]);  
set(AX(2), 'xlim', [0 3.5]);  
ylabel(AX(1), 'Angle [deg]');  
ylabel(AX(2),  'Percentage  Error [%]');  
filename='Torsiontesterrors.png'
print(filename)

figure
bar([deltaphideg, deltaphiequideg])
ylabel('Moment [Nm]')

