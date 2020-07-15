## Copyright (C) 2020 Pal
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} elemstiff (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Pal <pal@Fix>
## Created: 2020-07-13

function Ke=elemstiff(E,G,A,Iz,Iy,J,L,alpha,Cx,Cy,Cz,Cxz,vert)
    global ndof nnode;
    Ke=zeros(nnode*ndof); 
    Ke=[E*A/L,0,0,0,0,0,-E*A/L,0,0,0,0,0; ...
        0,12*E*Iz/L^3,0,0,0,6*E*Iz/L^2,0,-12*E*Iz/L^3,0,0,0,6*E*Iz/L^2; ...
        0,0,12*E*Iy/L^3,0,-6*E*Iy/L^2,0,0,0,-12*E*Iy/L^3,0,-6*E*Iy/L^2,0; ...
        0,0,0,J*G/L,0,0,0,0,0,-J*G/L,0,0; ...
        0,0,-6*E*Iy/L^2,0,4*E*Iy/L,0,0,0,6*E*Iy/L^2,0,2*E*Iy/L,0; ...
        0,6*E*Iz/L^2,0,0,0,4*E*Iz/L,0,-6*E*Iz/L^2,0,0,0,2*E*Iz/L; ...
        -E*A/L,0,0,0,0,0,E*A/L,0,0,0,0,0; ...
        0,-12*E*Iz/L^3,0,0,0,-6*E*Iz/L^2,0,12*E*Iz/L^3,0,0,0,-6*E*Iz/L^2; ...
        0,0,-12*E*Iy/L^3,0,6*E*Iy/L^2,0,0,0,12*E*Iy/L^3,0,6*E*Iy/L^2,0; ...
        0,0,0,-J*G/L,0,0,0,0,0,J*G/L,0,0; ...
        0,0,-6*E*Iy/L^2,0,2*E*Iy/L,0,0,0,6*E*Iy/L^2,0,4*E*Iy/L,0; ...
        0,6*E*Iz/L^2,0,0,0,2*E*Iz/L,0,-6*E*Iz/L^2,0,0,0,4*E*Iz/L];
    
    T=zeros(nnode*ndof);
    if ~vert
       r11=Cx;
       r12=Cy;
       r13=Cz;
       r21=(-Cx*Cy*cosd(alpha)-Cz*sind(alpha))/Cxz;
       r22=Cxz*cosd(alpha);
       r23=(-Cx*Cz*cosd(alpha)-Cx*sind(alpha))/Cxz;
       r31=(Cx*Cy*sind(alpha)-Cz*cosd(alpha))/Cxz, Cxz*sind(alpha);
       r32=Cxz*sind(alpha);
       r33=(Cx*Cy*sind(alpha)+Cx*cosd(alpha))/Cxz;
    else
        r11=0;
        r12=Cy;
        r13=0;
        r21=-Cy*cosd(alpha);
        r22=0;
        r23=sind(alpha);
        r31=Cy*sind(alpha);
        r32=0;
        r33=cosd(alpha);
    end
    T=[r11,r12,r13,0,0,0,0,0,0,0,0,0; ...
       r21,r22,r23,0,0,0,0,0,0,0,0,0; ...
       r31,r32,r33,0,0,0,0,0,0,0,0,0; ...
       0,0,0,r11,r12,r13,0,0,0,0,0,0; ...
       0,0,0,r21,r22,r23,0,0,0,0,0,0; ...
       0,0,0,r31,r32,r33,0,0,0,0,0,0; ...
       0,0,0,0,0,0,r11,r12,r13,0,0,0; ...
       0,0,0,0,0,0,r21,r22,r23,0,0,0; ...
       0,0,0,0,0,0,r31,r32,r33,0,0,0; ...
       0,0,0,0,0,0,0,0,0,r11,r12,r13; ...
       0,0,0,0,0,0,0,0,0,r21,r22,r23; ...
       0,0,0,0,0,0,0,0,0,r31,r32,r33];
   
    Ke=T'*Ke*T;
end
