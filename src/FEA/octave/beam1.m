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
## @deftypefn {} {@var{retval} =} beam1 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Pal <pal@Fix>
## Created: 2020-07-13

function nodedisp=beam1(nodes,elems,restraints,mats,sects,loads,prescribed)
global K Kff Kfr Krf Krr;
global Ff Fr deff defr;
global nnodes nnode nelems ndof;

Ke=zeros(nnode*ndof); 
K=zeros(nnodes*ndof);

% Evaluate element stiffness matrix
% and assemble overall matrix
for ielem=1:nelems
    E=mats(ielem,1);
    Poi=mats(ielem,2);
    G=E/(2*(1+Poi));
    A=sects(ielem,1);
    Iz=sects(ielem,2);
    Iy=sects(ielem,3);
    J=sects(ielem,4);
    alpha=sects(ielem,5);
    elnodes=elems(ielem,:);
    dx=nodes(elnodes(2),1)-nodes(elnodes(1),1);
    dy=nodes(elnodes(2),2)-nodes(elnodes(1),2);
    dz=nodes(elnodes(2),3)-nodes(elnodes(1),3);   
    L=sqrt(dx^2+dy^2+dz^2);
    Cx=dx/L;
    Cy=dy/L;
    Cz=dz/L;
    Cxz=sqrt(Cx^2+Cz^2);
    vert=false;
    if(dx==0 && dz==0)
       vert=true;
    end
    Ke=elemstiff(E,G,A,Iz,Iy,J,L,alpha,Cx,Cy,Cz,Cxz,vert);
    assem(Ke,elnodes);
end

% Obtain the new order for DOF
[nfree,order1,order2]=neworder(restraints);

totdof=nnodes*ndof;
Kff=zeros(nfree);
Kfr=zeros(nfree,totdof-nfree);
Krf=zeros(totdof-nfree,nfree);
Krr=zeros(totdof-nfree,totdof-nfree);
Ff=zeros(nfree,1);
Fr=zeros(nfree,1);
defr=zeros(totdof-nfree,1);


% reorder sttifness matrix and partition matrix
reorder(nfree,order2);
loading(nfree,loads,order2);
predisp(nfree,prescribed,order2);


% solve for deflections
deff=Kff\(Ff-Kfr*defr);

% Find reactions
Fr=Krf*deff+Krr*defr;

% Plot deformed shape
nodedisp=nodaldisp(nfree,order1);
%plotdisp(nodes,elems,nodedisp);
    

endfunction
