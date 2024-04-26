function beam1
global K Kff Kfr Krf Krr;
global M Mff Mfr Mrf Mrr;
global Ff Fr deff defr Felems;
global nnodes nnode nelems ndof;

% read in input data
[nodes,elems,restraints,mats,sects,loads,prescribed]=readindata();

% plot nodes and elements
plotmesh(nodes,elems);

Ke=zeros(nnode*ndof); 
K=zeros(nnodes*ndof);


Me=zeros(nnode*ndof);
M=zeros(nnodes*ndof);

% Evaluate element stiffness matrix
% and assemble overall matrix
for ielem=1:nelems
    E=mats(ielem,1);
    Poi=mats(ielem,2);
    Density=mats(ielem,3);
    G=E/(2*(1+Poi));
    A=sects(ielem,1);
    Iy=sects(ielem,2);
    Iz=sects(ielem,3);
    J=sects(ielem,4);
    alpha=sects(ielem,5);
    elnodes=elems(ielem,:);
    dx=nodes(elnodes(2),1)-nodes(elnodes(1),1);
    dy=nodes(elnodes(2),2)-nodes(elnodes(1),2);
    dz=nodes(elnodes(2),3)-nodes(elnodes(1),3);   
    L=sqrt(dx^2+dy^2+dz^2);
     if ielem==1
        minlen=L;
    else
        if L < minlen
            minlen = L;
        end
    end
    Cx=dx/L;
    Cy=dy/L;
    Cz=dz/L;
    Cxz=sqrt(Cx^2+Cz^2);
    vert=false;
    if(dx==0 && dz==0)
       vert=true;
    end
    [Ke,T]=elemstiff(E,G,A,Iz,Iy,J,L,alpha,Cx,Cy,Cz,Cxz,vert);
    Me=elemMass(T,Density,Poi,A,L,J);
    assem(Ke,elnodes);
    assemM(Me,elnodes);
end

% Obtain the new order for DOF
[nfree,order1,order2]=neworder(restraints);

totdof=nnodes*ndof;
Kff=zeros(nfree);
Mff=zeros(nfree);
Kfr=zeros(nfree,totdof-nfree);
Mfr=zeros(nfree,totdof-nfree);
Krf=zeros(totdof-nfree,nfree);
Mrf=zeros(totdof-nfree,nfree);
Krr=zeros(totdof-nfree,totdof-nfree);
Mrr=zeros(totdof-nfree,totdof-nfree);
Ff=zeros(nfree,1);
Fr=zeros(nfree,1);
defr=zeros(totdof-nfree,1);
Felems=zeros(ndof*nnode,nelems);

% reorder sttifness matrix and partition matrix
reorder(nfree,order2);
loading(nfree,loads,order2);
predisp(nfree,prescribed,order2);


% solve for deflections
deff=Kff\(Ff-Kfr*defr);

% Find reactions
%
% Plot deformed shape
nodedisp=nodaldisp(nfree,order1);
plotdisp(nodes,elems,nodedisp);


% Recover forces in elements
%recoverForces(sects,mats,elems,nodes,nodedisp)


% joining the beam and dynamic programs together.
k=Kff;
m=Mff;
minv=inv(m);
p=Ff;

%Obtain eigen values, only used to estimate damping matrix
%Rayleigh damping
eta=0.05;
[phi,D]=eig(k,m);
eigenvalues=zeros(nfree,1);
for i=1:nfree
    eigenvalues(i)=D(i,i);
end
%Order eigenvalues
[eigenvalues,order]=sort(eigenvalues,'ascend');
i=1;
j=3;
omegai=sqrt(eigenvalues(i));
omegaj=sqrt(eigenvalues(j));
a0=eta*((2*omegai*omegaj)/(omegai+omegaj));
a1=eta*2/(omegai+omegaj);
c=a0*m+a1*k;

% estimate time step
% estimated value is only used as guidence as yet
wavespeed=sqrt(E/(Density*(1-Poi^2)));
critdt=minlen/wavespeed;

%numerical explicit dynamics analyis.
%set time step
nt=10000;
dt=0.0000001;
%Inital calculation
u0=zeros(nfree,1);
udot0=zeros(nfree,1);
u=zeros(nfree,nt+1);
u(:,1)=u0;
udotdot0=minv*(p-c*udot0-k*u(:,1));
uminus=u(:,1)-dt*udot0+dt^2*udotdot0/2;

% Iterated over each time step
khat=m/(dt^2)+c/(2*dt);
khatinv=inv(khat);
a=m/(dt^2)-c/(2*dt);
b=k-2*m/dt^2;
time=zeros(1,nt+1);

for it=1:nt
    if it==1
        phat=p-a*uminus-b*u(:,it);
    else
        phat=p-a*u(:,it-1)-b*u(:,it);
    end
    u(:,it+1)=khatinv*phat;
    time(it+1)=time(it)+dt;
end
figure(3);
plot (time, u(22,:));
title('deflection at location u22');
xlabel('u5');
xlabel('Time sec');



    

end

function [Ke,T]=elemstiff(E,G,A,Iz,Iy,J,L,alpha,Cx,Cy,Cz,Cxz,vert)
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
%      r23=(-Cx*Cz*cosd(alpha)-Cx*sind(alpha))/Cxz;
       r23=(-Cy*Cz*cosd(alpha)+Cx*sind(alpha))/Cxz;
       r31=(Cx*Cy*sind(alpha)-Cz*cosd(alpha))/Cxz;
%      r32=Cxz*sind(alpha);
       r32=-Cxz*sind(alpha);
%       r33=(Cx*Cy*sind(alpha)+Cx*cosd(alpha))/Cxz;
       r33=(Cy*Cz*sind(alpha)+Cx*cosd(alpha))/Cxz;
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

function M=elemMass(T,Density,Poi,A,L,J)

    M=[1/3, 0, 0, 0, 0, 0, 1/6, 0, 0, 0, 0, 0; ...
       0, 13/35, 0, 0, 0 11*L/210, 0, 9/70, 0, 0, 0, -13*L/420; ...
       0, 0, 13/35, 0, -11*L/210, 0, 0, 0, 9/70, 0, 13*L/420, 0;...
       0, 0, 0, J/(3*A), 0, 0, 0, 0, 0, J/(6*A), 0, 0;...
       0, 0, -11*L/210, 0, L^2/105, 0, 0, 0, -13*L/420, 0, -L^2/140, 0; ...
       0, 11*L/210, 0, 0, 0, L^2/105, 0, 13*L/420, 0, 0, 0, -L^2/140; ...
       1/6, 0, 0, 0, 0, 0, 1/3, 0, 0, 0, 0, 0; ...
       0, 9/70, 0, 0, 0, 13*L/420, 0, 13/35, 0, 0, 0, -11*L/210; ...
       0, 0, 9/70, 0, -13*L/420, 0, 0, 0, 13/35, 0, 11*L/210, 0; ...
       0, 0, 0, J/(6*A), 0, 0, 0, 0, 0, J/(3*A), 0, 0; ...
       0, 0, 13*L/420, 0, -L^2/140, 0, 0, 0, 11*L/210, 0, L^2/105, 0; ...
       0, -13*L/420, 0, 0, 0, -L^2/140, 0, -11*L/210, 0, 0, 0, L^2/105];
    M=Density*A*L*M;
    M=T'*M*T;
end

function assem(Ke,elnodes)
    global K ndof nnode; 
    
    for innode=1:nnode
        ipos=(elnodes(innode)-1)*ndof;
        ipose=(innode-1)*ndof;
        for jnnode=1:nnode
            jpos=(elnodes(jnnode)-1)*ndof;
            jpose=(jnnode-1)*ndof;
            for idof=1:ndof
                for jdof=1:ndof
                    K(ipos+idof,jpos+jdof)=K(ipos+idof,jpos+jdof)+ ...
                        Ke(ipose+idof,jpose+jdof);
                end
            end
        end
    end
end

function assemM(Me,elnodes)
    global M ndof nnode; 
    
    for innode=1:nnode
        ipos=(elnodes(innode)-1)*ndof;
        ipose=(innode-1)*ndof;
        for jnnode=1:nnode
            jpos=(elnodes(jnnode)-1)*ndof;
            jpose=(jnnode-1)*ndof;
            for idof=1:ndof
                for jdof=1:ndof
                    M(ipos+idof,jpos+jdof)=M(ipos+idof,jpos+jdof)+ ...
                        Me(ipose+idof,jpose+jdof);
                end
            end
        end
    end
end

            
function [nfree,order1,order2]=neworder(restraints)
    global nnodes ndof;
    
    totdof=nnodes*ndof;
    order1=zeros(1,totdof);
    order2=zeros(1,totdof);
    irestraint=0;
    icount=0;
    restraintlist=zeros(1,totdof);
    for inode=1:nnodes;
        for idof=1:ndof
            icount=icount+1;
            ivalue=restraints(inode,idof);
            restraintlist(icount)=ivalue;
            if(ivalue == 1)
                irestraint=irestraint+1;
            end
        end
    end
    nrestraints=irestraint;
    nfree=totdof-nrestraints;
   
    irestraint=0;
    ifree=0;
    for idof=1:totdof;
        if(restraintlist(idof) == 1)
            irestraint=irestraint+1;
            order1(nfree+irestraint)=idof;
            order2(idof)=nfree+irestraint;
        else
            ifree=ifree+1;
            order1(ifree)=idof;
            order2(idof)=ifree;
        end
    end   
end

function reorder(nfree,order2)
global K Kff Kfr Krf Krr;
global M Mff Mfr Mrf Mrr;
global nnodes ndof;
totdof=nnodes*ndof;
    for idof=1:totdof
        for jdof=1:totdof
            ipos=order2(idof);
            jpos=order2(jdof);
            if (ipos <= nfree && jpos <=nfree)
                Kff(ipos,jpos)=K(idof,jdof);
                Mff(ipos,jpos)=M(idof,jdof);
            elseif (ipos > nfree && jpos <=nfree)
                Krf(ipos-nfree,jpos)=K(idof,jdof);
                Mrf(ipos-nfree,jpos)=M(idof,jdof);
            elseif (ipos <= nfree && jpos > nfree)
                Kfr(ipos,jpos-nfree)=K(idof,jdof);
                Mfr(ipos,jpos-nfree)=M(idof,jdof);
            else
                Krr(ipos-nfree,jpos-nfree)=K(idof,jdof);
                Mrr(ipos-nfree,jpos-nfree)=M(idof,jdof);
            end
        end
    end

end



function loading(nfree,loads,order2)
    global nnodes ndof;
    global Ff Fr;
    totdof=nnodes*ndof;
    loadvector=zeros(1,totdof);
    icount=0;
    for inode=1:nnodes;
        for idof=1:ndof
            icount=icount+1;
            ivalue=loads(inode,idof);
            loadvector(icount)=ivalue;
        end
    end
    for idof=1:totdof
        ipos=order2(idof);
        if (ipos <= nfree)
            Ff(ipos)=loadvector(idof);
        else 
            Fr(ipos-nfree)=loadvector(idof);
        end
    end    
end

function predisp(nfree,prescribed,order2)
    global nnodes ndof;
    global defr;
    totdof=nnodes*ndof;
    dispvector=zeros(1,totdof);
    icount=0;
    for inode=1:nnodes;
        for idof=1:ndof
            icount=icount+1;
            ivalue=prescribed(inode,idof);
            dispvector(icount)=ivalue;
        end
    end
    for idof=1:totdof
        ipos1=order2(idof);
        ipos2=order2(idof)-nfree;
        if (ipos1 > nfree)
            defr(ipos2)=dispvector(idof);
        end
    end    
end

function recoverForces(sects,mats,elems,nodes,nodedisp)
   global Felems ;
   global nnodes nnode nelems ndof;
   % calculate forces for each element 
   forces=zeros(ndof*2,1);
   def=zeros(ndof*2,1);
   for ielem=1:nelems
       forces=zeros(ndof*2,1);
       E=mats(ielem,1); 
       Poi=mats(ielem,2);
       G=E/(2*(1+Poi));
       A=sects(ielem,1);
       Iy=sects(ielem,2);
       Iz=sects(ielem,3);
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
       Ke=elemstiffForce(E,G,A,Iz,Iy,J,L,alpha,Cx,Cy,Cz,Cxz,vert);
       % retrive glogal deflections for element
       
       for inode=1:nnode
           for idof=1:ndof
              node=elnodes(inode);
              offset1=(node-1)*ndof;
              offset2=(inode-1)*ndof;
              def(offset2+idof)=nodedisp(offset1+idof);
           end
       end
       Felems(:,ielem)=Ke*def;
   end 
end


function Ke=elemstiffForce(E,G,A,Iz,Iy,J,L,alpha,Cx,Cy,Cz,Cxz,vert)
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
       r31=(Cx*Cy*sind(alpha)-Cz*cosd(alpha))/Cxz;
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
   
    Ke=Ke*T;
end

function [nodes,elems,restraints,mats,sects,loads,prescribed]=readindata()
   global nnodes nnode nelems ndof; 
   [filename,pathname]=uigetfile('*.inp');
   fid=fopen([pathname,filename],'r');
   textline=fgetl(fid);
   textline=fgetl(fid);
   [output,count]=sscanf(textline, '%d',3);
   nnode=output(1);
   ndof=output(2);
   ndim=output(3);
   textline=fgetl(fid);
   textline=fgetl(fid);
   [output,count]=sscanf(textline, '%d',1);
   nnodes=output;
   nodes=zeros(nnodes,ndim);
   textline=fgetl(fid);
   for i=1:nnodes
       textline=fgetl(fid);
       [output,count]=sscanf(textline,'%f',3);
      for j=1:ndim
         nodes(i,j)=output(j);
      end
   end
   textline=fgetl(fid);
   textline=fgetl(fid);
   [output,count]=sscanf(textline, '%d',1);
   nelems=output;
   elems=zeros(nelems,2);
   propno=zeros(nelems,2);
   textline=fgetl(fid);
   for i=1:nelems
       textline=fgetl(fid);
       [output,count]=sscanf(textline,'%d',4);
      for j=1:2
         elems(i,j)=output(j);
         propno(i,j)=output(j+2);
      end    
   end
   textline=fgetl(fid);
   textline=fgetl(fid);
   [output,count]=sscanf(textline, '%d',1);
   nmats=output;
   matvals=zeros(nmats,3);  
   textline=fgetl(fid);
   for i=1:nmats
       textline=fgetl(fid);
       [output,count]=sscanf(textline,'%f',3);
      for j=1:3
         matvals(i,j)=output(j);
      end    
   end
   textline=fgetl(fid);
   textline=fgetl(fid);
   [output,count]=sscanf(textline, '%d',1);
   nsects=output;
   sectvals=zeros(nsects,4);  
   textline=fgetl(fid);
   for i=1:nsects
       textline=fgetl(fid);
       [output,count]=sscanf(textline,'%f',5);
       sectvals(i,1)=output(1);
       sectvals(i,2)=output(2);
       sectvals(i,3)=output(3);
       sectvals(i,4)=output(4);
       sectvals(i,5)=output(5);
   end    
   textline=fgetl(fid);
   textline=fgetl(fid);
   [output,count]=sscanf(textline, '%d',1);
   nrestraints=output;
   restraintvals=zeros(nrestraints,13);  
   textline=fgetl(fid);
   for i=1:nrestraints
       textline=fgetl(fid);
       [output,count]=sscanf(textline,'%d %d %d %d %d %d %d %f %f %f %f %f %f',13);
      for j=1:13
         restraintvals(i,j)=output(j);
      end    
   end
   textline=fgetl(fid);
   textline=fgetl(fid);
   [output,count]=sscanf(textline, '%d',1);
   nloads=output;
   loadvals=zeros(nloads,6);
   loadnodes=zeros(nloads,1);
   textline=fgetl(fid);
   for i=1:nloads
       textline=fgetl(fid);
       [output,count]=sscanf(textline,'%d %f %f %f %f %f %f',7);
      loadnodes(i)=output(1);  
      for j=1:6
         loadvals(i,j)=output(j+1);
      end    
   end
   mats=zeros(nelems,2);
   sects=zeros(nelems,1);
   restraints=zeros(nnodes,6);
   prescribed=zeros(nnodes,6);
   loads=zeros(nnodes,6);
   for i=1:nelems
       imat=propno(i,1);
       isect=propno(i,2);
       mats(i,1)=matvals(imat,1);
       mats(i,2)=matvals(imat,2);
       mats(i,3)=matvals(imat,3);
       sects(i,1)=sectvals(isect,1);
       sects(i,2)=sectvals(isect,2);
       sects(i,3)=sectvals(isect,3);
       sects(i,4)=sectvals(isect,4);
       sects(i,5)=sectvals(isect,5);
   end
   for i=1:nrestraints
       irestraint=restraintvals(i,1);
       for j=1:ndof
          restraints(irestraint,j)=restraintvals(i,1+j);
          prescribed(irestraint,j)=restraintvals(i,7+j);
       end
   end
   for i=1:nloads
       iload=loadnodes(i);
       for j=1:ndof
          loads(iload,j)=loadvals(i,j);
       end
   end
   fclose(fid);
end

function nodedisp=nodaldisp(nfree,order1)
    global nnodes ndof deff;
    totdof=ndof*nnodes;
    nodedisp=zeros(totdof,1);
    for idof=1:nfree
        nodedisp(order1(idof))=deff(idof);  
    end
end

function plotdisp(nodes,elems,nodedisp)
    global nnodes nelems ndof;
    figure(2); 
    hold on;
    scale=20;
   
    for ielem=1:nelems
        node1=elems(ielem,1);
        node2=elems(ielem,2);
        ipos1=ndof*(node1-1);
        ipos2=ndof*(node2-1);
        x0(1)=nodes(node1,1);
        y0(1)=nodes(node1,2);
        z0(1)=nodes(node1,3);
        x0(2)=nodes(node2,1);
        y0(2)=nodes(node2,2);
        z0(2)=nodes(node2,3);
        x(1)=nodes(node1,1)+scale*nodedisp(ipos1+1);
        y(1)=nodes(node1,2)+scale*nodedisp(ipos1+2);
        z(1)=nodes(node1,3)+scale*nodedisp(ipos1+3);
        x(2)=nodes(node2,1)+scale*nodedisp(ipos2+1);
        y(2)=nodes(node2,2)+scale*nodedisp(ipos2+2);
        z(2)=nodes(node2,3)+scale*nodedisp(ipos2+3);
        plot3(x0,y0,z0,'b');
        plot3(x,y,z,'r');
    end
    title('displacements');
    hold off;
end

function plotmesh(nodes,members)
global nnodes nnode nelems ndof;
   figure(1); 
   hold on;
   for i=1:nelems
       node1=members(i,1);
       node2=members(i,2);
       linx=[nodes(node1,1),nodes(node2,1)];
       liny=[nodes(node1,2),nodes(node2,2)];
       linz=[nodes(node1,3),nodes(node2,3)];
       plot3(linx,liny,linz,'k');
   end
   for i=1:nnodes
       shift=0.05;
       xcord=nodes(i,1)+shift;
       ycord=nodes(i,2)+shift;
       zcord=nodes(i,3)+shift;
       label1=sprintf('%d',i);
       text(xcord,ycord,zcord,label1,'Color','red');
   end
   for i=1:nelems
       shift=0.05;
       node1=members(i,1);
       node2=members(i,2);
       xcord1=nodes(node1,1);
       xcord2=nodes(node2,1);
       xcord=(xcord1+xcord2)/2+shift;
       ycord1=nodes(node1,2);
       ycord2=nodes(node2,2);
       ycord=(ycord1+ycord2)/2+shift;
       zcord1=nodes(node1,2);
       zcord2=nodes(node2,2);
       zcord=(zcord1+zcord2)/2+shift;
       label1=sprintf('%d',i);
       text(xcord,ycord,zcord,label1,'Color','blue');
   end
   title('nodes and elements');
   hold off;
end

function plotSFBM(ielem)
   global Felems nelems ndof
   

   % SF or BM diagram
   choice = questdlg('Do you wnt shear force or bending moment?', ...
	'Selection Menu','Shear Force','Bending Moment','Bending Moment');
   % Handle response
   switch choice
      case 'Shear Force'
        type = 1;
      case 'Bending Moment'
        type = 2;
   end
   
   % retrieve forces
   if type == 1
     force1end1=Felems(3,ielem);
     force2end1=Felems(2,ielem);
     force1end2=Felems(9,ielem);
     force2end2=Felems(8,ielem);
   else
     force1end1=Felems(5,ielem);
     force2end1=Felems(6,ielem);
     force1end2=Felems(11,ielem);
     force2end2=Felems(12,ielem); 
   end
   figure(3)
   clf(3);
   view(3)
   daspect([1 1 1]);
   rotate3d on;
   len=10;
   scale=0.1; 
   vertices=[0,0,0; ...
             10,0,0; ...
             10,0,scale*force1end2; ...
             0,0,-scale*force1end1];  
   faces=[1,2,3,4];
   patch('Faces',faces,'Vertices',vertices,'FaceColor','yellow');
   hold on;
   offset=0.1;
   label=sprintf('%5.2f',-force1end1);
   text(offset,offset,vertices(4,3)+offset,label,'FontSize',14);
   label=sprintf('%5.2f',force1end2);
   text(10+offset,offset,vertices(3,3)+offset,label,'FontSize',14);
   
   vertices=[0,0,0; ...
             10,0,0; ...
             10,-scale*force2end2,0; ...
             0,scale*force2end1,0];
   
   if(type==1)
      vertices(3,2)=-vertices(3,2);
      vertices(4,2)=-vertices(4,2);     
   end 
   faces=[1,2,3,4];
   if(type==1)
      force2end1=-force2end1;
      force2end2=-force2end2;     
   end 
   patch('Faces',faces,'Vertices',vertices,'FaceColor','red');
   label=sprintf('%5.2f',force2end1);
   text(offset,vertices(4,2)+offset,offset,label,'FontSize',14);
   label=sprintf('%5.2f',-force2end2);
   text(10+offset,vertices(3,2)+offset,offset,label,'FontSize',14);
   
   offset=0.01;
   scale=1;
   p1 = [-1 -1 -1];                         
   p2 = [1 -1 -1];
   dp = p2-p1;
   quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),scale,'MaxHeadSize',0.5,'color','k');
   text(p2(1)+offset,p2(2)+offset,p2(3)+offset,'x','FontSize',14);
   p1 = [-1 -1 -1];                         
   p2 = [-1 2 -1];                         
   dp = p2-p1;
   quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),scale,'MaxHeadSize',0.5,'color','k');
   text(p2(1)+offset,p2(2)+offset,p2(3)+offset,'y','FontSize',14);
   p1 = [-1 -1 -1];                         
   p2 = [-1 -1 1];                         
   dp = p2-p1;
   quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),scale,'MaxHeadSize',0.5,'color','k');
   text(p2(1)+offset,p2(2)+offset,p2(3)+offset,'z','FontSize',14);
   
   if(type==1)
      label=sprintf('Shear Force Diagram \n  For element%3d',ielem); 
   else
      label=sprintf('Bending Moment Diagram \n  For element %5d',ielem);    
   end
   title(label);
   
   hold off

   
end