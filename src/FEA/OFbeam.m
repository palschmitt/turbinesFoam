function beam1(nodes,elems,restraints,mats,sects,loads,prescribed)
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
    

end

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
global nnodes ndof;
totdof=nnodes*ndof;
    for idof=1:totdof
        for jdof=1:totdof
            ipos=order2(idof);
            jpos=order2(jdof);
            if (ipos <= nfree && jpos <=nfree)
                Kff(ipos,jpos)=K(idof,jdof);
            elseif (ipos > nfree && jpos <=nfree)
                Krf(ipos-nfree,jpos)=K(idof,jdof);
            elseif (ipos <= nfree && jpos > nfree)
                Kfr(ipos,jpos-nfree)=K(idof,jdof);
            else
                Krr(ipos-nfree,jpos-nfree)=K(idof,jdof);
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

function nodedisp=nodaldisp(nfree,order1)
    global nnodes ndof deff;
    totdof=ndof*nnodes;
    nodedisp=zeros(totdof,1);
    for idof=1:nfree
        nodedisp(order1(idof))=deff(idof);  
    end
end
