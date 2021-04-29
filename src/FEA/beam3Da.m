close all 


global K Kff Kfr Krf Krr;
global Ff Fr deff defr;
global nnodes nnode nelems ndof;

% read in input data
[nodes,elems,restraints,mats,sects,loads,prescribed]=readindata();

% plot nodes and elements
plotmesh(nodes,elems);

Ke=zeros(nnode*ndof); 
K=zeros(nnodes*ndof);

% Evaluate element stiffness matrix
% and assemble overall matrix
for ielem=1:nelems
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
plotdisp(nodes,elems,nodedisp);
    



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
   matvals=zeros(nmats,2);  
   textline=fgetl(fid);
   for i=1:nmats
       textline=fgetl(fid);
       [output,count]=sscanf(textline,'%f',2);
      for j=1:2
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
        x(1)=nodes(node1,1)+scale*nodedisp(ipos1+1);
        y(1)=nodes(node1,2)+scale*nodedisp(ipos1+2);
        z(1)=nodes(node1,3)+scale*nodedisp(ipos1+3);
        x(2)=nodes(node2,1)+scale*nodedisp(ipos2+1);
        y(2)=nodes(node2,2)+scale*nodedisp(ipos2+2);
        z(2)=nodes(node2,3)+scale*nodedisp(ipos2+3);
      
        plot3(x,y,z);
    end
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
   hold off;
end