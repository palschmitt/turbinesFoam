/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FrameAnalysis.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
const arma::Mat<double> Foam::FrameAnalysis::List2Mat(const List<List<scalar>>& InputList)
{
arma::Mat<double> Mat;
Mat.zeros(InputList.size(),InputList[0].size());
//Must be more fficient way, fix later
for (int k=0; k<InputList.size(); k++)
	{
		for(int l =0; l<InputList[0].size();l++)
			{
			Mat(k,l)=InputList[k][l];
			};
		};
	
	
	return Mat;
}

const Foam::List<Foam::List<Foam::scalar>> Foam::FrameAnalysis::Mat2List(const arma::Mat<double>& Mat)
{
	List<List<Foam::scalar>> OutPut;
	OutPut.resize(Mat.n_rows);
	List<Foam::scalar> SubList;
	for (uint j=0;j<Mat.n_rows;j++)
	{
		SubList.resize(Mat.n_cols);
		for (uint i=0;i<Mat.n_cols;i++)
		{
			SubList[i]=Mat(j,i);
			}
		OutPut[j]=SubList;
	}
	return OutPut;
	}


const arma::Mat<int> Foam::FrameAnalysis::List2intMat(const List<List<int>>& InputList)
{
arma::Mat<int> Mat(InputList.size(),InputList[0].size());
//Must be more fficient way, fix later
forAll(InputList, k)
	{
		forAll(InputList[k],l)
			{
			Mat(k,l)=InputList[k][l];
			};
		};
	
	
	return Mat;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FrameAnalysis::FrameAnalysis()
    {}


Foam::FrameAnalysis::FrameAnalysis
(
const List<List<scalar>>& 	nodes,
const List<List<int>>& 		elems,
const List<List<int>>& 		restraints,
const List<List<scalar>>&  	mats,
const List<List<scalar>>&  	sects,
const List<List<scalar>>&  	loads,
const List<List<scalar>>&  	prescribed
)
{
//Create and set correspinding armadillo Matrices

nelems_=elems.size();
nnodes_=nodes.size();
nnode_=2; //Not variable at the moment, beam element has 6DoF
ndof_=6;
nodes_=List2Mat(nodes);
//std::cout << "nodes_: "<< nodes_ <<endl;
elems_=List2intMat(elems);
//std::cout << "elems_: "<< elems_ <<endl;
restraints_=List2intMat(restraints);
//std::cout << "restraints_: "<< restraints_ <<endl;
sects_=List2Mat(sects);
//std::cout << "sects_: "<< sects_ <<endl;
loads_=List2Mat(loads);
//std::cout << "loads_: "<< loads_ <<endl;
mats_=List2Mat(mats);
//std::cout << "mats_: "<< mats_ <<endl;
prescribed_=List2Mat(prescribed);
//std::cout << "prescribed_: "<< prescribed_ <<endl;

beam1();//Evaluate deformation
}


Foam::FrameAnalysis::FrameAnalysis(const FrameAnalysis&)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::FrameAnalysis>
Foam::FrameAnalysis::New()
{
    return autoPtr<FrameAnalysis>(new FrameAnalysis);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FrameAnalysis::~FrameAnalysis()
{
    Info<<"Called FrameAnalysis Destructor."<<endl;
    }


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const arma::Mat<double>&  Foam::FrameAnalysis::nodedisp()
{
	return nodedisp_;
	}
	
const Foam::List<Foam::List<Foam::scalar>> Foam::FrameAnalysis::nodedispList()
{
	return Mat2List(nodedisp_);
	}
const arma::Mat<double>&  Foam::FrameAnalysis::Felems()
{
	return Felems_;
	}
	
const Foam::List<Foam::List<Foam::scalar>> Foam::FrameAnalysis::FelemsList()
{
	return Mat2List(Felems_);
	}
	
void Foam::FrameAnalysis::beam1()
{
    Ke_.zeros(nnode_*ndof_,nnode_*ndof_); 
    K_.zeros(nnodes_*ndof_,nnodes_*ndof_);

// Evaluate element stiffness matrix
// and assemble overall matrix
  for(ielem_=0; ielem_ < nelems_;ielem_++)
  {
    E_=mats_(ielem_,0);
    Poi_=mats_(ielem_,1);
    G_=E_/(2*(1+Poi_));
    A_=sects_(ielem_,0);
    Iy_=sects_(ielem_,1);
    Iz_=sects_(ielem_,2);
    J_=sects_(ielem_,3);
    alpha_=sects_(ielem_,4);
    //Check for nicer way of doing this
    elnodes_.zeros(2,1) ;
    elnodes_(0)=elems_(ielem_,0);
    elnodes_(1)=elems_(ielem_,1);
    dx_=nodes_(elnodes_(1),0)-nodes_(elnodes_(0),0);
    dy_=nodes_(elnodes_(1),1)-nodes_(elnodes_(0),1);
    dz_=nodes_(elnodes_(1),2)-nodes_(elnodes_(0),2);   
    arma::vec Lv = {dx_, dy_, dz_};
    L_=arma::norm(Lv);
    
    if (fabs(dx_)<SMALL)
    Cx_=0.;
		else
    Cx_=dx_/L_;
    
     if (fabs(dy_)<SMALL)
    Cy_=0.;
    		else
    Cy_=dy_/L_;
    
     if (fabs(dz_)<SMALL)
    Cz_=0.;
    		else
    Cz_=dz_/L_;
    
    arma::vec L2={Cx_, Cz_};
    
    Cxz_=arma::norm(L2);
    vert_=false;
    if((fabs(dx_)<SMALL) && (fabs(dz_)<SMALL))
    {
       vert_=true;
    }
    elemstiff();
    //std::cout << "Ke_: "<< Ke_ <<endl;
    assem();
    //std::cout << "K_: "<< K_ <<endl;
}


// Obtain the new order for DOF
neworder();

totdof_=nnodes_*ndof_;
Kff_.zeros(nfree_,nfree_);
Kfr_.zeros(nfree_,totdof_-nfree_);
Krf_.zeros(totdof_-nfree_,nfree_);
Krr_.zeros(totdof_-nfree_,totdof_-nfree_);
Ff_.zeros(nfree_,1);
Fr_.zeros(totdof_-nfree_,1);
defr_.zeros(totdof_-nfree_,1);


// reorder sttifness matrix and partition matrix
reorder();
loading();
predisp();

//std::cout << "Kff_: "<<endl<< Kff_ <<endl;
//std::cout << "Kfr_: "<<endl<< Kfr_ <<endl;
//std::cout << "Ff: "<<endl<< Ff_ <<endl;
//std::cout << "defr: "<<endl<< defr_ <<endl;
//std::cout << "deff: "<<endl<< deff_ <<endl;

// solve for deflections
deff_=solve(Kff_,(Ff_-Kfr_*defr_),arma::solve_opts::no_approx);

// Find reactions
Fr_=Krf_*deff_+Krr_*defr_;


nodaldisp();
recoverforces();
	
	}

void Foam::FrameAnalysis::elemstiff()
{
//Ke=zeros(nnode*ndof); 
// A = { {1, 3, 5},
   //       {2, 4, 6} };
scalar L2=pow(L_,2);
scalar L3=pow(L_,3);

//Make E,A, L,,Iz...,Cx.. function paramaters and remove from global	
    Ke_={{E_*A_/L_,0,0,0,0,0,-E_*A_/L_,0,0,0,0,0},
        {0,12*E_*Iz_/L3,0,0,0,6*E_*Iz_/L2,0,-12*E_*Iz_/L3,0,0,0,6*E_*Iz_/L2},
        {0,0,12*E_*Iy_/L3,0,-6*E_*Iy_/L2,0,0,0,-12*E_*Iy_/L3,0,-6*E_*Iy_/L2,0},
        {0,0,0,J_*G_/L_,0,0,0,0,0,-J_*G_/L_,0,0},
        {0,0,-6*E_*Iy_/L2,0,4*E_*Iy_/L_,0,0,0,6*E_*Iy_/L2,0,2*E_*Iy_/L_,0},
        {0,6*E_*Iz_/L2,0,0,0,4*E_*Iz_/L_,0,-6*E_*Iz_/L2,0,0,0,2*E_*Iz_/L_},
        {-E_*A_/L_,0,0,0,0,0,E_*A_/L_,0,0,0,0,0},
        {0,-12*E_*Iz_/L3,0,0,0,-6*E_*Iz_/L2,0,12*E_*Iz_/L3,0,0,0,-6*E_*Iz_/L2},
        {0,0,-12*E_*Iy_/L3,0,6*E_*Iy_/L2,0,0,0,12*E_*Iy_/L3,0,6*E_*Iy_/L2,0},
        {0,0,0,-J_*G_/L_,0,0,0,0,0,J_*G_/L_,0,0},
        {0,0,-6*E_*Iy_/L2,0,2*E_*Iy_/L_,0,0,0,6*E_*Iy_/L2,0,4*E_*Iy_/L_,0},
        {0,6*E_*Iz_/L2,0,0,0,2*E_*Iz_/L_,0,-6*E_*Iz_/L2,0,0,0,4*E_*Iz_/L_}};
    
    
    
    arma::Mat<double> RotMat;
    RotMat.zeros(nnode_*ndof_, nnode_*ndof_);
    alpha_=degToRad(alpha_);
    
    scalar r11;
    scalar r12;
    scalar r13;
    scalar r21;
    scalar r22;
    scalar r23;
    scalar r31;
    scalar r32;
    scalar r33;
    
    
    if (!vert_)
    {
       r11=Cx_;
       r12=Cy_;
       r13=Cz_;
       r21=(-dot(dot(Cx_,Cy_),cos(alpha_))-dot(Cz_,sin(alpha_)))/Cxz_;
       r22=dot(Cxz_,cos(alpha_));
       r23=(-dot(dot(Cy_,Cz_),cos(alpha_))+dot(Cx_,sin(alpha_)))/Cxz_;
       r31=(dot(dot(Cx_,Cy_),sin(alpha_))-dot(Cz_,cos(alpha_)))/Cxz_;
       r32=dot(-Cxz_,sin(alpha_));
       r33=(Cy_*Cz_*sin(alpha_)+Cx_*cos(alpha_))/Cxz_;
   }
    else
    {
        r11=0;
        r12=Cy_;
        r13=0;
        r21=-Cy_*sin(alpha_);
        r22=0;
        r23=cos(alpha_);
        r31=Cy_*cos(alpha_);
        r32=0;
        r33=cos(alpha_);
	}
    
    RotMat={{r11,r12,r13,0,0,0,0,0,0,0,0,0},
       {r21,r22,r23,0,0,0,0,0,0,0,0,0},
       {r31,r32,r33,0,0,0,0,0,0,0,0,0},
       {0,0,0,r11,r12,r13,0,0,0,0,0,0},
       {0,0,0,r21,r22,r23,0,0,0,0,0,0},
       {0,0,0,r31,r32,r33,0,0,0,0,0,0},
       {0,0,0,0,0,0,r11,r12,r13,0,0,0},
       {0,0,0,0,0,0,r21,r22,r23,0,0,0},
       {0,0,0,0,0,0,r31,r32,r33,0,0,0},
       {0,0,0,0,0,0,0,0,0,r11,r12,r13},
       {0,0,0,0,0,0,0,0,0,r21,r22,r23}, 
       {0,0,0,0,0,0,0,0,0,r31,r32,r33}};
   
    Ke_=RotMat.t()*Ke_*RotMat;
    //std::cout <<" Ke_ "<<Ke_ <<endl;
	}

void Foam::FrameAnalysis::elemstiffforce()
{
	//Ke=zeros(nnode*ndof); 
	
	// A = { {1, 3, 5},
   //       {2, 4, 6} };
scalar L2=pow(L_,2);
scalar L3=pow(L_,3);

//Make E,A, L,,Iz...,Cx.. function paramaters and remove from global	
    Ke_={{E_*A_/L_,0,0,0,0,0,-E_*A_/L_,0,0,0,0,0},
        {0,12*E_*Iz_/L3,0,0,0,6*E_*Iz_/L2,0,-12*E_*Iz_/L3,0,0,0,6*E_*Iz_/L2},
        {0,0,12*E_*Iy_/L3,0,-6*E_*Iy_/L2,0,0,0,-12*E_*Iy_/L3,0,-6*E_*Iy_/L2,0},
        {0,0,0,J_*G_/L_,0,0,0,0,0,-J_*G_/L_,0,0},
        {0,0,-6*E_*Iy_/L2,0,4*E_*Iy_/L_,0,0,0,6*E_*Iy_/L2,0,2*E_*Iy_/L_,0},
        {0,6*E_*Iz_/L2,0,0,0,4*E_*Iz_/L_,0,-6*E_*Iz_/L2,0,0,0,2*E_*Iz_/L_},
        {-E_*A_/L_,0,0,0,0,0,E_*A_/L_,0,0,0,0,0},
        {0,-12*E_*Iz_/L3,0,0,0,-6*E_*Iz_/L2,0,12*E_*Iz_/L3,0,0,0,-6*E_*Iz_/L2},
        {0,0,-12*E_*Iy_/L3,0,6*E_*Iy_/L2,0,0,0,12*E_*Iy_/L3,0,6*E_*Iy_/L2,0},
        {0,0,0,-J_*G_/L_,0,0,0,0,0,J_*G_/L_,0,0},
        {0,0,-6*E_*Iy_/L2,0,2*E_*Iy_/L_,0,0,0,6*E_*Iy_/L2,0,4*E_*Iy_/L_,0},
        {0,6*E_*Iz_/L2,0,0,0,2*E_*Iz_/L_,0,-6*E_*Iz_/L2,0,0,0,4*E_*Iz_/L_}};
    
    
    
    arma::Mat<double> RotMat;
    RotMat.zeros(nnode_*ndof_, nnode_*ndof_);
    alpha_=degToRad(alpha_);
    
    scalar r11;
    scalar r12;
    scalar r13;
    scalar r21;
    scalar r22;
    scalar r23;
    scalar r31;
    scalar r32;
    scalar r33;
    
    
    if (!vert_)
    {
       r11=Cx_;
       r12=Cy_;
       r13=Cz_;
       r21=(-dot(dot(Cx_,Cy_),cos(alpha_))-dot(Cz_,sin(alpha_)))/Cxz_;
       r22=dot(Cxz_,cos(alpha_));
       r23=(-dot(dot(Cy_,Cz_),cos(alpha_))+dot(Cx_,sin(alpha_)))/Cxz_;
       r31=(dot(dot(Cx_,Cy_),sin(alpha_))-dot(Cz_,cos(alpha_)))/Cxz_;
       r32=dot(-Cxz_,sin(alpha_));
       r33=(Cy_*Cz_*sin(alpha_)+Cx_*cos(alpha_))/Cxz_;
   }
    else
    {
        r11=0;
        r12=Cy_;
        r13=0;
        r21=-Cy_*sin(alpha_);
        r22=0;
        r23=cos(alpha_);
        r31=Cy_*cos(alpha_);
        r32=0;
        r33=cos(alpha_);
	}
    
    RotMat={{r11,r12,r13,0,0,0,0,0,0,0,0,0},
       {r21,r22,r23,0,0,0,0,0,0,0,0,0},
       {r31,r32,r33,0,0,0,0,0,0,0,0,0},
       {0,0,0,r11,r12,r13,0,0,0,0,0,0},
       {0,0,0,r21,r22,r23,0,0,0,0,0,0},
       {0,0,0,r31,r32,r33,0,0,0,0,0,0},
       {0,0,0,0,0,0,r11,r12,r13,0,0,0},
       {0,0,0,0,0,0,r21,r22,r23,0,0,0},
       {0,0,0,0,0,0,r31,r32,r33,0,0,0},
       {0,0,0,0,0,0,0,0,0,r11,r12,r13},
       {0,0,0,0,0,0,0,0,0,r21,r22,r23}, 
       {0,0,0,0,0,0,0,0,0,r31,r32,r33}};
   
    Ke_=Ke_*RotMat;
    //std::cout <<" Ke_ "<<Ke_ <<endl;
	}
void Foam::FrameAnalysis::assem()
{
	scalar ipos,ipose,jpos,jpose;

	 for (int innode=0;innode< nnode_; innode++)
	 {
        ipos=(elnodes_(innode))*ndof_;
        ipose=(innode)*ndof_;
        for (int jnnode=0;jnnode<nnode_;jnnode++)
        {
            jpos=(elnodes_(jnnode))*ndof_;
            jpose=(jnnode)*ndof_;
            for (int idof=0; idof<ndof_;idof++ )
            {
                for (int jdof=0; jdof<ndof_; jdof++)
                {
                    K_(ipos+idof,jpos+jdof)=K_(ipos+idof,jpos+jdof)+Ke_(ipose+idof,jpose+jdof);
				}
			}
		}
	}
}

void Foam::FrameAnalysis::neworder()
{
		//ielem_=0; ielem_ < nelems_;ielem_++
	int 	ivalue;
	int nrestraints, ifree;
	totdof_=nnodes_*ndof_;
    order1_.zeros(1,totdof_);
    order2_.zeros(1,totdof_);
    int irestraint=0;
    int icount=0;
    arma::Mat<int> restraintlist;
    restraintlist.zeros(1,totdof_);
    for (int inode=0; inode<nnodes_;inode++)
    {
        for (int idof=0; idof<ndof_; idof++)
        {
            ivalue=restraints_(inode,idof);
            restraintlist(icount)=ivalue;
            icount=icount+1;
            if (ivalue== int(1))
            {
                irestraint=irestraint+1;
            }
        }
    }
    nrestraints=irestraint;
    nfree_=totdof_-nrestraints;
   
    irestraint=0;
    ifree=0;
    for (int idof=0;idof<totdof_; idof++)
    {
        if(restraintlist(idof) == 1)
        {
            order1_(nfree_+irestraint)=idof;
            order2_(idof)=nfree_+irestraint;
            irestraint=irestraint+1;
		}
        else
        {
            order1_(ifree)=idof;
            order2_(idof)=ifree;
            ifree=ifree+1;
        }
    }   
	//std::cout << "order1:\n" << order1_ << "\n";
	//std::cout << "order2:\n" << order2_ << "\n";
	}

void Foam::FrameAnalysis::reorder()
{
	int ipos,jpos;
	totdof_=nnodes_*ndof_;
    for (uint idof=0;idof<totdof_;idof++)
    {
        for (uint jdof=0;jdof<totdof_;jdof++)
        {
            ipos=order2_(idof);
            jpos=order2_(jdof);
            if (ipos+1 <= nfree_ && jpos+1 <=nfree_)
                Kff_(ipos,jpos)=K_(idof,jdof);
            else if (ipos+1 > nfree_ && jpos+1 <=nfree_)
                Krf_(ipos-nfree_,jpos)=K_(idof,jdof);
            else if (ipos+1 <= nfree_ && jpos+1 > nfree_)
                Kfr_(ipos,jpos-nfree_)=K_(idof,jdof);
            else
                Krr_(ipos-nfree_,jpos-nfree_)=K_(idof,jdof);
            
        }
    }

}

void  Foam::FrameAnalysis::loading()
{
	scalar ivalue;
	int ipos;
	totdof_=nnodes_*ndof_;
    arma::Mat<double>  loadvector;
    loadvector.zeros(1,totdof_);
    int icount=0;
    for (int inode=0; inode<nnodes_; inode++)
    {
        for (int idof=0;idof<ndof_; idof++)
        {
            ivalue=loads_(inode,idof);
            loadvector(icount)=ivalue;
            icount=icount+1;
        }
    }
    for (int idof=0; idof<totdof_;idof++)
    {
        ipos=order2_(idof);
        if (ipos+1 <= nfree_)
            Ff_(ipos)=loadvector(idof);
        else 
            Fr_(ipos-nfree_)=loadvector(idof);
        
    }    
	}


void Foam::FrameAnalysis::predisp()
{
	totdof_=nnodes_*ndof_;
	scalar ivalue;
	int ipos1,ipos2;
    arma::Mat<double>  dispvector;
    dispvector.zeros(1,totdof_);
    int icount=0;
    for (int inode=0;inode<nnodes_;inode++)
    {
        for (int idof=0;idof<ndof_;idof++)
        {
            ivalue=prescribed_(inode,idof);
            dispvector(icount)=ivalue;
            icount=icount+1;
        }
    }
    for (int idof=0; idof<totdof_;idof++)
    {
        ipos1=order2_(idof);
        ipos2=order2_(idof)-nfree_;
        if (ipos1+1 > nfree_)
            defr_(ipos2)=dispvector(idof);
        
    }    
	}

void Foam::FrameAnalysis::nodaldisp()
{

	double totdof_=ndof_*nnodes_;
    nodedisp_.zeros(totdof_,1);
    for (int idof=0;idof<nfree_;idof++)
    {
        nodedisp_(order1_(idof))=deff_(idof);  
	}
    //Reformat to get xyzrxryrz per node in list
    nodedisp_.reshape(6,nnodes_);
    inplace_trans(nodedisp_);
    //std::cout << "nodedisp_ in FEA: "<< nodedisp_ <<endl;

	}

void Foam::FrameAnalysis::recoverforces()
{
    //TODO: initilise Elementforce matrix, 
    //Remove double evaluationg of node displacement?
//Felems_
arma::Mat<double>  def(ndof_*2,1);
Felems_.set_size(nelems_,ndof_*2);
  for(ielem_=0; ielem_ < nelems_;ielem_++)
  {
    E_=mats_(ielem_,0);
    Poi_=mats_(ielem_,1);
    G_=E_/(2*(1+Poi_));
    A_=sects_(ielem_,0);
    Iy_=sects_(ielem_,1);
    Iz_=sects_(ielem_,2);
    J_=sects_(ielem_,3);
    alpha_=sects_(ielem_,4);
    //Check for nicer way of doing this
    elnodes_.zeros(2,1) ;
    elnodes_(0)=elems_(ielem_,0);
    elnodes_(1)=elems_(ielem_,1);
    dx_=nodes_(elnodes_(1),0)-nodes_(elnodes_(0),0);
    dy_=nodes_(elnodes_(1),1)-nodes_(elnodes_(0),1);
    dz_=nodes_(elnodes_(1),2)-nodes_(elnodes_(0),2);   
    arma::vec Lv = {dx_, dy_, dz_};
    L_=arma::norm(Lv);
    
    if (fabs(dx_)<SMALL)
    Cx_=0.;
		else
    Cx_=dx_/L_;
    
     if (fabs(dy_)<SMALL)
    Cy_=0.;
    		else
    Cy_=dy_/L_;
    
     if (fabs(dz_)<SMALL)
    Cz_=0.;
    		else
    Cz_=dz_/L_;
    
    arma::vec L2={Cx_, Cz_};
    
    Cxz_=arma::norm(L2);
    vert_=false;
    if((fabs(dx_)<SMALL) && (fabs(dz_)<SMALL))
    {
       vert_=true;
    }
    elemstiffforce();
     
     //retrive global deflections for element
    for (int inode=0;inode<nnode_;inode++)
    {
	for (int idof=0;idof<ndof_;idof++)
	{
	    int node=elnodes_(inode);
	    int offset1=(node)*ndof_;
	    int offset2=(inode)*ndof_;
	    def(offset2+idof)=nodedisp_(offset1+idof);
	}
	}
	Felems_.row(ielem_)=trans(Ke_*def);
    }

}
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//void Foam::FrameAnalysis::operator=(const FrameAnalysis& rhs)
//{
    //// Check for assignment to self
    //if (this == &rhs)
    //{
        //FatalErrorInFunction
            //<< "Attempted assignment to self"
            //<< abort(FatalError);
    //}
//}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
