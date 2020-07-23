//Dummy main program to test FEA solver
#include "FrameAnalysis.H"
int main () 
{ 
	
int 	nElements_=3;
//Create input data	
Info <<"Number of elements " <<nElements_ <<endl;
List<List<scalar>> nodes;
nodes.resize(nElements_+1);
List<List<int>>	elems;
elems.resize(nElements_);
List<List<int>> restraints;
restraints.resize(nElements_+1);
List<List<scalar>>  mats;
mats.resize(nElements_);
List<List<scalar>>  sects;
sects.resize(nElements_);
List<List<scalar>>  loads;
loads.resize(nElements_+1);
List<List<scalar>>  prescribed;
prescribed.resize(nElements_+1);	
	
FrameAnalysis FA(nodes,elems,restraints, mats,sects,loads,prescribed);
Info<< "Returned from FEA Analysis "<<FA.nodedispList()<< endl;
}
