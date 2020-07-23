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

Class
    Foam::FrameAnalysis

Description

SourceFiles
    FrameAnalysisI.H
    FrameAnalysis.C
    FrameAnalysisIO.C

\*---------------------------------------------------------------------------*/

#ifndef _Template_H
#define _Template_H

//For octave type matrix operations
#include <armadillo>

#include "List.H"
#include "dictionary.H"
#include "vector.H"
#include "unitConversion.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class FrameAnalysis;
Istream& operator>>(Istream&, FrameAnalysis&);
Ostream& operator<<(Ostream&, const FrameAnalysis&);


/*---------------------------------------------------------------------------*\
                         Class FrameAnalysis Declaration
\*---------------------------------------------------------------------------*/

class FrameAnalysis{
    // Private Data
        //- Individual blade moments about turbine origin
        List<vector> Input_;
        //- Individual blade moments about turbine origin
        List<vector> Output_;
		//Dummy to see armadillo comppiles
		arma::Mat<double> A;

		arma::Mat<double> K_, Kff_, Kfr_, Krf_, Krr_;
		arma::Mat<double> Ff_, Fr_, deff_, defr_,nodedisp_;

		double A_               ;
        double Cx_              ;
        double Cxz_             ;
        double Cy_              ;
        double Cz_              ;
        double E_               ;
        double G_               ;
        double Iy_              ;
        double Iz_              ;
        double J_               ;
        arma::Mat<double> Ke_  ;
        double L_               ;
        double Poi_             ;
        double alpha_           ;
        double dx_              ;
        double dy_              ;
        double dz_              ;
        arma::Mat<int> elems_   ;
        arma::Mat<int> elnodes_ ;
        int ielem_          ;
        arma::Mat<double> loads_  ;
        arma::Mat<double>  mats_  ;  
       int ndof_            ;
       int nelems_          ;
       int  nfree_           ;
       int nnode_           ;
       int nnodes_          ;
        arma::Mat<double> nodes_  ;      
        arma::Mat<double> order1_;       
        arma::Mat<double> order2_;       
        arma::Mat<double> prescribed_;   
        arma::Mat<int> restraints_;   
        arma::Mat<double> sects_;        
        double totdof_          ;
        bool vert_;





    // Private Member Functions


public:

    // Static Data Members

    // Constructors

        //- Construct null
        FrameAnalysis();

        //- Construct from components
        FrameAnalysis
			(
			const List<List<scalar>>& 		nodes,
			const List<List<int>>& 	elems,
			const List<List<int>>&  restraints,
			const List<List<scalar>>&  mats,
			const List<List<scalar>>&  sects,
			const List<List<scalar>>&  loads,
			const List<List<scalar>>&  prescribed
			);

        //- Construct from Istream
        FrameAnalysis(Istream&);

        //- Copy constructor
        FrameAnalysis(const FrameAnalysis&);

        //- Disallow default bitwise copy construction
        //FrameAnalysis(const FrameAnalysis&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<FrameAnalysis> New();


    //- Destructor
    ~FrameAnalysis();


    // Member Functions

        // Access 

		const arma::Mat<double>& nodedisp();
		//Return node diplacement as list
		const Foam::List<Foam::List<Foam::scalar>> nodedispList();
        // Check

        // Edit
        
        //Perofrm entire Frame Analysis
		void beam1();
        //Create Stiffness matrix for current element
        arma::Mat<double> elemstiff();
        //Assemble system Matrix
        void assem();
        //Reorder to consider restraints
        void neworder();
        
        void reorder();
        //Set loading
        void loading();
        //Set initial nodel displacement
        void predisp();
        //Evaluate nodal displacement
        void nodaldisp();
        
        
        //Create Double Matrix from nested OpenFOAM List 
        arma::Mat<double> List2Mat(const List<List<scalar>>& List);
        arma::Mat<int> List2intMat(const List<List<int>>& List);
        
        
        //And back
        const Foam::List<Foam::List<scalar>> Mat2List(const arma::Mat<double>& Mat);
        
        
        // Write


    // Member Operators

        //void operator=(const FrameAnalysis&);

        //- Disallow default bitwise assignment
        //void operator=(const FrameAnalysis&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, FrameAnalysis&);
        friend Ostream& operator<<(Ostream&, const FrameAnalysis&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "FrameAnalysisI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
