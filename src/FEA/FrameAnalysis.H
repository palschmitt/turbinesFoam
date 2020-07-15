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

    // Private Member Functions


public:

    // Static Data Members

    // Constructors

        //- Construct null
        FrameAnalysis();

        //- Construct from components
        FrameAnalysis(const List<vector>& Input_);

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

        // Check

        // Edit

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
