/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

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
    Foam::fv::actuatorFlexibleLineSource

Description
    Actuator line class, which is a collection of actuator line elements.

SourceFiles
    actuatorFlexibleLineSource.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorFlexibleLineSource_H
#define actuatorFlexibleLineSource_H

#include "List.H"
#include "dictionary.H"
#include "vector.H"
#include "actuatorLineSource.H"
#include "actuatorBernoulliLineElement.H"
#include "cellSetOption.H"
#include "volFieldsFwd.H"
#include <iomanip>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class actuatorFlexibleLineSource Declaration
\*---------------------------------------------------------------------------*/

class actuatorFlexibleLineSource
:
public virtual actuatorLineSource
{

protected:

    // Protected data

        //- Total force vector from all elements
        vector force_;

        //- Force field from all elements
        volVectorField forceField_;

        //- List of actuator line elements
        PtrList<actuatorBernoulliLineElement> elements_;

        //- Switch for writing performance
        bool writePerf_;

        //- Switch for writing VTK
        bool writeVTK_;

        //- VTK output directory name
        fileName vtkDir_;

        //- VTK file counter
        label vtkFileSequence_;

        //- VTK file ptr
        mutable autoPtr<OFstream> vtkFilePtr_;
        
        //- Time value to track whether to move
        scalar lastMotionTime_;
        
        //- Switch for correcting end effects
        bool endEffectsActive_;

    // Protected Member Functions

        //- Create initial actuator line elements
        void createInitialElements();

        //- Create actuator line elements from updated positions
        void reevaluateElements();

        //- Read dictionary
        bool read(const dictionary& dict);
        
        //- Create the performance output file
        virtual void createOutputFile();

        //- Create the VTK output directory
        void createOutputDir();
        
        //Evaluate Deformation from forces and update positions
        void evaluateDeformation();
        
        //Does what is sais
        void writeVTK();


public:

    //- Runtime type information
    TypeName("actuatorFlexibleLineSource");

    // Selectors

    //- Return a reference to the selected fvOption model
    static autoPtr<actuatorFlexibleLineSource> New
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Constructor
    actuatorFlexibleLineSource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~actuatorFlexibleLineSource();


    // Member functions



			//Evaluate foeppl functions
			scalar foeppl(scalar x,  scalar pos, scalar exp);
			
			//Deflection of cantilever beam
			scalar Cantileverdeflection(scalar x, scalar l, scalar F, scalar EI);

        // Source term addition

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Source term to turbulence scalars
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
