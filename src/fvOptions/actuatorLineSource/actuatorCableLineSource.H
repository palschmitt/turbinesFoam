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
    Foam::fv::actuatorCableLineSource

Description
    Actuator line source for tension-only cables, ropes, and mooring lines.

    Mirrors the structure of actuatorFlexibleLineSource but replaces
    FrameAnalysis (Euler-Bernoulli beam) with CableAnalysis (tension-only
    nonlinear cable FEA).  Key differences:

    - Elements are actuatorCableLineElement objects (3 DOF per node, no
      rotational DOF).
    - The structural model uses CableAnalysis, which performs a
      Newton-Raphson iteration for geometric nonlinearity and enforces
      tension-only behaviour (slack cable detection).
    - Each element contributes a hydrodynamic drag load in the normal-to-span
      direction; no lift is generated.
    - Pitching/torsion feedback is absent; span direction is updated purely
      from the deformed nodal geometry.
    - The recovered element tensions are written to VTK output.

    Dictionary format (elementGeometry entries per geometry point):
    \verbatim
    elementGeometry
    (
        (
            (x y z)         // [0] position
            (sx sy sz)      // [1] span direction
            (chord)         // [2] chord length (used as effective diameter)
            (crx cry crz)   // [3] chord reference direction (unused for cables,
                            //     kept for interface compatibility)
            (chordMount)    // [4] chord mount fraction (unused for cables)
            (pitch)         // [5] pitch (unused for cables)
            (EA pretension) // [6] cable material: axial stiffness, pretension
            (Cd)            // [7] cable section: drag coefficient
            (r1 r2 r3)      // [8] restraints per translational DOF (0 or 1)
        )
        ...
    );
    \endverbatim

SourceFiles
    actuatorCableLineSource.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorCableLineSource_H
#define actuatorCableLineSource_H

#include "List.H"
#include "dictionary.H"
#include "vector.H"
#include "actuatorLineSource.H"
#include "actuatorCableLineElement.H"
#include "CableAnalysis.H"
#include "volFieldsFwd.H"
#include <iomanip>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class actuatorCableLineSource Declaration
\*---------------------------------------------------------------------------*/

class actuatorCableLineSource
:
    public actuatorLineSource
{

protected:

    // Protected data

        //- Cable root position (first geometry point)
        vector rootLocation_;

        //- Cable tip / end anchor position (last geometry point)
        vector tipLocation_;
        
        List<vector> prevCableForce_;

        //- Total resultant hydrodynamic force on the cable
        vector force_;

        //- Volumetric force field contributed to the CFD momentum equation
        volVectorField forceField_;

        //- List of cable actuator elements
        PtrList<actuatorCableLineElement> elements_;

        //- Reference (undeformed) node positions for CableAnalysis.
        //  Set once in createInitialElements() and passed unchanged to every
        //  CableAnalysis call so that L0_ (unstressed element length) remains
        //  constant and is not corrupted by warm-starting from deformed nodes.
        //  Layout: nNodes = 2*nElements+1, matching the CANodes assembly in
        //  evaluateDeformation().
        List<vector> refNodePos_;


    // Protected Member Functions

        //- Create initial cable elements from elementGeometry
        void createInitialElements();

        //- Run CableAnalysis with current element loads and update positions
        void evaluateDeformation();

        //- Write VTK output file for the current time step
        void writeVTK();


public:

    //- Runtime type information
    TypeName("actuatorCableLineSource");


    // Selectors

        //- Return a reference to the selected fvOption model
        static autoPtr<actuatorCableLineSource> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Constructor
    actuatorCableLineSource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~actuatorCableLineSource();


    // Member Functions

        // Access

            //- Return reference to element pointer list
            PtrList<actuatorCableLineElement>& elements();


        // Source term addition

            //- Source term to momentum equation (incompressible)
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

#endif

// ************************************************************************* //
