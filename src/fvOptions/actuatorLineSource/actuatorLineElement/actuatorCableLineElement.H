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
    Foam::fv::actuatorCableLineElement

Description
    Actuator line element for tension-only cables/ropes/mooring lines.

    Derived from actuatorLineElement in the same way as
    actuatorBernoulliLineElement, but adapted for CableAnalysis:

    - No bending/torsion: only 3 translational DOF per node (x, y, z).
    - Structural data reduced to [E, A] (material) and pretension T0.
    - No rotational deformation is fed back to chord direction; instead the
      element span direction is updated from the deformed nodal positions.
    - The applied (hydrodynamic) drag force is computed from the local
      relative velocity projected onto the plane normal to the span axis,
      using a user-supplied drag coefficient Cd.  No lift is generated.
    - The structural tension recovered from CableAnalysis is stored and
      written to VTK for post-processing.
    - Buoyancy and self-weight are computed from the cable material density
      and cross-sectional area, using the local fluid density sampled from
      the 'rho' field (or a constant reference value for incompressible
      solvers) and the mesh gravity vector 'g'.  The net body force
      (buoyancy minus self-weight) is added to cableForce_ so that it is
      passed to CableAnalysis alongside the hydrodynamic drag.

    Input dictionary entries (in addition to those required by
    actuatorLineElement):
    \verbatim
        P1                  (x y z);  // element start node position
        P2                  (x y z);  // element end node position
        CableEA             value;    // axial stiffness  E*A  [N]
        CablePretension     value;    // initial pretension T0  [N]
        CableDragCoeff      value;    // normal drag coefficient Cd [-]
        CableRestraints     (r1 r2 r3); // 1=fixed, 0=free, per DOF
        CableDiameter       value;    // outer diameter d  [m]
        CableDensity        value;    // cable material density rho_c [kg/m3]
        CableFluidDensity   value;    // reference fluid density rho_f [kg/m3]
                                      //   (used for incompressible solvers;
                                      //    ignored if a 'rho' volScalarField
                                      //    is present in the mesh)
        CableGravity        (gx gy gz); // gravity vector [m/s2]
                                        //   default: read from mesh 'g' field
    \endverbatim

SourceFiles
    actuatorCableLineElement.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorCableLineElement_H
#define actuatorCableLineElement_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "actuatorLineElement.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class actuatorCableLineElement Declaration
\*---------------------------------------------------------------------------*/

class actuatorCableLineElement
:
    public actuatorLineElement
{

protected:

    // Protected data

        //- Element start point (node 1 of the CableAnalysis element)
        vector P1_;

        //- Element end point (node 2 of the CableAnalysis element)
        vector P2_;

        //- Accumulated total deformation applied to this element [m]
        vector deformation_;

        //- Nodal force applied to the structural model (fluid drag force) [N]
        vector cableForce_;

        //- Element tension recovered from CableAnalysis [N]
        scalar tension_;

        // --- CableAnalysis structural parameters ---

        //- Axial stiffness E*A  [N]
        scalar cableEA_;

        //- Initial pretension T0 [N]
        scalar cablePretension_;

        //- Normal (cross-flow) drag coefficient Cd [-]
        scalar cableDragCoeff_;

        //- Restraint flags per translational DOF (3 entries: x, y, z)
        //  1 = fixed, 0 = free
        List<int> cableRestraints_;

        // --- Buoyancy / self-weight parameters ---

        //- Outer diameter of the cable cross-section [m]
        //  Used to compute displaced volume = pi/4 * d^2 * spanLength
        scalar cableDiameter_;

        //- Cable material density [kg/m3]
        scalar cableDensity_;

        //- Reference fluid density [kg/m3]
        //  Used for incompressible solvers that have no 'rho' field.
        //  If a 'rho' volScalarField exists in the mesh it is sampled
        //  at the element position instead.
        scalar cableFluidDensity_;

        //- Gravity vector [m/s2] (positive direction = downward by convention)
        //  Read from the mesh 'g' uniformDimensionedVectorField if available,
        //  otherwise taken from dictionary entry CableGravity.
        vector gravity_;

        //- Net buoyancy force on element (buoyancy - self-weight) [N]
        //  Recomputed each call to calculateForce(); stored for VTK output.
        vector buoyancyForce_;

        //- Reference (unstressed) span length [m].
        //  Set once at construction from the initial P1–P2 distance and used
        //  for buoyancy volume and projected drag area in calculateForce().
        //  Using the deformed span length causes a divergent feedback loop:
        //  buoyancy grows with stretch, driving further displacement.
        scalar refSpanLength_;

        //- Reference (unstressed) span direction (unit vector along P2-P1).
        //  Set once after construction from the initial element geometry and
        //  used in calculateForce() to project the inflow velocity onto the
        //  cable axis.  Using the deformed spanDirection_ after a buoyancy-
        //  driven catenary deformation introduces a phantom drag component:
        //  the span rotates away from the flow direction, relNormal picks up
        //  the full inflow velocity, and drag diverges exponentially.
        vector refSpanDirection_;

        //- True when this element's position_ was found in the CFD mesh on
        //  the most recent calculateForce() call.  Set before applyForceField
        //  is called so that applyForceField can skip the base-class
        //  calcProjectionEpsilon() FatalError for out-of-domain nodes.
        bool positionInMesh_;


    // Protected Member Functions

        //- (none beyond base class)


public:

    //- Runtime type information
    TypeName("actuatorCableLineElement");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        actuatorCableLineElement,
        dictionary,
        (
            const dictionary& dict,
            const word& modelName
        ),
        (dict, modelName)
    );


    // Selectors

        //- Return a reference
        static autoPtr<actuatorCableLineElement> New(const dictionary& dict);


    //- Constructor
    actuatorCableLineElement
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~actuatorCableLineElement();


    // Member Functions

        // Access

            //- Return element start node position
            const vector& P1() const;

            //- Return element end node position
            const vector& P2() const;

            //- Return accumulated deformation
            const vector& deformation() const;

            //- Return the structural (drag) force vector
            const vector& cableForce() const;

            //- Return element tension from last CableAnalysis solve
            scalar tension() const;

            //- Return axial stiffness EA
            scalar cableEA() const;

            //- Return pretension T0
            scalar cablePretension() const;

            //- Return drag coefficient
            scalar cableDragCoeff() const;

            //- Return restraint flags (3 entries)
            const List<int>& cableRestraints() const;

            //- Return reference (unstressed) span length
            scalar refSpanLength() const;

            //- Return reference (unstressed) span direction
            const vector& refSpanDirection() const;

            //- Return cable outer diameter
            scalar cableDiameter() const;

            //- Return cable material density
            scalar cableDensity() const;

            //- Return reference fluid density
            scalar cableFluidDensity() const;

            //- Return gravity vector
            const vector& gravity() const;

            //- Return last computed net buoyancy force
            const vector& buoyancyForce() const;

            //- Return span direction (unit vector along deformed element axis)
            const vector& spanDirection() const;


        // Manipulation

            //- Set element start position
            void setP1(vector P1);

            //- Set element end position
            void setP2(vector P2);

            //- Set accumulated deformation (additive)
            void setDeformation(vector deformation);

            //- Set the element position (midpoint)
            void setPosition(vector newPosition);

            //- Set span direction (unit vector along deformed element axis)
            void setSpanDirection(vector spanDir);

            //- Set span length
            void setSpanLength(scalar spanLength);

            //- Set reference (unstressed) span length — call once after construction
            void setRefSpanLength(scalar refSpanLength);

            //- Set reference (unstressed) span direction — call once after construction
            void setRefSpanDirection(vector refSpanDir);

            //- Set element tension
            void setTension(scalar tension);

            //- Rotate element (positions, span direction, velocity)
            void rotate
            (
                vector rotationPoint,
                vector axis,
                scalar radians,
                bool rotateVelocity = true
            );

            void rotate
            (
                vector rotationPoint,
                vector axis,
                scalar radians
            );


        // Evaluation

            //- Read element data from dictionary
            void read();

            //- Calculate the total structural force on the element:
            //    cableForce_ = hydrodynamic drag (normal component)
            //                + net buoyancy (buoyancy - self-weight)
            //  The result is stored in forceVector_ (base class) and
            //  cableForce_ for CableAnalysis assembly.
            void calculateForce(const volVectorField& Uin);

            //- Source term to momentum equation.
            //  Overridden to skip applyForceField when position is outside
            //  the CFD domain (which is valid for cable anchor nodes).
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                volVectorField& forceField
            );

            //- Source term to compressible momentum equation.
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                volVectorField& forceField
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
