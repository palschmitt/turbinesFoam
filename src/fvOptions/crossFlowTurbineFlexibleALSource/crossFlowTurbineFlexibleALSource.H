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
    Foam::fv::crossFlowTurbineFlexibleALSource

Description
    Cross-flow turbine with flexible (Bernoulli beam) actuator line blades.

    Combines crossFlowTurbineALSource (turbine geometry, rotation, performance)
    with actuatorFlexibleLineSource (structural deformation via FrameAnalysis).

    The diamond-inheritance problem that arises from both base classes
    independently deriving from cellSetOption / option is resolved by:

      1. Declaring crossFlowTurbineALSource and actuatorFlexibleLineSource
         both with 'virtual' in their inheritance of cellSetOption (requires
         editing those headers — see note below).

      2. Providing explicit final overrides here for every virtual function
         that would otherwise have two equally-valid definitions.

    NOTE: If you cannot modify the base-class headers to add 'virtual', the
    alternative is single inheritance — inherit only from
    crossFlowTurbineALSource and hold an actuatorFlexibleLineSource by
    composition or replace createBlades() to produce flexible blades directly
    (see crossFlowTurbineFlexibleALSource::createBlades() in the .C file).
    The implementation below uses the compositional approach so that no changes
    to the existing headers are required.

SourceFiles
    crossFlowTurbineFlexibleALSource.C

\*---------------------------------------------------------------------------*/

#ifndef crossFlowTurbineFlexibleALSource_H
#define crossFlowTurbineFlexibleALSource_H

#include "crossFlowTurbineALSource.H"
#include "actuatorFlexibleLineSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class crossFlowTurbineFlexibleALSource Declaration
\*---------------------------------------------------------------------------*/

class crossFlowTurbineFlexibleALSource
:
    public crossFlowTurbineALSource   // single-inheritance path; flexible
                                      // blade objects are created inside
                                      // blades_ via createBlades() override
{

protected:

    // Protected data

        //- Flexible actuator line blades (owned here; stored in
        //  blades_ (base class PtrList) as actuatorLineSource pointers but
        //  the concrete objects are actuatorFlexibleLineSource instances so
        //  virtual dispatch works correctly via the base pointer)
        //  No separate list is needed; the override of createBlades() puts
        //  actuatorFlexibleLineSource* into the inherited blades_ PtrList.


    // Protected Member Functions

        //- Override createBlades() to instantiate actuatorFlexibleLineSource
        //  objects instead of plain actuatorLineSource objects.
        virtual void createBlades();

        //- Rotate override: delegates to base class.  Re-declared here so
        //  the compiler always picks this class's rotate() when called
        //  from addSup(), avoiding any lingering ambiguity.
        virtual void rotate(scalar radians);


public:

    //- Runtime type information
    TypeName("crossFlowTurbineFlexibleALSource");


    // Constructors

        //- Construct from components
        crossFlowTurbineFlexibleALSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~crossFlowTurbineFlexibleALSource();


    // Member Functions

        // Source term addition
        // These are explicitly declared here to resolve the "no unique final
        // overrider" errors.  Each one delegates to
        // crossFlowTurbineALSource's implementation which in turn calls
        // the virtual addSup on each blade (which dispatches to
        // actuatorFlexibleLineSource::addSup because the concrete type is
        // actuatorFlexibleLineSource).

            //- Source term to momentum equation (incompressible)
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Source term to turbulence model equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );

        // I-O

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

            //- Write data (resolves writeData ambiguity from option base)
            virtual void writeData(Ostream& os) const;

            //- Print dictionary values
            virtual void printCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
