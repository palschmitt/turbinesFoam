/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.  (GPL v3)

\*---------------------------------------------------------------------------*/

#include "crossFlowTurbineFlexibleALSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(crossFlowTurbineFlexibleALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        crossFlowTurbineFlexibleALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::crossFlowTurbineFlexibleALSource::createBlades()
{
    // ---------------------------------------------------------------------
    // This method mirrors crossFlowTurbineALSource::createBlades() exactly,
    // except that the concrete element type created is
    // actuatorFlexibleLineSource instead of actuatorLineSource.
    //
    // The blade structural data (material, section, restraints) is read from
    // three extra columns appended to each elementData row:
    //   col 6 : [E  Poisson]       – FEA material
    //   col 7 : [A Iz Iy J alpha]  – FEA cross-section
    //   col 8 : [r1 r2 r3 r4 r5 r6] – nodal restraints (0/1 per DOF)
    // ---------------------------------------------------------------------

    int nBlades = nBlades_;
    blades_.setSize(nBlades);

    // modelType string that matches actuatorFlexibleLineSource's TypeName
    word modelType = "actuatorFlexibleLineSource";

    List<scalar> frontalAreas(nBlades, 0.0);

    forAll(blades_, i)
    {
        word bladeName = bladeNames_[i];
        dictionary bladeSubDict = bladesDict_.subDict(bladeName);

        int nElements;
        bladeSubDict.lookup("nElements") >> nElements;

        List<List<scalar>> elementData;
        bladeSubDict.lookup("elementData") >> elementData;

        scalar azimuthalOffset =
            bladeSubDict.lookupOrDefault("azimuthalOffset", 0.0);

        bladeSubDict.add("freeStreamVelocity", freeStreamVelocity_);
        bladeSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
        bladeSubDict.add("profileData", profileData_);

        // -- Convert elementData into actuatorFlexibleLineSource geometry --
        // Each geometry point carries 9 sub-lists (indices 0-8); the first 6
        // match the standard actuatorLineSource format, the last 3 are the
        // structural data required by actuatorFlexibleLineSource.

        label nGeomPoints = elementData.size();
        List<List<List<scalar>>> elementGeometry(nGeomPoints);
        List<vector> initialVelocities(nGeomPoints, vector::zero);
        scalar frontalArea = 0.0;

        forAll(elementData, j)
        {
            scalar axialDistance  = elementData[j][0];
            scalar radius         = elementData[j][1];
            scalar azimuthDegrees = elementData[j][2] + azimuthalOffset;
            scalar azimuthRadians = degToRad(azimuthDegrees);
            scalar chordLength    = elementData[j][3];
            scalar chordMount     = elementData[j][4];
            scalar pitch          = elementData[j][5];

            // Frontal area contribution
            if (j > 0)
            {
                scalar deltaAxial  = axialDistance - elementData[j-1][0];
                scalar meanRadius  = (radius + elementData[j-1][1]) / 2.0;
                frontalArea       += mag(deltaAxial * meanRadius);
            }

            // Allocate 9 sub-lists per geometry point
            elementGeometry[j].setSize(9);
            elementGeometry[j][0].setSize(3); // position
            elementGeometry[j][1].setSize(3); // span direction
            elementGeometry[j][2].setSize(1); // chord length
            elementGeometry[j][3].setSize(3); // chord reference direction
            elementGeometry[j][4].setSize(1); // chord mount
            elementGeometry[j][5].setSize(1); // pitch
            // Structural sub-lists – sizes match actuatorFlexibleLineSource
            elementGeometry[j][6].setSize(2); // FEA material [E, Poisson]
            elementGeometry[j][7].setSize(5); // FEA section  [A Iz Iy J alpha]
            elementGeometry[j][8].setSize(6); // FEA restraints [r1..r6]

            // Position
            vector point = origin_;
            point += axialDistance * axis_;
            scalar chordDisplacement = (chordMount - 0.25) * chordLength;
            point -= chordDisplacement * freeStreamDirection_;
            point += radius * radialDirection_;

            // Initial velocity at quarter chord
            scalar radiusCorr = sqrt
            (
                magSqr((chordMount - 0.25)*chordLength) + magSqr(radius)
            );
            vector initialVelocity = -freeStreamDirection_ * omega_ * radiusCorr;
            scalar velAngle = atan2((chordMount - 0.25)*chordLength, radius);
            rotateVector(initialVelocity, vector::zero, axis_, velAngle);
            initialVelocities[j] = initialVelocity;

            rotateVector(point, origin_, axis_, azimuthRadians);
            rotateVector(initialVelocities[j], vector::zero, axis_, azimuthRadians);

            // [0] position
            elementGeometry[j][0][0] = point.x();
            elementGeometry[j][0][1] = point.y();
            elementGeometry[j][0][2] = point.z();

            // [1] span direction (along turbine axis for blades)
            elementGeometry[j][1][0] = axis_.x();
            elementGeometry[j][1][1] = axis_.y();
            elementGeometry[j][1][2] = axis_.z();

            // [2] chord length
            elementGeometry[j][2][0] = chordLength;

            // [3] chord reference direction
            vector chordDirection = -freeStreamDirection_;
            rotateVector(chordDirection, vector::zero, axis_, azimuthRadians);
            elementGeometry[j][3][0] = chordDirection.x();
            elementGeometry[j][3][1] = chordDirection.y();
            elementGeometry[j][3][2] = chordDirection.z();

            // [4] chord mount
            elementGeometry[j][4][0] = chordMount;

            // [5] pitch
            elementGeometry[j][5][0] = pitch;

            // [6] FEA material – read from elementData columns 6-7
            // Expected: elementData[j][6] = E, elementData[j][7] = Poisson
            elementGeometry[j][6][0] =
                (elementData[j].size() > 6) ? elementData[j][6] : 2.1e11;
            elementGeometry[j][6][1] =
                (elementData[j].size() > 7) ? elementData[j][7] : 0.3;

            // [7] FEA section – read from elementData columns 8-12
            // Expected: A, Iz, Iy, J, alpha
            for (int k = 0; k < 5; k++)
                elementGeometry[j][7][k] =
                    (elementData[j].size() > 8 + k) ? elementData[j][8 + k] : 0.0;

            // [8] FEA restraints – read from elementData columns 13-18
            // Expected: 6 integers (0 or 1) per DOF
            for (int k = 0; k < 6; k++)
                elementGeometry[j][8][k] =
                    (elementData[j].size() > 13 + k) ? elementData[j][13 + k] : 0.0;
        }

        frontalAreas[i] = frontalArea;

        bladeSubDict.add("elementGeometry",  elementGeometry);
        bladeSubDict.add("initialVelocities", initialVelocities);
        bladeSubDict.add("dynamicStall",      dynamicStallDict_);
        bladeSubDict.add
        (
            "addedMass",
            coeffs_.lookupOrDefault("addedMass", false)
        );
        bladeSubDict.add
        (
            "velocitySampleRadius",
            coeffs_.lookupOrDefault("velocitySampleRadius", 0.0)
        );
        bladeSubDict.add
        (
            "nVelocitySamples",
            coeffs_.lookupOrDefault("nVelocitySamples", 20)
        );
        bladeSubDict.add("selectionMode", coeffs_.lookup("selectionMode"));
        bladeSubDict.add("cellSet",       coeffs_.lookup("cellSet"));

        dictionary fcDict = coeffs_.subOrEmptyDict("flowCurvature");
        fcDict.lookupOrAddDefault("active", true);
        word defaultFCModel = "Goude";
        fcDict.lookupOrAddDefault("flowCurvatureModel", defaultFCModel);
        bladeSubDict.add("flowCurvature", fcDict);

        bladeSubDict.lookupOrAddDefault("writeForceField", false);

        dictionary dict;
        dict.add("actuatorLineSourceCoeffs", bladeSubDict);
        dict.add("type", modelType);
        dict.add("active", dict_.lookup("active"));

        // Create the concrete flexible blade and store it as an
        // actuatorLineSource pointer (base class PtrList type).
        // Virtual dispatch through the pointer will call
        // actuatorFlexibleLineSource::addSup at run time.
        actuatorLineSource* blade = new actuatorFlexibleLineSource
        (
            name_ + "." + bladeName,
            modelType,
            dict,
            mesh_
        );

        blades_.set(i, blade);

        if (debug)
        {
            Info<< "Created flexible blade: " << bladeName << nl
                << "  frontalArea: " << frontalArea << endl;
        }
    }

    frontalArea_ = 2.0 * max(frontalAreas);
    Info<< "Frontal area of " << name_ << ": " << frontalArea_ << endl;
}


void Foam::fv::crossFlowTurbineFlexibleALSource::rotate(scalar radians)
{
    // Identical to crossFlowTurbineALSource::rotate() – re-implemented here
    // so this class always provides the unique final overrider.
    if (debug)
        Info<< "Rotating " << name_ << " by " << radians << " rad" << endl;

    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, axis_, radians);
        blades_[i].setSpeed(origin_, axis_, omega_);
    }

    if (hasStruts_)
    {
        forAll(struts_, i)
        {
            struts_[i].rotate(origin_, axis_, radians);
            struts_[i].setSpeed(origin_, axis_, omega_);
        }
    }

    if (hasShaft_)
    {
        shaft_->rotate(origin_, axis_, radians);
        shaft_->setSpeed(origin_, axis_, omega_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineFlexibleALSource::crossFlowTurbineFlexibleALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    // cellSetOption is the common virtual base; it must be initialised first
    // and only once, regardless of how many paths reach it.
    cellSetOption(name, modelType, dict, mesh),
    crossFlowTurbineALSource(name, modelType, dict, mesh)
    // Note: actuatorFlexibleLineSource is NOT a base class here.
    // Flexible blades are instantiated inside createBlades() as concrete
    // actuatorFlexibleLineSource objects stored through the inherited
    // blades_ PtrList<actuatorLineSource>.  This avoids the diamond problem
    // entirely without requiring virtual inheritance changes in other headers.
{
    // crossFlowTurbineALSource constructor calls read(), createCoordinateSystem()
    // and createBlades() – but its createBlades() creates rigid blades.
    // We override createBlades() here and re-run it to replace those blades
    // with flexible ones.
    //
    // Destroy the rigid blades created by the base constructor:
    blades_.clear();

    // Re-run the full turbine setup using this class's createBlades():
    createBlades();

    // Struts and shaft remain as rigid actuatorLineSource objects because
    // structural flexibility of struts / shaft is typically not required.
    // Override createStruts() / createShaft() if needed in future.
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineFlexibleALSource::~crossFlowTurbineFlexibleALSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fv::crossFlowTurbineFlexibleALSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Delegate entirely to the base-class implementation.
    // Because blades_ holds actuatorFlexibleLineSource objects,
    // the virtual dispatch inside the base addSup() calls
    // actuatorFlexibleLineSource::addSup() for each blade automatically.
    crossFlowTurbineALSource::addSup(eqn, fieldI);
}


void Foam::fv::crossFlowTurbineFlexibleALSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    crossFlowTurbineALSource::addSup(rho, eqn, fieldI);
}


void Foam::fv::crossFlowTurbineFlexibleALSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    crossFlowTurbineALSource::addSup(eqn, fieldI);
}


void Foam::fv::crossFlowTurbineFlexibleALSource::writeData(Ostream& os) const
{
    // Resolve the writeData ambiguity.  Delegate to cellSetOption which
    // provides the canonical implementation through the option hierarchy.
    cellSetOption::writeData(os);
}


bool Foam::fv::crossFlowTurbineFlexibleALSource::read(const dictionary& dict)
{
    return crossFlowTurbineALSource::read(dict);
}


void Foam::fv::crossFlowTurbineFlexibleALSource::printCoeffs() const
{
    crossFlowTurbineALSource::printCoeffs();
    Info<< "Blade type: flexible (FrameAnalysis)" << endl;
}


// ************************************************************************* //
