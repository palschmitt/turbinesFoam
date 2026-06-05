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

\*---------------------------------------------------------------------------*/

#include "crossFlowTurbineFlexibleALSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

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
    int nBlades = turbineFALSource::nBlades_;
    turbineFALSource::blades_.setSize(nBlades);
    int nElements;
    List<List<scalar> > elementData;
    List<List<scalar> > profileData;
    word modelType = "actuatorFlexibleLineSource";
    List<scalar> frontalAreas(nBlades); // Frontal area from each blade

    forAll(turbineFALSource::blades_, i)
    {
        word bladeName = turbineFALSource::bladeNames_[i];
        // Create dictionary items for this blade
        dictionary bladeSubDict = turbineFALSource::bladesDict_.subDict(bladeName);
        bladeSubDict.lookup("nElements") >> nElements;
        bladeSubDict.lookup("elementData") >> elementData;
        scalar azimuthalOffset = bladeSubDict.lookupOrDefault
        (
            "azimuthalOffset",
            0.0
        );

        bladeSubDict.add("freeStreamVelocity", turbineFALSource::freeStreamVelocity_);
        bladeSubDict.add("fieldNames", turbineFALSource::coeffs_.lookup("fieldNames"));
        bladeSubDict.add("profileData", turbineFALSource::profileData_);

        if (debug)
        {
            Info<< "Creating flexible actuator line blade " << bladeName << endl;
            Info<< "Blade has " << nElements << " elements" << endl;
            Info<< "Element data:" << endl;
            Info<< elementData << endl << endl;
        }

        // Convert element data into actuator line element geometry
        label nGeomPoints = elementData.size();
        List<List<List<scalar> > > elementGeometry(nGeomPoints);
        List<vector> initialVelocities(nGeomPoints, vector::zero);
        // Frontal area for this blade
        scalar frontalArea = 0.0;
        forAll(elementData, j)
        {
            // Read CFTAL dict data
            scalar axialDistance = elementData[j][0];
            scalar radius = elementData[j][1];
            scalar azimuthDegrees = elementData[j][2] + azimuthalOffset;
            scalar azimuthRadians = degToRad(azimuthDegrees);
            scalar chordLength = elementData[j][3];
            scalar chordMount = elementData[j][4];
            scalar pitch = elementData[j][5];

            // Compute frontal area contribution from this geometry segment
            if (j > 0)
            {
                scalar deltaAxial = axialDistance - elementData[j-1][0];
                scalar meanRadius = (radius + elementData[j-1][1])/2;
                frontalArea += mag(deltaAxial*meanRadius);
            }

            // Set sizes for actuatorFlexibleLineSource elementGeometry lists
            elementGeometry[j].setSize(9);
            elementGeometry[j][0].setSize(3);
            elementGeometry[j][1].setSize(3);
            elementGeometry[j][2].setSize(1);
            elementGeometry[j][3].setSize(3);
            elementGeometry[j][4].setSize(1);
            elementGeometry[j][5].setSize(1);
            elementGeometry[j][6].setSize(2);//E mu
            elementGeometry[j][7].setSize(5);//A        Iz       Iy          J alpha
            elementGeometry[j][8].setSize(6);//restraints
            
            // Create geometry point for AL source at origin
            vector point = turbineFALSource::origin_;
            // Move along axis
            point += axialDistance*turbineFALSource::axis_;
            // Move along chord according to chordMount
            scalar chordDisplacement = (chordMount - 0.25)*chordLength;
            point -= chordDisplacement*turbineFALSource::freeStreamDirection_;
            // Move along radial direction
            point += radius*turbineFALSource::radialDirection_;
            // Set initial velocity of quarter chord
            scalar radiusCorr = sqrt(magSqr((chordMount - 0.25)*chordLength)
                                     + magSqr(radius));
            vector initialVelocity = -turbineFALSource::freeStreamDirection_*turbineFALSource::omega_*radiusCorr;
            scalar velAngle = atan2(((chordMount - 0.25)*chordLength), radius);
            turbineFALSource::rotateVector(initialVelocity, vector::zero, turbineFALSource::axis_, velAngle);
            initialVelocities[j] = initialVelocity;
            // Rotate point and initial velocity according to azimuth value
            turbineFALSource::rotateVector(point, turbineFALSource::origin_, turbineFALSource::axis_, azimuthRadians);
            turbineFALSource::rotateVector
            (
                initialVelocities[j],
                vector::zero,
                turbineFALSource::axis_,
                azimuthRadians
            );

            // Set point coordinates for AL source
            elementGeometry[j][0][0] = point.x(); // x location of geom point
            elementGeometry[j][0][1] = point.y(); // y location of geom point
            elementGeometry[j][0][2] = point.z(); // z location of geom point

            // Set span directions for AL source
            elementGeometry[j][1][0] = turbineFALSource::axis_.x(); // x component of span dir
            elementGeometry[j][1][1] = turbineFALSource::axis_.y(); // y component of span dir
            elementGeometry[j][1][2] = turbineFALSource::axis_.z(); // z component of span dir

            // Set chord length
            elementGeometry[j][2][0] = chordLength;

            // Set chord reference direction
            vector chordDirection = -turbineFALSource::freeStreamDirection_;
            turbineFALSource::rotateVector(chordDirection, vector::zero, turbineFALSource::axis_, azimuthRadians);
            elementGeometry[j][3][0] = chordDirection.x();
            elementGeometry[j][3][1] = chordDirection.y();
            elementGeometry[j][3][2] = chordDirection.z();

            // Set chord mount
            elementGeometry[j][4][0] = chordMount;

            // Set pitch
            elementGeometry[j][5][0] = pitch;
                        //Element mats
            //E
            elementGeometry[j][6][0] = elementData[j][6];
            elementGeometry[j][6][1] = elementData[j][7];
            //elementGeometry[j][6][1] = mu;
            //A        Iz       Iy          J alpha
            //Needs rotation?
            elementGeometry[j][7][0] = elementData[j][8];
            elementGeometry[j][7][1] = elementData[j][9];
            elementGeometry[j][7][2] = elementData[j][10];
            elementGeometry[j][7][3] = elementData[j][11];
            elementGeometry[j][7][4] = elementData[j][12];
            //Restraints, fix hub
            
            if (j<1)
            {
                elementGeometry[j][8] = 1;
                }
            else
            {
            elementGeometry[j][8] = 0.;
        }
            
            
        }

        // Add frontal area to list
        frontalAreas[i] = frontalArea;

        if (debug)
        {
            Info<< "Converted element geometry:" << endl << elementGeometry
                << endl;
            Info<< "Frontal area from " << bladeName << ": " << frontalArea
                << endl;
        }

        bladeSubDict.add("elementGeometry", elementGeometry);
        bladeSubDict.add("initialVelocities", initialVelocities);
        bladeSubDict.add("dynamicStall", turbineFALSource::dynamicStallDict_);
        bladeSubDict.add
        (
            "addedMass",
            turbineFALSource::coeffs_.lookupOrDefault("addedMass", false)
        );
        bladeSubDict.add
        (
            "velocitySampleRadius",
            turbineFALSource::coeffs_.lookupOrDefault("velocitySampleRadius", 0.0)
        );
        bladeSubDict.add
        (
            "nVelocitySamples",
            turbineFALSource::coeffs_.lookupOrDefault("nVelocitySamples", 20)
        );
        bladeSubDict.add("selectionMode", turbineFALSource::coeffs_.lookup("selectionMode"));
        bladeSubDict.add("cellSet", turbineFALSource::coeffs_.lookup("cellSet"));

        // Lookup or create flowCurvature subDict
        dictionary fcDict = turbineFALSource::coeffs_.subOrEmptyDict("flowCurvature");
        fcDict.lookupOrAddDefault("active", true);
        word defaultFCModel = "Goude";
        fcDict.lookupOrAddDefault
        (
            "flowCurvatureModel",
            defaultFCModel
        );
        bladeSubDict.add("flowCurvature", fcDict);

        // Do not write force from individual actuator line unless specified
        bladeSubDict.lookupOrAddDefault("writeForceField", false);

        dictionary dict;
        dict.add("actuatorFlexibleLineSourceCoeffs", bladeSubDict);
        dict.add("type", "actuatorFlexibleLineSource");
        dict.add("active", turbineFALSource::dict_.lookup("active"));

        actuatorFlexibleLineSource* blade = new actuatorFlexibleLineSource
        (
            turbineFALSource::name_ + "." + bladeName,
            modelType,
            dict,
            turbineFALSource::mesh_
        );

        turbineFALSource::blades_.set(i, blade);
    }

    // Frontal area is twice the maximum blade frontal area
    turbineFALSource::frontalArea_ = 2*max(frontalAreas);
    Info<< "Frontal area of " << turbineFALSource::name_ << ": " << turbineFALSource::frontalArea_ << endl;
}





// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineFlexibleALSource::crossFlowTurbineFlexibleALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    crossFlowTurbineALSource(name, modelType, dict, mesh),
    turbineFALSource(name, modelType, dict, mesh)
{
    read(dict);
    turbineFALSource::createCoordinateSystem();
    turbineFALSource::createBlades();
    if (hasStruts_)
    {
        createStruts();
    }
    if (hasShaft_)
    {
        createShaft();
    }
    turbineFALSource::createOutputFile();

    // Rotate turbine to azimuthalOffset if necessary
    scalar azimuthalOffset = turbineFALSource::coeffs_.lookupOrDefault("azimuthalOffset", 0.0);
    rotate(degToRad(azimuthalOffset));

    if (debug)
    {
        Info<< "crossFlowTurbineFlexibleALSource created at time = " << turbineFALSource::time_.value()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineFlexibleALSource::~crossFlowTurbineFlexibleALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




bool Foam::fv::crossFlowTurbineFlexibleALSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        turbineFALSource::read(dict);

        // Get struts information
        strutsDict_ = turbineFALSource::coeffs_.subOrEmptyDict("struts");
        if (strutsDict_.keys().size() > 0)
        {
            turbineFALSource::hasStruts_ = true;
        }

        // Get shaft information
        shaftDict_ = turbineFALSource::coeffs_.subOrEmptyDict("shaft");
        if (shaftDict_.keys().size() > 0)
        {
            turbineFALSource::hasShaft_ = true;
        }

        if (debug)
        {
            Info<< "Debugging on" << endl;
            Info<< "Cross-flow turbine properties:" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
