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

#include "actuatorFlexibleLineSource.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorFlexibleLineSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuatorFlexibleLineSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::fv::actuatorFlexibleLineSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {

        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Look up information in dictionary
        coeffs_.lookup("elementProfiles") >> elementProfiles_;
        profileData_ = coeffs_.subDict("profileData");
        coeffs_.lookup("elementGeometry") >> elementGeometry_;
        coeffs_.lookup("nElements") >> nElements_;
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
        endEffectsActive_ = coeffs_.lookupOrDefault("endEffects", false);

        // Read harmonic pitching parameters if present
        dictionary pitchDict = coeffs_.subOrEmptyDict("harmonicPitching");
        harmonicPitchingActive_ = pitchDict.lookupOrDefault("active", false);
        reducedFreq_ = pitchDict.lookupOrDefault("reducedFreq", 0.0);
        pitchAmplitude_ = pitchDict.lookupOrDefault("amplitude", 0.0);

        // Read option for writing forceField
        bool writeForceField = coeffs_.lookupOrDefault
        (
            "writeForceField",
            true
        );
        if (not writeForceField)
        {
            forceField_.writeOpt() = IOobject::NO_WRITE;
        }

        if (debug)
        {
            Info<< "Debugging for actuatorFlexibleLineSource on" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::actuatorFlexibleLineSource::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/actuatorFlexibleLines"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/actuatorFlexibleLines"
            / mesh_.time().timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_<< "time,x,y,z,rel_vel_mag,alpha_deg,alpha_geom_deg,cl,cd,cm"
                << endl;
}

void Foam::fv::actuatorFlexibleLineSource::createOutputDir()
{
    if (Pstream::parRun())
    {
        vtkDir_ = mesh_.time().path()/"../postProcessing/actuatorFlexibleLines"
            / mesh_.time().timeName()/"VTK";
    }
    else
    {
        vtkDir_ = mesh_.time().path()/"postProcessing/actuatorFlexibleLines"
            / mesh_.time().timeName()/"VTK";
    }

    if (not isDir(vtkDir_))
    {
        Foam::mkDir(vtkDir_);
    }
}
void Foam::fv::actuatorFlexibleLineSource::createInitialElements()
{
    elements_.setSize(nElements_);

    label nGeometryPoints = elementGeometry_.size();
    label nGeometrySegments = nGeometryPoints - 1;
    label nElementsPerSegment = nElements_/nGeometrySegments;
    if (nElements_ % nGeometrySegments)
    {
        // Need to have integer number of elements per geometry segment
        FatalErrorIn("void actuatorFlexibleLineSource::createElements()")
            << "Number of actuator line elements must be multiple of the "
            << "number of actuator line geometry segments"
            << abort(FatalError);
    }
    List<vector> points(nGeometryPoints);
    List<vector> spanDirs(nGeometryPoints);
    List<scalar> chordLengths(nGeometryPoints);
    List<scalar> spanLengths(nGeometrySegments);
    List<vector> chordRefDirs(nGeometryPoints);
    List<vector> stiffnesses(nGeometryPoints);
    List<scalar> pitches(nGeometryPoints);
    List<scalar> chordMounts(nGeometryPoints);
    totalLength_ = 0.0;
    chordLength_ = 0.0;

    forAll(points, i)
    {
        // Extract geometry point
        scalar x = elementGeometry_[i][0][0];
        scalar y = elementGeometry_[i][0][1];
        scalar z = elementGeometry_[i][0][2];
        points[i] = vector(x, y, z);
               
        if (i > 0)
        {
            spanLengths[i - 1] = mag(points[i] - points[i-1]);
            totalLength_ += spanLengths[i - 1];
        }
        // Read span direction
        x = elementGeometry_[i][1][0];
        y = elementGeometry_[i][1][1];
        z = elementGeometry_[i][1][2];
        spanDirs[i] = vector(x, y, z);
        // Read chord length
        chordLengths[i] = elementGeometry_[i][2][0];
        chordLength_ += chordLengths[i];
        // Read chord ref dir
        x = elementGeometry_[i][3][0];
        y = elementGeometry_[i][3][1];
        z = elementGeometry_[i][3][2];
        chordRefDirs[i] = vector(x, y, z);
        // Read chord mount
        chordMounts[i] = elementGeometry_[i][4][0];
        // Read pitch
        pitches[i] = elementGeometry_[i][5][0];
        stiffnesses[i] = vector(elementGeometry_[i][6][0],elementGeometry_[i][6][1],elementGeometry_[i][6][2]);
    }

    // Store blade root and tip locations for distance calculations
    vector rootLocation = points[0];
    vector tipLocation = points[nGeometryPoints - 1];

    // Compute average chord length
    chordLength_ /= nGeometryPoints;

    // Compute aspect ratio
    aspectRatio_ = totalLength_/chordLength_;

    // Lookup initial element velocities if present
    List<vector> initialVelocities(nGeometryPoints, vector::zero);
    coeffs_.readIfPresent("initialVelocities", initialVelocities);

    if (debug)
    {
        Info<< "Total length: " << totalLength_ << endl;
        Info<< "Elements per geometry segment: " << nElementsPerSegment
            << endl;
        Info<< "Points:" << endl << points << endl;
        Info<< "Span directions:" << endl << spanDirs << endl;
        Info<< "Span lengths: " << endl << spanLengths << endl;
        Info<< "Chord lengths:" << endl << chordLengths << endl;
        Info<< "Pitches:" << endl << pitches << endl;
        Info<< "Root location: " << rootLocation << endl;
        Info<< "Tip location: " << tipLocation << endl;
    }

    forAll(elements_, i)
    {
        std::stringstream ss;
        ss << i;
        string str = ss.str();
        const word name = name_ + ".element" + str;

        // Actuator point geometry to be calculated from elementGeometry
        label geometrySegmentIndex = i/nElementsPerSegment;
        label pointIndex = i % nElementsPerSegment;
        label elementProfileIndex = i*elementProfiles_.size()/nElements_;
        word profileName = elementProfiles_[elementProfileIndex];
        vector position;
        scalar chordLength;
        vector chordDirection;
        vector chordRefDirection;
        //Need matrix?
        vector stiffness;
        scalar spanLength = spanLengths[geometrySegmentIndex];
        spanLength /= nElementsPerSegment;
        vector spanDirection;
        scalar pitch;
        scalar chordMount;
        vector initialVelocity;

        // Linearly interpolate position
        vector point1 = points[geometrySegmentIndex];
        vector point2 = points[geometrySegmentIndex + 1];
        vector segment = point2 - point1;
        position = point1
                 + segment/nElementsPerSegment*pointIndex
                 + segment/nElementsPerSegment/2;

		//Linearly interpolate stiffness
        vector stiffness1 = stiffnesses[geometrySegmentIndex];
        vector stiffness2 = stiffnesses[geometrySegmentIndex + 1];
        vector deltaStiffnessTotal = stiffness2 - stiffness1;
        stiffness = stiffness1
                    + deltaStiffnessTotal/nElementsPerSegment*pointIndex
                    + deltaStiffnessTotal/nElementsPerSegment/2;

        // Linearly interpolate chordLength
        scalar chordLength1 = chordLengths[geometrySegmentIndex];
        scalar chordLength2 = chordLengths[geometrySegmentIndex + 1];
        scalar deltaChordTotal = chordLength2 - chordLength1;
        chordLength = chordLength1
                    + deltaChordTotal/nElementsPerSegment*pointIndex
                    + deltaChordTotal/nElementsPerSegment/2;

        // Linearly interpolate spanDirection
        vector spanDir1 = spanDirs[geometrySegmentIndex];
        vector spanDir2 = spanDirs[geometrySegmentIndex + 1];
        vector deltaSpanTotal = spanDir2 - spanDir1;
        spanDirection = spanDir1
                      + deltaSpanTotal/nElementsPerSegment*pointIndex
                      + deltaSpanTotal/nElementsPerSegment/2;

        // Linearly interpolate section pitch
        scalar pitch1 = pitches[geometrySegmentIndex];
        scalar pitch2 = pitches[geometrySegmentIndex + 1];
        scalar deltaPitchTotal = pitch2 - pitch1;
        pitch = pitch1
              + deltaPitchTotal/nElementsPerSegment*pointIndex
              + deltaPitchTotal/nElementsPerSegment/2;

        // Linearly interpolate chord mount
        scalar cm1 = chordMounts[geometrySegmentIndex];
        scalar cm2 = chordMounts[geometrySegmentIndex + 1];
        scalar deltaCmTotal = cm2 - cm1;
        chordMount = cm1 + deltaCmTotal/nElementsPerSegment*pointIndex
                   + deltaCmTotal/nElementsPerSegment/2;

        // Linearly interpolate element velocity
        vector vel1 = initialVelocities[geometrySegmentIndex];
        vector vel2 = initialVelocities[geometrySegmentIndex + 1];
        vector deltaVelTotal = vel2 - vel1;
        initialVelocity = vel1
                        + deltaVelTotal/nElementsPerSegment*pointIndex
                        + deltaVelTotal/nElementsPerSegment/2;

        // Linearly interpolate chordDirection
        vector chordDir1 = chordRefDirs[geometrySegmentIndex];
        vector chordDir2 = chordRefDirs[geometrySegmentIndex + 1];
        vector deltaChordDirTotal = chordDir2 - chordDir1;
        chordDirection = chordDir1
                       + deltaChordDirTotal/nElementsPerSegment*pointIndex
                       + deltaChordDirTotal/nElementsPerSegment/2;

        // Chord reference direction (before pitching)
        chordRefDirection = chordDirection;
        
        // Calculate nondimensional root distance
        scalar rootDistance = mag(position - rootLocation)/totalLength_;

        // Create a dictionary for this actuatorFlexibleLineElement
        dictionary dict;
        dict.add("position", position);
        dictionary profileDataDict = profileData_.subDict(profileName);
        dict.add("profileData", profileDataDict);
        dict.add("profileName", profileName);
        dict.add("chordLength", chordLength);
        dict.add("stiffness", stiffness);
        dict.add("chordDirection", chordDirection);
        dict.add("chordRefDirection", chordRefDirection);
        dict.add("spanLength", spanLength);
        dict.add("spanDirection", spanDirection);
        dict.add("freeStreamVelocity", freeStreamVelocity_);
        dict.add("chordMount", chordMount);
        dict.add("rootDistance", rootDistance);
        dict.add("addedMass", coeffs_.lookupOrDefault("addedMass", false));
        dict.add
        (
            "velocitySampleRadius",
            coeffs_.lookupOrDefault("velocitySampleRadius", 0.0)
        );
        dict.add
        (
            "nVelocitySamples",
            coeffs_.lookupOrDefault("nVelocitySamples", 20)
        );
        if (coeffs_.found("dynamicStall"))
        {
            dictionary dsDict = coeffs_.subDict("dynamicStall");
            dsDict.add("chordLength", chordLength);
            dict.add("dynamicStall", dsDict);
        }
        dictionary fcDict = coeffs_.subOrEmptyDict("flowCurvature");
        dict.add("flowCurvature", fcDict);
        bool writeElementPerf
        (
            coeffs_.lookupOrDefault("writeElementPerf", false)
        );
        dict.add("writePerf", writeElementPerf);

        if (debug)
        {
            Info<< "Creating actuatorLineElement: " << name << endl;
            Info<< "Geometry segment index: " << geometrySegmentIndex << endl;
            Info<< "Position: " << position << endl;
            Info<< "Chord length: " << chordLength << endl;
            Info<< "Chord direction (before pitching): " << chordDirection
                << endl;
            Info<< "Pitch (degrees): " << pitch << endl;
            Info<< "Span length: " << spanLength << endl;
            Info<< "Span direction: " << spanDirection << endl;
            Info<< "Stiffness: " << stiffness << endl;
            Info<< "Profile name index: " << elementProfileIndex << endl;
            Info<< "Profile name: " << profileName << endl;
            Info<< "writePerf: " << writeElementPerf << endl;
            Info<< "Root distance (nondimensional): " << rootDistance << endl;
        }

        actuatorBernoulliLineElement* element = new actuatorBernoulliLineElement
        (
            name, dict, mesh_
        );
        elements_.set(i, element);
        pitch = Foam::degToRad(pitch);
        elements_[i].pitch(pitch);
        elements_[i].setVelocity(initialVelocity);
    }
}







scalar Foam::fv::actuatorFlexibleLineSource::Cantileverdeflection(scalar x, scalar l, scalar F, scalar EI)
{
				scalar chordflexcontrib =0.;
		    
			chordflexcontrib=F*pow((l-x),3)/(6*EI)*(2-3*(l-x)/l+pow((l-x),3)/pow(l,2));
			
			if (x>l)
			{
				chordflexcontrib+=(x-l)*(F*pow(l,2)/(2*EI));
			}
	return chordflexcontrib;
	}

void Foam::fv::actuatorFlexibleLineSource::evaluateDeformation()
{

	//Nothing if no forces to avoid segfaults at startup
	if (mag(force())>SMALL)
	{
	// Move geometric positions for elements according to some beam theory or similar
	//Assuming 0DoF for element 0!!!!
		forAll(elements_, i)
		{
			//Ignoring variation of EI along length
			//Deformation at position of element i due to all other elements forces
			vector delta(0,0,0);
			//Spanwise distance
			scalar x = mag(elements_[i].position()-elements_[0].position());	
			
			forAll(elements_,j)  
			{
				if ((i>0) && (j>0))
				{		
				//Linear superpositions of Deformation contributions from all forces on all elements
				vector ResForce=elements_[j].force() - elements_[j].structforce();			
				//Contribution in chord direction	
				scalar Chorddirectionforce = mag(elements_[j].chordDirection()
				* (ResForce & elements_[j].chordDirection())
				/ magSqr(elements_[j].chordDirection()));
				//Contribution in normal to chord direction	
				//Info << "planformNormal"<<elements_[j].planformNormal()<<endl;
				//Info << "ResForce"<<ResForce<<endl;
				//Info << "Force"<<elements_[j].force()<<endl;
				//Info << "EI1"<<elements_[j].stiffness()[0]<<endl;
				//Info << "EI2"<<elements_[j].stiffness()[1]<<endl;
				
				scalar ChordNormalforce = mag(elements_[j].planformNormal()
				* (ResForce & elements_[j].planformNormal())
				/ magSqr(elements_[j].planformNormal()));
				
				//Stiffness is not integrated correctly along blade, rather use FEA?
				scalar l = mag(elements_[j].position()-elements_[0].position());	
				scalar chordflexcontrib = Cantileverdeflection(x,l,Chorddirectionforce,elements_[j].stiffness()[0]);	
				scalar chordnormalflexcontrib = Cantileverdeflection(x,l,ChordNormalforce,elements_[j].stiffness()[1]);	
				
				//Transform back into global coordinates
				//Shorten code  if forces are normalised?
				vector contrib=	elements_[j].planformNormal()*chordnormalflexcontrib/mag(elements_[j].planformNormal())+
								elements_[j].chordDirection()*chordflexcontrib/mag(elements_[j].chordDirection());
		
				delta+=contrib;
				
				//Save total force applied to each element, alternative store original position?
				elements_[i].setStructForce(elements_[i].structforce()-elements_[i].force());
				}

			
			}
			if (debug)
			{
				Info<<"Delta for element "<<i <<" is "<<delta<<endl;
				Info<<"elements_[i].structforce()"<<elements_[i].structforce()<<endl;
				Info<<"elements_[i].force()"<<elements_[i].force()<<endl;
				}
			elements_[i].translate(delta);	
		}
	}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::fv::actuatorFlexibleLineSource::actuatorFlexibleLineSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
actuatorLineSource(name,    modelType,    dict,    mesh),
    cellSetOption(name, modelType, dict, mesh),
    force_(vector::zero),
    forceField_
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "force",
            dimForce/dimVolume,
            vector::zero
        )
    ),
    writePerf_(coeffs_.lookupOrDefault("writePerf", false)),
    writeVTK_(coeffs_.lookupOrDefault("writeVTK", false)),
    vtkFileSequence_(0),
    vtkFilePtr_(NULL),
    lastMotionTime_(mesh.time().value()),
    endEffectsActive_(false)
{
    read(dict_);
    //First iteration to create elements and evaluate forces
    createInitialElements();
    //Maybe later in loop to find some equilibrium
    //Evaluate deformation from forces and change geometry
    evaluateDeformation();
    //Recreate Elements 
    if (writePerf_)
    {
        createOutputFile();
    }
    if (writeVTK_)
    {
        createOutputDir();
    }
    if (forceField_.writeOpt() == IOobject::AUTO_WRITE)
    {
        forceField_.write();
    }
    // Calculate end effects
    if (endEffectsActive_)
    {
        calcEndEffects();
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorFlexibleLineSource::~actuatorFlexibleLineSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar Foam::fv::actuatorFlexibleLineSource::foeppl
(
    scalar x,
    scalar pos,
    scalar exp
)
{
if (x>pos)
	{
	return	0;
	}
	else
	{
	return	pow(x-pos,exp);
	}
}


void Foam::fv::actuatorFlexibleLineSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // If harmonic pitching is active, do harmonic pitching
    if (harmonicPitchingActive_)
    {
        harmonicPitching();
    }
	evaluateDeformation();

    // Zero out force field
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);

    // Zero the total force vector
    force_ = vector::zero;

    forAll(elements_, i)
    {
        elements_[i].addSup(eqn, forceField_);
        force_ += elements_[i].force();
    }

    Info<< "Force (per unit density) on " << name_ << ": "
        << endl << force_ << endl << endl;

    // Check dimensions on force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    // Add source to eqn
    eqn += forceField_;

    // Write performance to file
    if (writePerf_ and Pstream::master())
    {
        writePerf();
    }

    // Write the VTK file
    if
    (
        writeVTK_ &&
        mesh_.time().outputTime() &&
        Pstream::master()
    )
    {
        writeVTK();
    }
}


void Foam::fv::actuatorFlexibleLineSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // If harmonic pitching is active, do harmonic pitching
    if (harmonicPitchingActive_)
    {
        harmonicPitching();
    }
	evaluateDeformation();

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    word fieldName = fieldNames_[fieldI];

    Info<< endl << "Adding " << fieldName << " from " << name_ << endl << endl;
    forAll(elements_, i)
    {
        elements_[i].calculateForce(U);
        elements_[i].addTurbulence(eqn, fieldName);
    }
}


void Foam::fv::actuatorFlexibleLineSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // If harmonic pitching is active, do harmonic pitching
    if (harmonicPitchingActive_)
    {
        harmonicPitching();
    }
	evaluateDeformation();
    // Zero out force field
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);

    // Zero the total force vector
    force_ = vector::zero;

    forAll(elements_, i)
    {
        elements_[i].addSup(rho, eqn, forceField_);
        force_ += elements_[i].force();
    }

    Info<< "Force on " << name_ << ": " << endl << force_ << endl << endl;

    // Check dimensions of force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    // Add source to eqn
    eqn += forceField_;

    // Write performance to file
    if (writePerf_ and Pstream::master())
    {
        writePerf();
    }

    // Write the VTK file
    if
    (
        writeVTK_ &&
        mesh_.time().outputTime() &&
        Pstream::master()
    )
    {
        writeVTK();
    }

}

void Foam::fv::actuatorFlexibleLineSource::writeVTK()
{

    fileName vtkFileName;

    // Pad the integer name for the VTK reader
    std::ostringstream cfc;
    cfc     << std::setw(12)
        << std::setfill('0')
        << vtkFileSequence_;

    // Construct file name
    vtkFileName = vtkDir_+"/"+name_+"_"+cfc.str()+".vtk";

    // Reset the file pointer to the soon to be written vtk
    vtkFilePtr_.reset
    (
        new OFstream(vtkFileName)
    );

    // Write header and time
    vtkFilePtr_()
        << "# vtk DataFile Version 3.0" << nl
        << "actuator line "<<name_<< nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    // Write Points
    vtkFilePtr_()
        << "POINTS "<<elements_.size()<<" double"<< nl;

        forAll(elements_, i)
        {
            vector ePosition (elements_[i].position());
            vtkFilePtr_()
                << ePosition[0]
                << " "
                << ePosition[1]
                << " "
                << ePosition[2]
                << nl;
        }


    // Write lines connecting nodes
    vtkFilePtr_()<< "LINES 1 "<<elements_.size()+1<< nl;
    vtkFilePtr_()<<elements_.size()<<" ";

    for (int i = 0; i<elements_.size();i++)
    {
        vtkFilePtr_()<<i<<" ";
    }
    vtkFilePtr_() << nl;
    vtkFilePtr_() << endl;

    // Tell VTK there is element data next
    vtkFilePtr_()
        << nl
        << "POINT_DATA "<<elements_.size()<<nl;

    // Write element velocity
    vtkFilePtr_()
        << "VECTORS Velocity double "<<nl;

        forAll(elements_, i)
        {
            vector eVel (elements_[i].velocity());
            vtkFilePtr_()
                << eVel[0]
                << " "
                << eVel[1]
                << " "
                << eVel[2]
                << nl;
        }

    vtkFilePtr_() << endl;

    // Write element force
    vtkFilePtr_()
        << "VECTORS Force double "<<nl;

        forAll(elements_, i)
        {
            vector eForce (elements_[i].force());
            vtkFilePtr_()
                << eForce[0]
                << " "
                << eForce[1]
                << " "
                << eForce[2]
                << nl;
        }


    // Write element stiffness
    vtkFilePtr_()
        << "VECTORS Stiffness double "<<nl;

        forAll(elements_, i)
        {
            vector eStiffness (elements_[i].stiffness());
            vtkFilePtr_()
                << eStiffness[0]
                << " "
                << eStiffness[1]
                << " "
                << eStiffness[2]
                << nl;
        }
    vtkFilePtr_() << endl;

    // Add to the VTK sequence counter
    vtkFileSequence_++;
}


// ************************************************************************* //
