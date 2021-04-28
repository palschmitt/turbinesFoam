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
    List<scalar> pitches(nGeometryPoints);
    List<scalar> chordMounts(nGeometryPoints);
    totalLength_ = 0.0;
    chordLength_ = 0.0;
	//FEA Data
    List<List<scalar>> materials(nGeometryPoints);
    List<List<scalar>> sections(nGeometryPoints);
    List<List<scalar>> FEArestraints(nGeometryPoints);


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
        materials[i] = elementGeometry_[i][6];
        sections[i] = elementGeometry_[i][7];
		FEArestraints[i] = elementGeometry_[i][8]; 
    }

    // Store blade root and tip locations for distance calculations
    rootLocation_ = points[0];
    tipLocation_ = points[nGeometryPoints - 1];

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
        Info<< "Root location: " << rootLocation_ << endl;
        Info<< "Tip location: " << tipLocation_ << endl;
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
        List<scalar> material(2);
        List<scalar> section(5);
        List<int> restraint(6,0.);
	
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

		//Linearly interpolate List element by element
		//Include dimension check?
		//material //scalar
		for (int j=0.; j<2;j++)
		{
		scalar material1 = materials[geometrySegmentIndex][j];
        scalar material2 = materials[geometrySegmentIndex + 1][j];
        scalar deltamaterialTotal = material1 - material2;
        material[j] = material1
                    + deltamaterialTotal/nElementsPerSegment*pointIndex
                    + deltamaterialTotal/nElementsPerSegment/2;
			}
			
		//Include dimension check?
		//Section Data //scalar
		for (int j=0.; j<5;j++)
		{
		scalar section1 = sections[geometrySegmentIndex][j];
        scalar section2 = sections[geometrySegmentIndex + 1][j];
        scalar deltasectionTotal = section1 - section2;
        section[j] = section1
                    + deltasectionTotal/nElementsPerSegment*pointIndex
                    + deltasectionTotal/nElementsPerSegment/2;
			}
			
		//Include dimension check? 
		//Do not interpolate restraint! Simply apply first
		//restraint int
		if (i==0)
			{
			for (int j=0.; j<6;j++)
			{
				restraint[j] = FEArestraints[geometrySegmentIndex][j];
			}						
			}
	
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
        scalar rootDistance = mag(position - rootLocation_)/totalLength_;

        // Create a dictionary for this actuatorFlexibleLineElement
        dictionary dict;
        dict.add("position", position);
        dictionary profileDataDict = profileData_.subDict(profileName);
        dict.add("profileData", profileDataDict);
        dict.add("profileName", profileName);
        dict.add("chordLength", chordLength);
        dict.add("chordDirection", chordDirection);
        dict.add("chordRefDirection", chordRefDirection);
        dict.add("spanLength", spanLength);
        dict.add("spanDirection", spanDirection);
        dict.add("freeStreamVelocity", freeStreamVelocity_);
        dict.add("chordMount", chordMount);
        dict.add("rootDistance", rootDistance);
        dict.add("FEAmaterial", material);
        dict.add("FEAsection", section);
        dict.add("FEArestraints", restraint);
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
            Info<< "Creating actuatorBernoulliLineElement: " << name << endl;
            Info<< "Geometry segment index: " << geometrySegmentIndex << endl;
            Info<< "Position: " << position << endl;
            Info<< "Chord length: " << chordLength << endl;
            Info<< "Chord direction (before pitching): " << chordDirection
                << endl;
            Info<< "Pitch (degrees): " << pitch << endl;
            Info<< "Span length: " << spanLength << endl;
            Info<< "Span direction: " << spanDirection << endl;
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



void Foam::fv::actuatorFlexibleLineSource::writePerf()
{
    scalar time = mesh_.time().value();
    scalar totalArea = 0.0;
    scalar x = 0.0;
    scalar y = 0.0;
    scalar z = 0.0;
    scalar relVelMag = 0.0;
    scalar alphaDeg = 0.0;
    scalar alphaGeom = 0.0;
    scalar cl = 0.0;
    scalar cd = 0.0;
    scalar cm = 0.0;

    forAll(elements_, i)
    {
        scalar area = elements_[i].chordLength()*elements_[i].spanLength();
        totalArea += area;
        vector pos = elements_[i].position();
        x += pos[0]; y += pos[1]; z += pos[2];
        relVelMag += mag(elements_[i].relativeVelocity())*area;
        alphaDeg += elements_[i].angleOfAttack()*area;
        alphaGeom += elements_[i].angleOfAttackGeom()*area;
        cl += elements_[i].liftCoefficient()*area;
        cd += elements_[i].dragCoefficient()*area;
        cm += elements_[i].momentCoefficient()*area;
    }

    x /= nElements_; y /= nElements_; z /= nElements_;
    relVelMag /= totalArea;
    alphaDeg /= totalArea;
    alphaGeom /= totalArea;
    cl /= totalArea; cd /= totalArea; cm /= totalArea;

    // write time,x,y,z,rel_vel_mag,alpha_deg,alpha_geom_deg,cl,cd,cm
    *outputFile_<< time << "," << x << "," << y << "," << z << "," << relVelMag
                << "," << alphaDeg << "," << alphaGeom << "," << cl << ","
                << cd << "," << cm << endl;
}


void Foam::fv::actuatorFlexibleLineSource::calcEndEffects()
{
    if (debug)
    {
        Info<< "Calculating end effects for " << name_ << endl;
    }

    scalar pi = Foam::constant::mathematical::pi;
    List<scalar> c(nElements_, 1.0); // Chord lengths
    List<scalar> alpha(nElements_, 0.1); // Geometric AoA in radians
    List<scalar> theta(nElements_); // Span distance rescaled on [0, pi]
    List<scalar> relVelMag(nElements_, 1.0);
    simpleMatrix<scalar> D(nElements_, 0.0, 0.1);
    List<scalar> A(nElements_); // Fourier coefficients
    List<scalar> circulation(nElements_);
    List<scalar> cl(nElements_);

    // Create lists from element parameters
    forAll(elements_, n)
    {
        theta[n] = elements_[n].rootDistance()*pi;
        c[n] = elements_[n].chordLength();
        //~ alpha[n] = Foam::degToRad(elements_[n].angleOfAttackGeom());
        //~ relVelMag[n] = mag(elements_[n].relativeVelocityGeom());
    }

    // Create D matrix
    forAll(elements_, i)
    {
        scalar n = i + 1;
        forAll(elements_, m)
        {
            D[m][i] = 2.0*totalLength_/(pi*c[m])*sin(n*theta[m])
                    + n*sin(n*theta[m]) / sin(theta[m]);
        }
        D.source()[i] = alpha[i];
    }
    A = D.solve();

    forAll(elements_, m)
    {
        scalar sumA = 0.0;
        forAll(elements_, i)
        {
            scalar n = i + 1;
            sumA += A[i]*sin(n*theta[m]);
        }
        circulation[m] = 2*totalLength_*relVelMag[m]*sumA;
        cl[m] = circulation[m]/(0.5*c[m]*relVelMag[m]);
    }

    // Set endEffectFactor for all elements
    List<scalar> factors = cl/Foam::max(cl);
    forAll(elements_, i)
    {
        elements_[i].setEndEffectFactor(factors[i]);
    }

    if (debug == 2)
    {
        Info<< "Debug output from actuatorFlexibleLineSource::calcEndEffects:" << endl;
        Info<< "theta: " << theta << endl;
        Info<< "A: " << A << endl;
        Info<< "c: " << c << endl;
        Info<< "D.source: " << D.source() << endl;
        Info<< "D: " << D << endl;
        Info<< "cl: " << cl << endl;
        Info<< "factors:" << factors << endl;
    }
}


void Foam::fv::actuatorFlexibleLineSource::harmonicPitching()
{
    // Pitch the actuator line if time has changed
    scalar t = mesh_.time().value();
    if (t != lastMotionTime_)
    {
        scalar omega = reducedFreq_*2*mag(freeStreamVelocity_)/chordLength_;
        scalar dt = mesh_.time().deltaT().value();
        scalar deltaPitch = degToRad(pitchAmplitude_)*(sin(omega*t)
                          - sin(omega*(t - dt)));
        pitch(deltaPitch);
        lastMotionTime_ = t;
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
//Run FEA analysis of the actuator line if time has changed
scalar t = mesh_.time().value();
if (t != lastMotionTime_)
{
/*///////////////Debugging Data///////////////////////////////////////

*///////////////////////////////////////////////////////////////////////			
//Create input data for FEA Analysis
Info <<"Number of elements " <<nElements_ <<endl;
List<List<scalar>> FEAnodes;
FEAnodes.resize(2*nElements_+1);
List<List<int>>	FEAelems;
FEAelems.resize(2*nElements_);
List<List<int>> FEArestraints;
FEArestraints.resize(2*nElements_+1);
List<List<scalar>>  FEAmats;
FEAmats.resize(2*nElements_);
List<List<scalar>>  FEAsects;
FEAsects.resize(2*nElements_);
List<List<scalar>>  FEAloads;
FEAloads.resize(2*nElements_+1);
List<List<scalar>>  FEAprescribed;
FEAprescribed.resize(2*nElements_+1);
List<scalar> SubList;//Ugly workaround
List<int> SubiList;//Ugly workaround

//Better set FEA material as material/rho? Or get rhoInf from ObjReg?
scalar rho=1000; //Density missing for incompressible cases? rhoref?
//
//AL Line element positions are nodes of FEA model
//Get position +- spanwidth*spandirection as nodes positions
//Info<< "AL Element position "<<elements_[0].position() <<endl;
//Info<< "AL Element spanDirection "<<elements_[0].spanDirection() <<endl;
//Info<< "AL Element spanLength "<<elements_[0].spanLength() <<endl;

vector Position=elements_[0].position()-0.5*elements_[0].spanLength()*elements_[0].spanDirection();
Info<< "First Node Pos "<<Position <<endl;
//vector Position=elements_[0].position();
SubList.clear();
SubList.resize(3);
SubList[0]=Position.x();
SubList[1]=Position.y();
SubList[2]=Position.z();
FEAnodes[0]=(SubList);

//No force at first node
SubList.clear();
SubList.resize(6);
SubList=0.;
FEAloads[0]=SubList;//Fluid force


//Using restraind defined at first and/or last node, impossible to define restrained in between!
FEArestraints[0]=elements_[0].FEArestraints();
FEArestraints[2*nElements_]=elements_[nElements_-1].FEArestraints();

SubList.resize(6);
SubList=0.;
FEAprescribed[2*nElements_]=SubList;

forAll(elements_, i)
{

//Data per node	
//Positions	
//Info<< "AL Element position "<<elements_[i].position() <<endl;
//Info<< "AL Element spanDirection "<<elements_[i].spanDirection() <<endl;
//Info<< "AL Element spanlength "<<elements_[i].spanLength() <<endl;
vector Position=elements_[i].position();
SubList.clear();
SubList.resize(3);
SubList=0.;
SubList[0]=Position.x();
SubList[1]=Position.y();
SubList[2]=Position.z();
FEAnodes[2*i+1]=SubList;

Position=elements_[i].position()+0.5*elements_[i].spanLength()*elements_[i].spanDirection();
SubList.clear();
SubList.resize(3);
SubList=0.;
SubList[0]=Position.x();
SubList[1]=Position.y();
SubList[2]=Position.z();
FEAnodes[2*i+2]=SubList;

//List of forces in XYZ
//Applying fluid force to FEAnode at center of element

SubList.clear();
SubList.resize(6);
SubList=0.;
SubList[0]=(elements_[i].force().x()-elements_[i].structforce().x())*rho;
SubList[1]=(elements_[i].force().y()-elements_[i].structforce().y())*rho;
SubList[2]=(elements_[i].force().z()-elements_[i].structforce().z())*rho;
//Dummy force for debugging
////vector DummyForce=vector(0.5, 0., 0.);
//vector DummyForce=vector(0.1, 0.1*sin(elements_[i].omega()*t), 0.1*cos(elements_[i].omega()*t));
//Info<< "Dummy Force: " <<DummyForce <<endl;
//SubList[0]=(DummyForce.x()-elements_[i].structforce().x())*rho;
//SubList[1]=(DummyForce.y()-elements_[i].structforce().y())*rho;
//SubList[2]=(DummyForce.z()-elements_[i].structforce().z())*rho;

FEAloads[2*i+1]=SubList;//Fluid force

//Empty node
SubList=0.;
FEAloads[2*i+2]=SubList;


SubList.clear();
SubList.resize(6);
SubList=0.;
FEAprescribed[2*i]=SubList;//Apply previous deformation?
FEAprescribed[2*i+1]=SubList;//Apply previous deformation?

//Save old force to only apply difference causing additional deformation
elements_[i].setStructForce(elements_[i].force());

//elements_[i].setStructForce(DummyForce);

//Data below per element
//Element definitions
SubiList.clear();
SubiList=0.;
SubiList.resize(2);
SubiList[0]=2*i;
SubiList[1]=2*i+1;
FEAelems[2*i]=SubiList;//FEA Element connections are simply Node(N) to Node(N+1);

SubiList[0]=2*i+1;
SubiList[1]=2*i+2;
FEAelems[2*i+1]=SubiList;//FEA Element connections are simply Node(N) to Node(N+1);

FEAmats[2*i]=elements_[i].FEAmaterial();// E  Poisson
FEAmats[2*i+1]=elements_[i].FEAmaterial();// E  Poisson
//Info <<"Material of element "<< i << " is "<<elements_[i].FEAmaterial() <<endl;
FEAsects[2*i]=elements_[i].FEAsects();//Section data A        Iz       Iy          J        alpha
FEAsects[2*i+1]=elements_[i].FEAsects();//Section data A        Iz       Iy          J        alpha
//Info <<"e Element Sections " <<elements_[i].FEAsects() <<endl;
}	
//*////////////////////////////////////////////////////////////////////////		
//Create FA and apply returned discplacement
Info<< "Input for Frame Analysis: "<< endl;
Info<< "FEAnodes: "<<FEAnodes<< endl;
//Info<< "FEAelems: "<<FEAelems<< endl;
//Info<< "FEArestraints: "<<FEArestraints<< endl;
//Info<< "FEAmats: "<<FEAmats<< endl;
//Info<< "FEAsects: "<<FEAsects<< endl;
Info<< "FEAloads: "<<FEAloads<< endl;
//Info<< "FEAprescribed: "<<FEAprescribed<< endl;


//Execute FEA simulation
FrameAnalysis FA(FEAnodes,FEAelems,FEArestraints, FEAmats,FEAsects,FEAloads,FEAprescribed);
Info<< "Deformation from FEA Analysis "<<FA.nodedispList()<< endl;

//Ugly data transfer, needs cleaning and proper access to FA data
List<List<scalar>> FEADeformation=FA.nodedispList();

//Create Vector List of new nodepositions
List<vector> NewFEANodepositions;
NewFEANodepositions.resize(2*nElements_+1);

forAll(FEADeformation,i)
{
NewFEANodepositions[i][0]=FEAnodes[i][0]+FEADeformation[i][0];
NewFEANodepositions[i][1]=FEAnodes[i][1]+FEADeformation[i][1];
NewFEANodepositions[i][2]=FEAnodes[i][2]+FEADeformation[i][2];
}

forAll(elements_,i)
{
//Update ALelement  positions to center FEA node

elements_[i].setPosition(NewFEANodepositions[2*i+1]);

//Update ALelement  spandirection and length
vector P1=vector(NewFEANodepositions[2*i]);
vector P2=vector(NewFEANodepositions[2*i+2]);
vector spanLength=P2-P1;
//Info<< "Updating Spanlength to "<<spanLength<< endl;

elements_[i].setSpanLength(mag(spanLength));
elements_[i].setSpanDirection(spanLength/mag(spanLength));
//Info<< "AoA "<< elements_[i].angleOfAttack()<< endl;
}
lastMotionTime_ = t;
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

void Foam::fv::actuatorFlexibleLineSource::printCoeffs() const
{
    // Print turbine properties
    Info<< "Actuator line properties:" << endl;
    Info<< "Profile data:" << endl;
    Info<< profileData_ << endl;
    Info<< "First item of element geometry:" << endl;
    Info<< elementGeometry_[0] << endl;
}


void Foam::fv::actuatorFlexibleLineSource::rotate
(
    vector rotationPoint,
    vector axis,
    scalar radians
)
{
    forAll(elements_, i)
    {
        elements_[i].rotate(rotationPoint, axis, radians, true);
    }
}

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

void Foam::fv::actuatorFlexibleLineSource::pitch(scalar radians)
{
    forAll(elements_, i)
    {
        elements_[i].pitch(radians);
    }
}


void Foam::fv::actuatorFlexibleLineSource::pitch(scalar radians, scalar chordFraction)
{
    forAll(elements_, i)
    {
        elements_[i].pitch(radians, chordFraction);
    }
}


void Foam::fv::actuatorFlexibleLineSource::translate(vector translationVector)
{
    forAll(elements_, i)
    {
        elements_[i].translate(translationVector);
    }
}


void Foam::fv::actuatorFlexibleLineSource::setSpeed
(
    vector point,
    vector axis,
    scalar omega
)
{
    forAll(elements_, i)
    {
        elements_[i].setSpeed(point, axis, omega);
    }
}


void Foam::fv::actuatorFlexibleLineSource::scaleVelocity(scalar scale)
{
    forAll(elements_, i)
    {
        elements_[i].scaleVelocity(scale);
    }
}


void Foam::fv::actuatorFlexibleLineSource::setOmega(scalar omega)
{
    forAll(elements_, i)
    {
        elements_[i].setOmega(omega);
    }
}


const Foam::vector& Foam::fv::actuatorFlexibleLineSource::force()
{
    return force_;
}


const Foam::volVectorField& Foam::fv::actuatorFlexibleLineSource::forceField()
{
    return forceField_;
}


PtrList<Foam::fv::actuatorBernoulliLineElement>& Foam::fv::actuatorFlexibleLineSource::elements()
{
    return elements_;
}


Foam::vector Foam::fv::actuatorFlexibleLineSource::moment(vector point)
{
    vector moment(vector::zero);
    forAll(elements_, i)
    {
        moment += elements_[i].moment(point);
    }

    if (debug)
    {
        Info<< "Moment on " << name_ << " about " << point << ": " << moment
            << endl;
    }

    return moment;
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
    // Write element Displacement
    vtkFilePtr_()
        << "VECTORS displacement double "<<nl;

        forAll(elements_, i)
        {
            vector eDisp (elements_[i].displacement());
            vtkFilePtr_()
                << eDisp[0]
                << " "
                << eDisp[1]
                << " "
                << eDisp[2]
                << nl;
        }

    // Add to the VTK sequence counter
    vtkFileSequence_++;
}


// ************************************************************************* //
