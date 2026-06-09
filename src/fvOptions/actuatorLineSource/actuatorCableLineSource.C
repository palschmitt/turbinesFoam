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

#include "actuatorCableLineSource.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorCableLineSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuatorCableLineSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::actuatorCableLineSource::createInitialElements()
{
    elements_.setSize(nElements_);

    label nGeometryPoints   = elementGeometry_.size();
    label nGeometrySegments = nGeometryPoints - 1;
    label nElementsPerSegment = nElements_ / nGeometrySegments;

    if (nElements_ % nGeometrySegments)
    {
        FatalErrorIn("actuatorCableLineSource::createInitialElements()")
            << "Number of actuator line elements must be a multiple of the "
            << "number of geometry segments."
            << abort(FatalError);
    }

    // --- Extract geometry-point data ---
    List<vector> points(nGeometryPoints);
    List<vector> spanDirs(nGeometryPoints);
    List<scalar> chordLengths(nGeometryPoints);   // used as effective diameter
    List<scalar> spanLengths(nGeometrySegments);

    // Cable-specific per-point data
    // elementGeometry_[i][6] = (EA, pretension)
    // elementGeometry_[i][7] = (Cd)
    // elementGeometry_[i][8] = (r1 r2 r3)  restraints
    List<List<scalar>> cableMats(nGeometryPoints);  // [EA, T0]
    List<List<scalar>> cableSects(nGeometryPoints); // [Cd]
    List<List<scalar>> cableRestraints(nGeometryPoints); // [r1 r2 r3]

    totalLength_ = 0.0;
    chordLength_ = 0.0;

    forAll(points, i)
    {
        scalar x = elementGeometry_[i][0][0];
        scalar y = elementGeometry_[i][0][1];
        scalar z = elementGeometry_[i][0][2];
        points[i] = vector(x, y, z);

        if (i > 0)
        {
            spanLengths[i-1] = mag(points[i] - points[i-1]);
            totalLength_ += spanLengths[i-1];
        }

        x = elementGeometry_[i][1][0];
        y = elementGeometry_[i][1][1];
        z = elementGeometry_[i][1][2];
        spanDirs[i] = vector(x, y, z);

        chordLengths[i]    = elementGeometry_[i][2][0];
        chordLength_      += chordLengths[i];

        cableMats[i]       = elementGeometry_[i][6]; // [EA, T0]
        cableSects[i]      = elementGeometry_[i][7]; // [Cd]
        cableRestraints[i] = elementGeometry_[i][8]; // [r1 r2 r3]
    }

    rootLocation_ = points[0];
    tipLocation_  = points[nGeometryPoints - 1];
    chordLength_ /= nGeometryPoints;
    aspectRatio_  = totalLength_ / chordLength_;

    List<vector> initialVelocities(nGeometryPoints, vector::zero);
    coeffs_.readIfPresent("initialVelocities", initialVelocities);

    if (debug)
    {
        Info<< "Cable total length  : " << totalLength_ << endl;
        Info<< "Elements per segment: " << nElementsPerSegment << endl;
        Info<< "Root location       : " << rootLocation_ << endl;
        Info<< "Tip  location       : " << tipLocation_ << endl;
    }

    // --- Create individual elements ---
    forAll(elements_, i)
    {
        std::stringstream ss;
        ss << i;
        const word elemName = name_ + ".element" + ss.str();

        label seg  = i / nElementsPerSegment;
        label pt   = i % nElementsPerSegment;

        // Interpolated nodal positions
        vector point1 = points[seg];
        vector point2 = points[seg + 1];
        vector segment = point2 - point1;

        vector position = point1
                        + segment / nElementsPerSegment * pt
                        + segment / nElementsPerSegment / 2.0;
        vector P1 = point1 + segment / nElementsPerSegment * pt;
        vector P2 = point1 + segment / nElementsPerSegment * (pt + 1);

        // Interpolated chord length (effective diameter)
        scalar cl1 = chordLengths[seg];
        scalar cl2 = chordLengths[seg + 1];
        scalar chordLength = cl1
                           + (cl2 - cl1) / nElementsPerSegment * pt
                           + (cl2 - cl1) / nElementsPerSegment / 2.0;

        // Interpolated span direction
        vector sd1 = spanDirs[seg];
        vector sd2 = spanDirs[seg + 1];
        vector spanDirection = sd1
                             + (sd2 - sd1) / nElementsPerSegment * pt
                             + (sd2 - sd1) / nElementsPerSegment / 2.0;

        scalar spanLength = spanLengths[seg] / nElementsPerSegment;

        // Interpolated cable material [EA, T0]
        List<scalar> cableMat(2, 0.0);
        for (int j = 0; j < 2; j++)
        {
            scalar v1 = cableMats[seg][j];
            scalar v2 = cableMats[seg + 1][j];
            cableMat[j] = v1
                        + (v2 - v1) / nElementsPerSegment * pt
                        + (v2 - v1) / nElementsPerSegment / 2.0;
        }

        // Interpolated cable section [Cd]
        scalar Cd1 = cableSects[seg][0];
        scalar Cd2 = cableSects[seg + 1][0];
        scalar Cd  = Cd1
                   + (Cd2 - Cd1) / nElementsPerSegment * pt
                   + (Cd2 - Cd1) / nElementsPerSegment / 2.0;

        // Restraints: only at first and last node, zero elsewhere
        List<int> restraint(3, 0);
        if (i == 0)
            for (int j = 0; j < 3; j++)
                restraint[j] = int(cableRestraints[seg][j]);
        if (i == (nElements_ - 1))
            for (int j = 0; j < 3; j++)
                restraint[j] = int(cableRestraints[seg + 1][j]);

        // Interpolated element velocity
        vector vel1 = initialVelocities[seg];
        vector vel2 = initialVelocities[seg + 1];
        vector initialVelocity = vel1
                               + (vel2 - vel1) / nElementsPerSegment * pt
                               + (vel2 - vel1) / nElementsPerSegment / 2.0;

        // Build dictionary for this element
        dictionary dict;
        dict.add("position",        position);
        dict.add("P1",              P1);
        dict.add("P2",              P2);
        dict.add("chordLength",     chordLength);
        dict.add("chordDirection",  spanDirection); // dummy for base class
        dict.add("chordRefDirection", spanDirection);
        dict.add("spanLength",      spanLength);
        dict.add("spanDirection",   spanDirection);
        dict.add("freeStreamVelocity", freeStreamVelocity_);
        dict.add("rootDistance",
            mag(position - rootLocation_) / totalLength_);
        dict.add("CableEA",           cableMat[0]);
        dict.add("CablePretension",   cableMat[1]);
        dict.add("CableDragCoeff",    Cd);
        // cableSects[seg][1] = cable diameter (if provided); fall back
        // to chord length (effective diameter) if not set.
        scalar elemDiam = (cableSects[seg].size() > 1 && cableSects[seg][1] > VSMALL)
                        ? (cableSects[seg][1]
                           + (cableSects[seg+1][1] - cableSects[seg][1])
                             / nElementsPerSegment * pt
                           + (cableSects[seg+1][1] - cableSects[seg][1])
                             / nElementsPerSegment / 2.0)
                        : chordLength;
        // cableSects[seg][2] = cable material density [kg/m3]
        scalar elemRhoC = (cableSects[seg].size() > 2)
                        ? cableSects[seg][2] : 7850.0;
        // cableSects[seg][3] = reference fluid density [kg/m3]
        scalar elemRhoF = (cableSects[seg].size() > 3)
                        ? cableSects[seg][3] : 1025.0;
        dict.add("CableDiameter",    elemDiam);
        dict.add("CableDensity",     elemRhoC);
        dict.add("CableFluidDensity", elemRhoF);
        dict.add("CableRestraints",  restraint);
        dict.add("addedMass", false);
        dict.add("velocitySampleRadius",
            coeffs_.lookupOrDefault("velocitySampleRadius", 0.0));
        dict.add("nVelocitySamples",
            coeffs_.lookupOrDefault("nVelocitySamples", 20));

        bool writeElementPerf
        (
            coeffs_.lookupOrDefault("writeElementPerf", false)
        );
        dict.add("writePerf", writeElementPerf);

        if (debug)
        {
            Info<< "Creating actuatorCableLineElement: " << elemName << nl
                << "  position     : " << position << nl
                << "  P1           : " << P1 << nl
                << "  P2           : " << P2 << nl
                << "  spanDirection: " << spanDirection << nl
                << "  spanLength   : " << spanLength << nl
                << "  EA           : " << cableMat[0] << nl
                << "  T0           : " << cableMat[1] << nl
                << "  Cd           : " << Cd << endl;
        }

        actuatorCableLineElement* elem =
            new actuatorCableLineElement(elemName, dict, mesh_);
        elements_.set(i, elem);
        elements_[i].setVelocity(initialVelocity);
    }
}


void Foam::fv::actuatorCableLineSource::evaluateDeformation()
{
    scalar t = mesh_.time().value();
    if (t == lastMotionTime_) return;

    // ------------------------------------------------------------------
    // Assemble CableAnalysis input
    //
    // Node layout: one node at P1 of element 0, then one at the midpoint
    // and one at P2 for each element – identical scheme to the Bernoulli
    // version so the index arithmetic is unchanged.
    // Total nodes = 2*nElements + 1
    // Total elems = 2*nElements
    // ------------------------------------------------------------------

    const int nNodes = 2*nElements_ + 1;
    const int nElems = 2*nElements_;

    List<List<scalar>> CANodes(nNodes);
    List<List<int>>    CAElems(nElems);
    List<List<int>>    CARestraints(nNodes);
    List<List<scalar>> CAMats(nElems);      // [E_eff, A_eff] per FEA element
    List<List<scalar>> CALoads(nNodes);
    List<List<scalar>> CAPrescribed(nNodes);
    List<List<scalar>> CAPretension(nElems);

    List<scalar> sv(3, 0.0);
    List<int>    iv(3, 0);

    // --- First node (P1 of element 0) ---
    vector P1first = elements_[0].P1();
    sv[0] = P1first.x(); sv[1] = P1first.y(); sv[2] = P1first.z();
    CANodes[0]      = sv;
    CALoads[0]      = List<scalar>(3, 0.0);
    CAPrescribed[0] = List<scalar>(3, 0.0);

    forAll(elements_, i)
    {
        // --- Mid-node (element centroid) ---
        vector midPos = elements_[i].position();
        sv[0] = midPos.x(); sv[1] = midPos.y(); sv[2] = midPos.z();
        CANodes[2*i + 1] = sv;

        // --- End-node (P2) ---
        vector P2pos = elements_[i].P2();
        sv[0] = P2pos.x(); sv[1] = P2pos.y(); sv[2] = P2pos.z();
        CANodes[2*i + 2] = sv;

        // --- Loads at mid-node: incremental hydrodynamic drag ---
        // Use (current fluid force - previously applied structural force)
        // so the FEA sees only the increment, consistent with the Bernoulli
        // implementation.
        vector dF = elements_[i].cableForce() - elements_[i].structforce();
        sv[0] = dF.x(); sv[1] = dF.y(); sv[2] = dF.z();
        CALoads[2*i + 1] = sv;
        CALoads[2*i + 2] = List<scalar>(3, 0.0); // zero at end-node

        // Save current force as the new "structural" baseline
        elements_[i].setStructForce(elements_[i].cableForce());

        // --- Prescribed displacements: zero everywhere ---
        CAPrescribed[2*i + 1] = List<scalar>(3, 0.0);
        CAPrescribed[2*i + 2] = List<scalar>(3, 0.0);

        // --- Element connectivity ---
        List<int> con(2);
        con[0] = 2*i;     con[1] = 2*i + 1;
        CAElems[2*i]     = con;
        con[0] = 2*i + 1; con[1] = 2*i + 2;
        CAElems[2*i + 1] = con;

        // --- Material for sub-elements ---
        // CableAnalysis expects [E, A]; we store EA as a product, so we
        // pass E=cableEA (treating A=1 m2 as a normalisation convention).
        // Alternatively split: here we pass EA directly as E with A=1.
        List<scalar> mat(2);
        mat[0] = elements_[i].cableEA(); // E*A stored under E, A=1
        mat[1] = 1.0;
        CAMats[2*i]     = mat;
        CAMats[2*i + 1] = mat;

        // --- Pretension ---
        List<scalar> t0(1);
        t0[0] = elements_[i].cablePretension();
        CAPretension[2*i]     = t0;
        CAPretension[2*i + 1] = t0;

        // --- Restraints: zero on interior nodes ---
        CARestraints[2*i]     = List<int>(3, 0);
        CARestraints[2*i + 1] = List<int>(3, 0);
    }

    // Apply boundary restraints at first and last nodes
    CARestraints[0]          = elements_[0].cableRestraints();
    CARestraints[nNodes - 1] = elements_.last().cableRestraints();

    // ------------------------------------------------------------------
    // Solve
    // ------------------------------------------------------------------
    CableAnalysis CA
    (
        CANodes,
        CAElems,
        CARestraints,
        CAMats,
        CALoads,
        CAPrescribed,
        CAPretension
    );

    List<List<scalar>> deformations = CA.nodedispList();
    List<scalar>       tensions     = CA.tensionsList();

    // ------------------------------------------------------------------
    // Update element positions from deformed nodes
    // ------------------------------------------------------------------

    // Compute new absolute node positions
    List<vector> newNodePos(nNodes);
    forAll(newNodePos, k)
    {
        newNodePos[k] = vector
        (
            CANodes[k][0] + deformations[k][0],
            CANodes[k][1] + deformations[k][1],
            CANodes[k][2] + deformations[k][2]
        );
    }

    forAll(elements_, i)
    {
        elements_[i].setP1      (newNodePos[2*i]);
        elements_[i].setPosition(newNodePos[2*i + 1]);
        elements_[i].setP2      (newNodePos[2*i + 2]);

        // Accumulate total deformation at mid-node
        elements_[i].setDeformation
        (
            elements_[i].deformation()
          + vector
            (
                deformations[2*i + 1][0],
                deformations[2*i + 1][1],
                deformations[2*i + 1][2]
            )
        );

        // Update span geometry from deformed P1 -> P2
        vector span = elements_[i].P2() - elements_[i].P1();
        scalar len  = mag(span);
        if (len > VSMALL)
        {
            elements_[i].setSpanLength   (len);
            elements_[i].setSpanDirection(span / len);
        }

        // Store recovered tension (average of the two sub-elements
        // that share this actuator element)
        scalar T = 0.5*(tensions[2*i] + tensions[2*i + 1]);
        elements_[i].setTension(T);
    }

    lastMotionTime_ = t;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorCableLineSource::actuatorCableLineSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    actuatorLineSource(name, modelType, dict, mesh),
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
    )
{
    read(dict_);
    createInitialElements();
    evaluateDeformation();

    if (writePerf_)   createOutputFile();
    if (writeVTK_)    createOutputDir();

    if (forceField_.writeOpt() == IOobject::AUTO_WRITE)
        forceField_.write();

    if (endEffectsActive_)
        calcEndEffects();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::actuatorCableLineSource::~actuatorCableLineSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::fv::actuatorCableLineElement>&
Foam::fv::actuatorCableLineSource::elements()
{
    return elements_;
}


void Foam::fv::actuatorCableLineSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    evaluateDeformation();

    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ = vector::zero;

    forAll(elements_, i)
    {
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        elements_[i].calculateForce(U);
        elements_[i].addSup(eqn, forceField_);
        force_ += elements_[i].force();
    }

    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);

    eqn += forceField_;

    if (writePerf_ && Pstream::master())  writePerf();
    if (writeVTK_ && mesh_.time().outputTime() && Pstream::master()) writeVTK();
}


void Foam::fv::actuatorCableLineSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    evaluateDeformation();

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    forAll(elements_, i)
    {
        elements_[i].calculateForce(U);
        elements_[i].addTurbulence(eqn, fieldNames_[fieldI]);
    }
}


void Foam::fv::actuatorCableLineSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    evaluateDeformation();

    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ = vector::zero;

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    forAll(elements_, i)
    {
        elements_[i].calculateForce(U);
        elements_[i].addSup(rho, eqn, forceField_);
        force_ += elements_[i].force();
    }

    Info<< "Force on cable " << name_ << ": " << force_ << endl;

    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);

    eqn += forceField_;

    if (writePerf_ && Pstream::master())  writePerf();
    if (writeVTK_ && mesh_.time().outputTime() && Pstream::master()) writeVTK();
}


void Foam::fv::actuatorCableLineSource::writeVTK()
{
    fileName vtkFileName;
    std::ostringstream cfc;
    cfc << std::setw(12) << std::setfill('0') << vtkFileSequence_;
    vtkFileName = vtkDir_ + "/" + name_ + "_" + cfc.str() + ".vtk";

    vtkFilePtr_.reset(new OFstream(vtkFileName));

    vtkFilePtr_()
        << "# vtk DataFile Version 3.0" << nl
        << "actuator cable line " << name_ << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    // Points
    vtkFilePtr_() << "POINTS " << elements_.size() << " double" << nl;
    forAll(elements_, i)
    {
        vector p = elements_[i].position();
        vtkFilePtr_() << p[0] << " " << p[1] << " " << p[2] << nl;
    }

    // Connectivity
    vtkFilePtr_() << "LINES 1 " << elements_.size() + 1 << nl;
    vtkFilePtr_() << elements_.size() << " ";
    forAll(elements_, i) vtkFilePtr_() << i << " ";
    vtkFilePtr_() << nl << endl;

    // Point data
    vtkFilePtr_() << "POINT_DATA " << elements_.size() << nl;

    // Velocity
    vtkFilePtr_() << "VECTORS Velocity double" << nl;
    forAll(elements_, i)
    {
        vector v = elements_[i].velocity();
        vtkFilePtr_() << v[0] << " " << v[1] << " " << v[2] << nl;
    }
    vtkFilePtr_() << endl;

    // Hydrodynamic force
    vtkFilePtr_() << "VECTORS Force double" << nl;
    forAll(elements_, i)
    {
        vector f = elements_[i].force();
        vtkFilePtr_() << f[0] << " " << f[1] << " " << f[2] << nl;
    }
    vtkFilePtr_() << endl;

    // Deformation
    vtkFilePtr_() << "VECTORS Deformation double" << nl;
    forAll(elements_, i)
    {
        vector d = elements_[i].deformation();
        vtkFilePtr_() << d[0] << " " << d[1] << " " << d[2] << nl;
    }
    vtkFilePtr_() << endl;

    // Element tension (scalar)
    vtkFilePtr_() << "SCALARS Tension double 1" << nl
                  << "LOOKUP_TABLE default" << nl;
    forAll(elements_, i)
        vtkFilePtr_() << elements_[i].tension() << nl;
    vtkFilePtr_() << endl;

    // Net buoyancy force
    vtkFilePtr_() << "VECTORS BuoyancyForce double" << nl;
    forAll(elements_, i)
    {
        vector b = elements_[i].buoyancyForce();
        vtkFilePtr_() << b[0] << " " << b[1] << " " << b[2] << nl;
    }
    vtkFilePtr_() << endl;

    // Span direction
    vtkFilePtr_() << "VECTORS SpanDirection double" << nl;
    forAll(elements_, i)
    {
        vector s = elements_[i].spanDirection();
        vtkFilePtr_() << s[0] << " " << s[1] << " " << s[2] << nl;
    }
    vtkFilePtr_() << endl;

    vtkFileSequence_++;
}


// ************************************************************************* //
