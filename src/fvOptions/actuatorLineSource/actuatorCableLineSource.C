/*---------------------------------------------------------------------------*\\
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
#include <cmath>
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
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
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
        // Entries required by actuatorLineElement::read() in the base class.
        // Cables have no pitch/chordMount concept; safe dummy values are used
        // so the base class initialises without error.
        dict.add("chordMount",     0.5);          // half-chord mount (unused)
        dict.add("pitch",          scalar(0.0));  // no pitch for cables
        dict.add("flowCurvature",  dictionary()); // empty sub-dict
        // Do NOT add 'dynamicStall' to the element dict.
        // The base class actuatorLineElement only constructs a dynamic stall
        // model when the key is present (mirrors actuatorFlexibleLineSource
        // behaviour).  All registered models require full polar data; cables
        // have none, so omitting the entry disables the model entirely.
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
        // profileName and profileData are required by the actuatorLineElement
        // base class read().  Cables do not use polar data so we inject the
        // dummy 'cable' profile that must be present in the fvOptions dict.
        word profileName = "cable";
        if (elementProfiles_.size() > 0)
            profileName = elementProfiles_
                [i * elementProfiles_.size() / nElements_];
        dict.add("profileName", profileName);
        if (profileData_.found(profileName))
            dict.add("profileData", profileData_.subDict(profileName));
        else
        {
            // Build a minimal flat polar so the element does not fatal
            dictionary dummyProfile;
            List<List<scalar>> dummyData(2);
            dummyData[0] = {-180.0, 0.0, 0.0};
            dummyData[1] = {  180.0, 0.0, 0.0};
            dummyProfile.add("data", dummyData);
            dict.add("profileData", dummyProfile);
        }
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
    // Store the undeformed reference node positions once.
    // These are passed to CableAnalysis every time step as the mesh nodes
    // so that L0_ (unstressed element length) stays constant regardless of
    // how much the cable has deformed.
    const int nNodesRef = 2*nElements_ + 1;
    refNodePos_.setSize(nNodesRef);
    refNodePos_[0] = elements_[0].P1();
    forAll(elements_, i)
    {
        refNodePos_[2*i + 1] = elements_[i].position();
        refNodePos_[2*i + 2] = elements_[i].P2();
        // Freeze the reference span length and direction on each element.
        // calculateForce() uses these (not the deformed values) for buoyancy
        // volume, projected drag area, and velocity projection — preventing
        // the divergent feedback loops caused by deformed-geometry coupling.
        scalar refLen = mag(elements_[i].P2() - elements_[i].P1());
        if (refLen > VSMALL)
        {
            elements_[i].setRefSpanLength(refLen);
            elements_[i].setRefSpanDirection
            (
                (elements_[i].P2() - elements_[i].P1()) / refLen
            );
        }
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
    // Pre-initialise all restraint rows to zero (free DOF).
    // Boundary nodes are set explicitly after the loop so they cannot be
    // clobbered by interior-node assignments inside the loop.
    for (int k = 0; k < nNodes; k++)
        CARestraints[k] = List<int>(3, 0);
    // Pre-initialise all load and prescribed rows to zero.
    // The element loop uses += to accumulate force contributions from
    // adjacent elements onto shared nodes, so a clean zero baseline is
    // required before the loop begins.
    for (int k = 0; k < nNodes; k++)
    {
        CALoads[k]      = List<scalar>(3, 0.0);
        CAPrescribed[k] = List<scalar>(3, 0.0);
    }
    // --- Assemble CableAnalysis node positions from undeformed reference ---
    // IMPORTANT: always pass the original undeformed geometry as CANodes so
    // that L0_ inside CableAnalysis remains constant across time steps.
    // The warm-start u0 carries the accumulated displacement from all previous
    // steps; Newton-Raphson then finds only the *correction* needed for the
    // current load state.
    for (int k = 0; k < nNodes; k++)
    {
        sv[0] = refNodePos_[k].x();
        sv[1] = refNodePos_[k].y();
        sv[2] = refNodePos_[k].z();
        CANodes[k] = sv;
    }
// Runtime-switchable debug output for structural state handoff.
const bool debugCableState =
    coeffs_.lookupOrDefault<bool>("debugCableState", false);

auto finiteScalar = [](const scalar s) -> bool
{
    return std::isfinite(s);
};

auto finiteVector = [&](const vector& v) -> bool
{
    return finiteScalar(v.x()) && finiteScalar(v.y()) && finiteScalar(v.z());
};

// Apply boundary restraints early so warm-start reconstruction can honour
// restrained DOFs.
CARestraints[0]          = elements_[0].cableRestraints();
CARestraints[nNodes - 1] = elements_.last().cableRestraints();

// ------------------------------------------------------------------
// Reconstruct previous nodal positions consistently from the CURRENT
// actuator-element geometry (P1, midpoint position, P2).
// ------------------------------------------------------------------
auto reconstructPrevNodePos = [&]() -> List<vector>
{
    List<vector> nodePos(nNodes, vector::zero);
    List<label>  nContrib(nNodes, 0);

    forAll(elements_, i)
    {
        nodePos[2*i]     += elements_[i].P1();
        nContrib[2*i]++;
        nodePos[2*i + 1] += elements_[i].position();
        nContrib[2*i + 1]++;
        nodePos[2*i + 2] += elements_[i].P2();
        nContrib[2*i + 2]++;
    }

    for (label k = 0; k < nNodes; ++k)
    {
        if (nContrib[k] > 0)
        {
            nodePos[k] /= scalar(nContrib[k]);
        }
        else
        {
            nodePos[k] = refNodePos_[k];
        }

        for (label d = 0; d < 3; ++d)
        {
            if (CARestraints[k][d])
            {
                nodePos[k][d] = refNodePos_[k][d] + CAPrescribed[k][d];
            }
        }

        if (!finiteVector(nodePos[k]))
        {
            FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
                << "Non-finite prevNodePos[" << k << "] = " << nodePos[k] << nl
                << "  refNodePos=" << refNodePos_[k] << nl
                << "  restraints=(" << CARestraints[k][0] << ' '
                << CARestraints[k][1] << ' ' << CARestraints[k][2] << ")"
                << abort(FatalError);
        }
    }

    return nodePos;
};

List<vector> prevNodePos = reconstructPrevNodePos();

if (debugCableState)
{
    Info<< name_ << ": evaluateDeformation() start" << nl
        << "  time   = " << t << nl
        << "  nNodes = " << nNodes << nl
        << "  nElems = " << nElems << endl;
    for (label k = 0; k < nNodes; ++k)
    {
        Info<< "  refNodePos[" << k << "]=" << refNodePos_[k]
            << "  prevNodePos[" << k << "]=" << prevNodePos[k] << endl;
    }
}

// --- Build warm-start displacement u0 directly from reconstructed nodal positions ---
List<List<scalar>> CAu0(nNodes);
for (label k = 0; k < nNodes; ++k)
{
    CAu0[k] = List<scalar>(3, 0.0);

    vector dPrev = prevNodePos[k] - refNodePos_[k];

    for (label d = 0; d < 3; ++d)
    {
        if (CARestraints[k][d])
        {
            CAu0[k][d] = CAPrescribed[k][d];
        }
        else
        {
            CAu0[k][d] = dPrev[d];
        }
    }

    if (debugCableState)
    {
        Info<< "  CAu0[" << k << "] = ("
            << CAu0[k][0] << ' ' << CAu0[k][1] << ' ' << CAu0[k][2] << ")"
            << endl;
    }
}

    // ------------------------------------------------------------------
    // Taut cable: clamp per-element load passed to FEA so that no
    // element receives a transverse force exceeding cableLoadCapFraction
    // times its axial stiffness EA/L.  This prevents a near-slack element
    // from being driven to a large displacement in the first FEA step after
    // the mode switches from kinematic to taut.
    // The capped forces are used only for FEA assembly; the original
    // cableForce() values are preserved for the CFD body force application.
    // ------------------------------------------------------------------
    scalar loadCapFraction =
        coeffs_.lookupOrDefault<scalar>("cableLoadCapFraction", 0.1);
    List<vector> cappedForce(nElements_);
    forAll(elements_, i)
    {
        scalar elemL = max(elements_[i].spanLength(), VSMALL);
        scalar maxF  = loadCapFraction * elements_[i].cableEA() / elemL;
        vector F     = elements_[i].cableForce();
        scalar Fmag  = mag(F);
        cappedForce[i] = (Fmag > maxF && Fmag > VSMALL)
                       ? F * (maxF / Fmag)
                       : F;
    }
    forAll(elements_, i)
    {
        // Use the load-capped force for FEA assembly (prevents near-slack
        // elements from being driven to large displacements).  The original
        // cableForce() is preserved for the CFD body force application.
        vector F = cappedForce[i];
        sv[0] = 0.5*F.x(); sv[1] = 0.5*F.y(); sv[2] = 0.5*F.z();
        CALoads[2*i + 1] = sv;
        // End-node accumulates contributions from adjacent elements; add here
        // (it was initialised to zero above; the next element's loop will also
        // add its 0.5*F contribution to node 2*i+2 as its own mid-node entry,
        // so no double-counting occurs for interior nodes).
        sv[0] = 0.5*F.x(); sv[1] = 0.5*F.y(); sv[2] = 0.5*F.z();
        List<scalar>& endLoad = CALoads[2*i + 2];
        endLoad[0] += sv[0];
        endLoad[1] += sv[1];
        endLoad[2] += sv[2];
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
    }
    // Boundary restraints were already applied before reconstructing CAu0.

    if (debugCableState)
    {
        Info<< name_ << ": assembled structural input" << nl;
        forAll(elements_, i)
        {
            Info<< "  element " << i
                << " cableForce=" << elements_[i].cableForce()
                << " cappedForce=" << cappedForce[i]
                << " tension=" << elements_[i].tension()
                << " spanLength=" << elements_[i].spanLength()
                << " refSpanLength=" << elements_[i].refSpanLength()
                << endl;
        }
        for (label k = 0; k < nNodes; ++k)
        {
            Info<< "  CALoads[" << k << "] = ("
                << CALoads[k][0] << " " << CALoads[k][1] << " " << CALoads[k][2] << ")"
                << "  CARestraints=("
                << CARestraints[k][0] << " " << CARestraints[k][1] << " " << CARestraints[k][2] << ")"
                << endl;
        }
    }
    // ------------------------------------------------------------------
    // Kinematic (low-tension) mode
    //
    // When cable tension is low the geometric stiffness T/L is near-zero,
    // Kff_ is nearly singular, and Newton-Raphson produces unbounded
    // displacements in response to any transverse load (buoyancy, drag).
    //
    // Trigger condition: the MEAN tension across all elements is below
    // slackTensionThreshold.  This is deliberately broader than the old
    // all-elements-AND-pretension test, which failed to trigger for buoyant
    // cables whose pretension was non-zero but whose solved tension was low.
    //
    // Kinematic solution: fit a parabolic catenary to the distributed
    // transverse load.  For a cable pinned at both ends with uniform
    // transverse load w [N/m] and horizontal span H:
    //
    //   y(s) = (w / (2*T_h)) * s * (H - s)       (sag parabola)
    //
    // where T_h is the horizontal tension component, estimated as
    //   T_h = w * H^2 / (8 * sag)
    // with sag capped at 10 % of the chord length to keep elements
    // inside the mesh.
    //
    // If the net transverse load is near zero (cable hangs straight),
    // the chord (straight line) is used, which is the correct solution.
    //
    // The catenary plane is defined by the chord direction and the net
    // transverse load direction.  Lateral drag is accounted for by
    // rotating the catenary plane to align with the resultant transverse
    // force vector, so the x/y divergence caused by drag-induced lateral
    // displacement feeding back into the projected area is suppressed.
    // ------------------------------------------------------------------
    scalar slackTol = coeffs_.lookupOrDefault<scalar>("slackTensionThreshold", 10.0);
    scalar meanTension = 0.0;
    forAll(elements_, i)
        meanTension += elements_[i].tension();
    meanTension /= max(scalar(nElements_), scalar(1));

    if (debugCableState)
    {
        Info<< name_ << ": meanTension=" << meanTension
            << "  slackTol=" << slackTol << endl;
    }

    if (!finiteScalar(meanTension))
    {
        FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
            << "Non-finite meanTension before slack/FEA branch: " << meanTension
            << abort(FatalError);
    }
    /*
    if (meanTension < slackTol)
    {
        // Anchor positions: use restrained end-nodes.
        // If neither end is restrained (unusual), just use refNodePos_ ends.
        vector anchor0 = refNodePos_[0];
        vector anchor1 = refNodePos_[nNodes - 1];

        if (!finiteVector(anchor0) || !finiteVector(anchor1))
        {
            FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
                << "Non-finite catenary anchor(s): anchor0=" << anchor0
                << " anchor1=" << anchor1 << abort(FatalError);
        }

        vector chord   = anchor1 - anchor0;
        scalar chordLen = mag(chord);

        if (debugCableState)
        {
            Info<< name_ << ": catenary mode anchors" << nl
                << "  anchor0 = " << anchor0 << nl
                << "  anchor1 = " << anchor1 << nl
                << "  chord   = " << chord << nl
                << "  chordLen= " << chordLen << endl;
        }

        if (!finiteVector(chord) || !finiteScalar(chordLen))
        {
            FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
                << "Non-finite chord data in catenary mode: chord=" << chord
                << " chordLen=" << chordLen << abort(FatalError);
        }
        if (chordLen < VSMALL)
        {
            // Degenerate: zero-length cable — leave positions unchanged
            lastMotionTime_ = t;
            return;
        }
        vector chordUnit = chord / chordLen;
        // Net transverse load per unit arc-length from all elements
        // (buoyancy dominates; drag may add a lateral component)
        vector totalTransverse = vector::zero;
        forAll(elements_, i)
        {
            vector F = elements_[i].cableForce();
            // Remove the axial (chord-parallel) component
            vector Faxial = (F & chordUnit) * chordUnit;
            totalTransverse += (F - Faxial);
        }
        // Transverse load per unit length (distributed, uniform approximation)
        scalar arcLen = max(totalLength_, VSMALL);
        vector wVec   = totalTransverse / arcLen; // [N/m], direction = sag plane
        scalar wMag   = mag(wVec);
        // Parabolic catenary: maximum sag d = w*H^2/(8*T_h)
        // We choose sag = min(maxSagFraction * chordLen, w*H^2/(8*T_est))
        // where T_est is a floor tension that keeps the cable inside the mesh.
        scalar maxSagFraction =
            coeffs_.lookupOrDefault<scalar>("catMaxSagFraction", 0.05);
        scalar maxSag = maxSagFraction * chordLen;
        // Sag direction: unit vector perpendicular to chord in the plane of wVec
        vector sagDir = vector::zero;
        if (wMag > VSMALL)
        {
            vector wUnit = wVec / wMag;
            // Remove chord component to get the true transverse direction
            sagDir = wUnit - (wUnit & chordUnit) * chordUnit;
            scalar sagDirMag = mag(sagDir);
            if (sagDirMag > VSMALL)
                sagDir /= sagDirMag;
            else
                sagDir = vector::zero; // load is purely axial — no sag
        }
        // Actual sag: capped at maxSag
        scalar sag = 0.0;
        if (wMag > VSMALL && mag(sagDir) > VSMALL)
            sag = min(wMag * sqr(chordLen) / max(8.0*slackTol, VSMALL), maxSag);
        // Distribute nodes along the catenary parabola
        forAll(elements_, i)
        {
            auto nodeFrac = [&](scalar idx) -> vector
            {
                scalar s = idx / scalar(nNodes - 1);     // 0 … 1
                vector pt = anchor0 + s * chord;         // chord point
                scalar sagVal = sag * 4.0 * s * (1.0 - s); // parabola
                return pt + sagVal * sagDir;
            };
            vector newP1  = nodeFrac(scalar(2*i));
            vector newMid = nodeFrac(scalar(2*i + 1));
            vector newP2  = nodeFrac(scalar(2*i + 2));

            if (!finiteVector(newP1) || !finiteVector(newMid) || !finiteVector(newP2))
            {
                FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
                    << "Non-finite catenary node positions for element " << i << nl
                    << "  newP1=" << newP1 << nl
                    << "  newMid=" << newMid << nl
                    << "  newP2=" << newP2 << abort(FatalError);
            }

            if (debugCableState)
            {
                Info<< "  catenary element " << i
                    << " newP1=" << newP1
                    << " newMid=" << newMid
                    << " newP2=" << newP2 << endl;
            }
            elements_[i].setP1      (newP1);
            elements_[i].setPosition(newMid);
            elements_[i].setP2      (newP2);
            // Deformation relative to reference
            elements_[i].setDeformation(newMid - refNodePos_[2*i + 1]);
            elements_[i].setTension(0.0);
            // Do NOT update spanDirection_ or spanLength_ from the deformed
            // catenary geometry.  calculateForce() uses refSpanDirection_ for
            // velocity projection (not spanDirection_), so spanDirection_ only
            // needs to be current for structural output / VTK.  Updating it
            // here from the sagged geometry would rotate it away from the
            // flow reference frame, causing phantom drag at the next step.
        }
        if (debug)
            Info<< name_ << ": kinematic (catenary) mode — "
                << "meanTension=" << meanTension
                << " < slackTol=" << slackTol
                << "  sag=" << sag << "  sagDir=" << sagDir << endl;
        lastMotionTime_ = t;
        return;
    }
    */
    // ------------------------------------------------------------------
    // Solve (taut cable — FEA Newton-Raphson)
    // ------------------------------------------------------------------
    int    caMaxIter = coeffs_.lookupOrDefault<int>   ("CableMaxIter",    200);
    scalar caTol     = coeffs_.lookupOrDefault<scalar>("CableTolerance", 1e-8);
    CableAnalysis CA
    (
        CANodes,
        CAElems,
        CARestraints,
        CAMats,
        CALoads,
        CAPrescribed,
        CAPretension,
        caMaxIter,
        caTol,
        CAu0
    );
    List<List<scalar>> deformations = CA.nodedispList();
    List<scalar>       tensions     = CA.tensionsList();
// ------------------------------------------------------------------
// Update element positions from deformed nodes
// ------------------------------------------------------------------
// deformations[] from CableAnalysis is the TOTAL displacement from the
// undeformed reference geometry (because CANodes = refNodePos_).
//
// Use the SAME previous nodal displacement field that was used to build
// CAu0, so timestep handoff is consistent.
scalar maxDispFraction =
    coeffs_.lookupOrDefault<scalar>("maxDispFraction", 0.5);

List<vector> newNodePos(nNodes);

for (label k = 0; k < nNodes; ++k)
{
    vector dTotal
    (
        deformations[k][0],
        deformations[k][1],
        deformations[k][2]
    );

    if (!finiteVector(dTotal))
    {
        FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
            << "Non-finite deformations[" << k << "] = " << dTotal
            << abort(FatalError);
    }

    vector dPrev = prevNodePos[k] - refNodePos_[k];
    vector dStep = dTotal - dPrev;

    label ownerElem = (k == 0 ? 0 : (k - 1)/2);
    if (ownerElem > nElements_ - 1)
    {
        ownerElem = nElements_ - 1;
    }

    scalar maxStep = maxDispFraction*elements_[ownerElem].spanLength();
    scalar stepMag = mag(dStep);

    if (!finiteScalar(maxStep) || !finiteScalar(stepMag) || !finiteVector(dPrev) || !finiteVector(dStep))
    {
        FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
            << "Non-finite nodal update state at node " << k << nl
            << "  dPrev   = " << dPrev << nl
            << "  dTotal  = " << dTotal << nl
            << "  dStep   = " << dStep << nl
            << "  maxStep = " << maxStep << nl
            << "  stepMag = " << stepMag << nl
            << "  ownerElem=" << ownerElem
            << abort(FatalError);
    }

    if (stepMag > maxStep && stepMag > VSMALL)
    {
        dStep *= maxStep/stepMag;
    }

    vector dNew = dPrev + dStep;

    for (label d = 0; d < 3; ++d)
    {
        if (CARestraints[k][d])
        {
            dNew[d] = CAPrescribed[k][d];
        }
    }

    newNodePos[k] = refNodePos_[k] + dNew;

    if (!finiteVector(newNodePos[k]))
    {
        FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
            << "Non-finite newNodePos[" << k << "] = " << newNodePos[k] << nl
            << "  refNodePos = " << refNodePos_[k] << nl
            << "  dPrev      = " << dPrev << nl
            << "  dStep      = " << dStep << nl
            << "  dNew       = " << dNew
            << abort(FatalError);
    }

    if (debugCableState)
    {
        Info<< "  node " << k
            << " dPrev=" << dPrev
            << " dTotal=" << dTotal
            << " dStep=" << dStep
            << " newNodePos=" << newNodePos[k] << endl;
    }
}

    forAll(elements_, i)
    {
        elements_[i].setP1      (newNodePos[2*i]);
        elements_[i].setPosition(newNodePos[2*i + 1]);
        elements_[i].setP2      (newNodePos[2*i + 2]);
        // Store the total displacement of this element's mid-node.
        elements_[i].setDeformation
        (
            vector
            (
                newNodePos[2*i + 1].x() - refNodePos_[2*i + 1].x(),
                newNodePos[2*i + 1].y() - refNodePos_[2*i + 1].y(),
                newNodePos[2*i + 1].z() - refNodePos_[2*i + 1].z()
            )
        );
        // Update span geometry from deformed P1 -> P2
        vector span = elements_[i].P2() - elements_[i].P1();
        scalar len  = mag(span);

        if (!finiteVector(span) || !finiteScalar(len))
        {
            FatalErrorIn("actuatorCableLineSource::evaluateDeformation()")
                << "Non-finite span for element " << i << nl
                << "  P1=" << elements_[i].P1() << nl
                << "  P2=" << elements_[i].P2() << nl
                << "  span=" << span << nl
                << "  len=" << len
                << abort(FatalError);
        }

        if (debugCableState)
        {
            Info<< "  updated element " << i
                << " P1=" << elements_[i].P1()
                << " mid=" << elements_[i].position()
                << " P2=" << elements_[i].P2()
                << " span=" << span
                << " len=" << len
                << " tension(avg)=" << 0.5*(tensions[2*i] + tensions[2*i + 1])
                << endl;
        }

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
    // Reset the motion-time guard so evaluateDeformation() always runs.
    // The guard exists to prevent duplicate solves when addSup is called
    // multiple times within the same time step (e.g. PIMPLE outer loops).
    // We control the call order here, so we reset it to -GREAT to guarantee
    // the structural solve runs with the freshly computed forces.
    lastMotionTime_ = -GREAT;
    // Step 1: compute hydrodynamic forces from current velocity field.
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    forAll(elements_, i)
        elements_[i].calculateForce(U);
    // Step 2: structural solve with current forces.
    evaluateDeformation();
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ = vector::zero;
    forAll(elements_, i)
    {
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
    // Forces must be current before the structural solve.
    lastMotionTime_ = -GREAT;
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    forAll(elements_, i)
        elements_[i].calculateForce(U);
    evaluateDeformation();
    forAll(elements_, i)
        elements_[i].addTurbulence(eqn, fieldNames_[fieldI]);
}
void Foam::fv::actuatorCableLineSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    lastMotionTime_ = -GREAT;
    // Step 1: compute hydrodynamic forces from current velocity field.
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    forAll(elements_, i)
        elements_[i].calculateForce(U);
    // Step 2: structural solve with current forces.
    evaluateDeformation();
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ = vector::zero;
    forAll(elements_, i)
    {
        elements_[i].addSup(rho, eqn, forceField_);
        // After the cable-element fix, force() is the CFD feedback force
        // (hydrodynamic drag only), not the full structural load.
        force_ += elements_[i].force();
    }
    scalar sumCableForceMag = 0.0;
    vector sumCableForce(vector::zero);
    vector sumBuoyancyForce(vector::zero);
    vector sumFeedbackForce(vector::zero);
    forAll(elements_, i)
    {
        sumCableForce += elements_[i].cableForce();
        sumCableForceMag += mag(elements_[i].cableForce());
        sumBuoyancyForce += elements_[i].buoyancyForce();
        sumFeedbackForce += elements_[i].force();
    }
    Info<< "Force on cable " << name_
        << ": appliedFeedback=" << force_
        << "  structuralInput(sum cableForce)=" << sumCableForce
        << "  sumBuoyancyForce=" << sumBuoyancyForce
        << "  sumFeedbackForce=" << sumFeedbackForce
        << "  sum|cableForce|=" << sumCableForceMag << endl;
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
