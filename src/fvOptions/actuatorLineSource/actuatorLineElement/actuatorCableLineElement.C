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

#include "actuatorCableLineElement.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorCableLineElement, 0);
    defineRunTimeSelectionTable(actuatorCableLineElement, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::actuatorCableLineElement::read()
{
    dict_.lookup("position")         >> position_;
    dict_.lookup("P1")               >> P1_;
    dict_.lookup("P2")               >> P2_;
    dict_.lookup("CableEA")          >> cableEA_;
    dict_.lookup("CablePretension")  >> cablePretension_;
    dict_.lookup("CableDragCoeff")   >> cableDragCoeff_;
    dict_.lookup("CableRestraints")  >> cableRestraints_;
    dict_.lookup("CableDiameter")    >> cableDiameter_;
    dict_.lookup("CableDensity")     >> cableDensity_;
    dict_.lookup("CableFluidDensity") >> cableFluidDensity_;
    if (dict_.found("CableGravity"))  dict_.lookup("CableGravity") >> gravity_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorCableLineElement::actuatorCableLineElement
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    actuatorLineElement(name, dict, mesh),
    P1_(vector::zero),
    P2_(vector::zero),
    deformation_(vector::zero),
    cableForce_(vector::zero),
    tension_(0.0),
    cableEA_(dict.lookupOrDefault<scalar>("CableEA", 1.0e6)),
    cablePretension_(dict.lookupOrDefault<scalar>("CablePretension", 0.0)),
    cableDragCoeff_(dict.lookupOrDefault<scalar>("CableDragCoeff", 1.0)),
    cableRestraints_(3, 0),
    cableDiameter_(dict.lookupOrDefault<scalar>("CableDiameter", 0.05)),
    cableDensity_(dict.lookupOrDefault<scalar>("CableDensity", 7850.0)),
    cableFluidDensity_(dict.lookupOrDefault<scalar>("CableFluidDensity", 1025.0)),
    gravity_(dict.lookupOrDefault<vector>("CableGravity", vector(0, 0, -9.81))),
    buoyancyForce_(vector::zero),
    refSpanLength_(dict.lookupOrDefault<scalar>("spanLength", 1.0)),
    refSpanDirection_(dict.lookupOrDefault<vector>("spanDirection", vector(0,1,0))),
    positionInMesh_(true)
{
    // Override defaults with dictionary values if present
    if (dict.found("P1"))               dict.lookup("P1") >> P1_;
    if (dict.found("P2"))               dict.lookup("P2") >> P2_;
    if (dict.found("CableRestraints"))  dict.lookup("CableRestraints") >> cableRestraints_;

    // Attempt to read gravity from the mesh's 'g' field.
    // This is the standard OpenFOAM uniformDimensionedVectorField written
    // to the constant/ directory by solvers that use gravity.
    // If the field is absent (e.g., pure incompressible with no gravity)
    // the value set in the MIL / dictionary is kept.
    if (mesh.foundObject<uniformDimensionedVectorField>("g"))
    {
        gravity_ = mesh.lookupObject<uniformDimensionedVectorField>("g").value();
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::fv::actuatorCableLineElement::~actuatorCableLineElement()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// --- Access ---

const Foam::vector& Foam::fv::actuatorCableLineElement::P1() const
{ return P1_; }

const Foam::vector& Foam::fv::actuatorCableLineElement::P2() const
{ return P2_; }

const Foam::vector& Foam::fv::actuatorCableLineElement::deformation() const
{ return deformation_; }

const Foam::vector& Foam::fv::actuatorCableLineElement::cableForce() const
{ return cableForce_; }

const Foam::vector& Foam::fv::actuatorCableLineElement::dragForce() const
{ return dragForce_; }

Foam::scalar Foam::fv::actuatorCableLineElement::tension() const
{ return tension_; }

Foam::scalar Foam::fv::actuatorCableLineElement::cableEA() const
{ return cableEA_; }

Foam::scalar Foam::fv::actuatorCableLineElement::cablePretension() const
{ return cablePretension_; }

Foam::scalar Foam::fv::actuatorCableLineElement::cableDragCoeff() const
{ return cableDragCoeff_; }

const Foam::List<int>& Foam::fv::actuatorCableLineElement::cableRestraints() const
{ return cableRestraints_; }

Foam::scalar Foam::fv::actuatorCableLineElement::refSpanLength() const
{ return refSpanLength_; }

const Foam::vector& Foam::fv::actuatorCableLineElement::refSpanDirection() const
{ return refSpanDirection_; }

Foam::scalar Foam::fv::actuatorCableLineElement::cableDiameter() const
{ return cableDiameter_; }

Foam::scalar Foam::fv::actuatorCableLineElement::cableDensity() const
{ return cableDensity_; }

Foam::scalar Foam::fv::actuatorCableLineElement::cableFluidDensity() const
{ return cableFluidDensity_; }

const Foam::vector& Foam::fv::actuatorCableLineElement::gravity() const
{ return gravity_; }

const Foam::vector& Foam::fv::actuatorCableLineElement::buoyancyForce() const
{ return buoyancyForce_; }

const Foam::vector& Foam::fv::actuatorCableLineElement::spanDirection() const
{ return spanDirection_; }


// --- Manipulation ---

void Foam::fv::actuatorCableLineElement::setP1(vector P1)
{
    if (debug)
        Info<< "Changing P1 of " << name_ << " from " << P1_ << " to " << P1 << endl;
    P1_ = P1;
}

void Foam::fv::actuatorCableLineElement::setP2(vector P2)
{
    if (debug)
        Info<< "Changing P2 of " << name_ << " from " << P2_ << " to " << P2 << endl;
    P2_ = P2;
}

void Foam::fv::actuatorCableLineElement::setDeformation(vector d)
{
    if (debug)
        Info<< "Changing deformation of " << name_ << " from " << deformation_ << " to " << d << endl;
    deformation_ = d;
}

void Foam::fv::actuatorCableLineElement::setPosition(vector newPos)
{
    if (debug)
        Info<< "Changing position of " << name_ << " from " << position_ << " to " << newPos << endl;
    position_ = newPos;
}

void Foam::fv::actuatorCableLineElement::setSpanDirection(vector spanDir)
{
    spanDirection_ = spanDir;
}

void Foam::fv::actuatorCableLineElement::setSpanLength(scalar spanLength)
{
    spanLength_ = spanLength;
}

void Foam::fv::actuatorCableLineElement::setRefSpanLength(scalar refSpanLength)
{
    refSpanLength_ = refSpanLength;
}

void Foam::fv::actuatorCableLineElement::setRefSpanDirection(vector refSpanDir)
{
    refSpanDirection_ = refSpanDir;
}

void Foam::fv::actuatorCableLineElement::setTension(scalar t)
{
    tension_ = t;
}


void Foam::fv::actuatorCableLineElement::rotate
(
    vector rotationPoint,
    vector axis,
    scalar radians,
    bool rotateVelocity
)
{
    // Rodriguez rotation matrix (same as actuatorBernoulliLineElement)
    tensor RM;
    scalar c = Foam::cos(radians);
    scalar s = Foam::sin(radians);
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * c;
    RM.xy() = axis.x()*axis.y()*(1.0 - c) - axis.z()*s;
    RM.xz() = axis.x()*axis.z()*(1.0 - c) + axis.y()*s;
    RM.yx() = axis.x()*axis.y()*(1.0 - c) + axis.z()*s;
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * c;
    RM.yz() = axis.y()*axis.z()*(1.0 - c) - axis.x()*s;
    RM.zx() = axis.x()*axis.z()*(1.0 - c) - axis.y()*s;
    RM.zy() = axis.y()*axis.z()*(1.0 - c) + axis.x()*s;
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * c;

    // Rotate midpoint position
    vector pt = position_ - rotationPoint;
    position_ = (RM & pt) + rotationPoint;

    // Rotate P1
    pt = P1_ - rotationPoint;
    P1_ = (RM & pt) + rotationPoint;

    // Rotate P2
    pt = P2_ - rotationPoint;
    P2_ = (RM & pt) + rotationPoint;

    // Rotate span direction
    spanDirection_ = RM & spanDirection_;

    // Rotate structural data
    deformation_       = RM & deformation_;
    cableForce_        = RM & cableForce_;

    if (rotateVelocity)
        velocity_ = RM & velocity_;
}


void Foam::fv::actuatorCableLineElement::rotate
(
    vector rotationPoint,
    vector axis,
    scalar radians
)
{
    rotate(rotationPoint, axis, radians, true);
}


// --- Evaluation ---

void Foam::fv::actuatorCableLineElement::calculateForce
(
    const volVectorField& Uin
)
{
    // ------------------------------------------------------------------
    // Cable force calculation.
    //
    // IMPORTANT: we do NOT call the base class actuatorLineElement::
    // calculateForce() because it computes lift/drag from polar data and
    // performs angle-of-attack arithmetic that produces FPE (divide-by-
    // zero) when the relative velocity or chord length is near-zero.
    // All force physics for a cable are handled here explicitly.
    // ------------------------------------------------------------------

    // 1. Sample inflow velocity at element position.
    //    calculateInflowVelocity internally calls calcProjectionEpsilon,
    //    which issues a FatalError if the element position is not found in
    //    any mesh cell.  For a cable this is a legitimate condition: anchor
    //    nodes or elements that have deformed outside the CFD domain should
    //    contribute zero hydrodynamic force while still participating in the
    //    structural solve.  Check mesh membership first and return zero force
    //    if the element is outside the domain on all processors.
    {
        label cellI = mesh_.findCell(position_);
        label cellIGlobal = cellI;
        reduce(cellIGlobal, maxOp<label>());
        positionInMesh_ = (cellIGlobal >= 0);
        if (!positionInMesh_)
        {
            forceVector_   = vector::zero;
            cableForce_    = vector::zero;
            buoyancyForce_ = vector::zero;
            return;
        }
    }

    calculateInflowVelocity(Uin);

    relativeVelocity_ = inflowVelocity_ - velocity_;

    // ------------------------------------------------------------------
    // 2. Hydrodynamic drag (cross-flow, Morison-style)
    // ------------------------------------------------------------------

    // Unit span vector for velocity projection.
    // IMPORTANT: use the REFERENCE (undeformed) span direction, not the
    // current deformed spanDirection_.  After a buoyancy-driven catenary
    // deformation the deformed span is no longer aligned with the flow
    // reference frame; projecting relative velocity onto it introduces a
    // phantom normal-velocity component equal to the full inflow speed,
    // producing a drag force O(rho*Cd*A*U^2) at every step regardless of
    // the true cross-flow velocity.  That phantom drag then drives further
    // lateral deformation, rotating the span further, compounding each step
    // until the force diverges.  The reference span direction is the correct
    // basis for drag because drag is a fluid-dynamic quantity and the fluid
    // does not know or care that the cable has sagged.
    scalar magRefSpan = mag(refSpanDirection_);
    vector spanUnit = (magRefSpan > VSMALL)
                    ? refSpanDirection_ / magRefSpan
                    : vector(0, 0, 1);

    // Velocity component normal to span
    vector relNormal = relativeVelocity_
                     - spanUnit * (relativeVelocity_ & spanUnit);
    scalar magRelNormal = mag(relNormal);

    // Effective diameter: prefer explicit cableDiameter_, fall back to
    // chordLength_ (set from elementGeometry col [2])
    scalar diameter = (cableDiameter_ > VSMALL)
                    ? cableDiameter_
                    : max(chordLength_, VSMALL);

    // Use the REFERENCE (unstressed) span length for all area/volume
    // calculations.  Using the deformed span length creates a divergent
    // feedback loop: buoyancy and drag grow with stretch, driving further
    // displacement, further stretch, and so on.  The cable material volume
    // is conserved, so reference length is the physically correct choice.
    scalar spanLen = max(refSpanLength_, VSMALL);
    scalar projectedArea = diameter * spanLen;

    // Local fluid density: sample 'rho' field if present, else use
    // the user-supplied reference value
    scalar rhoFluid = max(cableFluidDensity_, VSMALL);
    if (Uin.mesh().foundObject<volScalarField>("rho"))
    {
        const volScalarField& rhoField =
            Uin.mesh().lookupObject<volScalarField>("rho");
        label cellI = Uin.mesh().findCell(position_);
        if (cellI >= 0)
            rhoFluid = max(rhoField[cellI], VSMALL);
    }

    // Drag: only non-zero when there is a meaningful relative velocity.
    dragForce_ = vector::zero;
    if (magRelNormal > VSMALL)
    {
        scalar drag = 0.5 * rhoFluid * cableDragCoeff_
                    * projectedArea * magSqr(relNormal);
        dragForce_ = drag * (relNormal / magRelNormal);
    }
	forceVector_=dragForce_;
    // ------------------------------------------------------------------
    // 3. Net buoyancy (buoyancy - self-weight)
    //    F_net = -(rho_f - rho_c) * g * V
    //    V = pi/4 * d^2 * L_ref
    //
    // Buoyancy is computed unconditionally — it does not depend on flow
    // velocity.  The previous guard (return zero when |U|~0) suppressed
    // buoyancy at startup, then it appeared as a step-change.
    // ------------------------------------------------------------------

    scalar pi = Foam::constant::mathematical::pi;
    scalar elemVolume = (pi / 4.0) * magSqr(diameter) * spanLen;

    buoyancyForce_ = -(rhoFluid - cableDensity_) * gravity_ * elemVolume;

    // ------------------------------------------------------------------
    // 4. Combine
    // ------------------------------------------------------------------

	// Feed back ONLY the hydrodynamic drag to the flow field; keep the full
	// external load (drag + buoyancy) for the structural solve.
	cableForce_  = dragForce_ + buoyancyForce_;

    
    if (debug)
    {
        Info<< "actuatorCableLineElement " << name_ << ":" << nl
            << "  position        : " << position_ << nl
            << "  spanUnit        : " << spanUnit << nl
            << "  inflowVelocity  : " << inflowVelocity_ << nl
            << "  relNormal       : " << relNormal << nl
            << "  rhoFluid        : " << rhoFluid << nl
            << "  diameter        : " << diameter << nl
            << "  spanLen         : " << spanLen << nl
            << "  drag force      : " << dragForce_ << nl
            << "  buoyancy force  : " << buoyancyForce_ << nl
            << "  total cableForce: " << cableForce_ << nl
            << "  tension         : " << tension_ << endl;
    }
}


// --- addSup overrides ---
// calculateForce() is now called by actuatorCableLineSource::addSup()
// BEFORE these per-element addSup() calls, so the force state is already
// current.  We only need to apply the force to the field here.
// applyForceField calls calcProjectionEpsilon, which issues a FatalError
// if the element position is not in any mesh cell.  positionInMesh_ was
// already set correctly by the preceding calculateForce() call, so we
// simply gate on it.


scalar Foam::fv::actuatorCableLineElement::calcProjectionEpsilon()
{
    // Lookup Gaussian coeffs from profileData dict if present
    dictionary GaussianCoeffs = profileData_.dict().subOrEmptyDict
    (
        "GaussianCoeffs"
    );
    scalar chordFactor = GaussianCoeffs.lookupOrDefault("chordFactor", 0.25);
    scalar dragFactor = GaussianCoeffs.lookupOrDefault("dragFactor", 1.0);
    scalar meshFactor = GaussianCoeffs.lookupOrDefault("meshFactor", 2.0);

    // Provide ideal epsilon target for lift based on chord length
    scalar epsilonLift = chordFactor*cableDiameter_;

    // Epsilon based on drag/momentum thickness
    scalar epsilonDrag = dragFactor*dragCoefficient_*cableDiameter_/2.0;

    // Threshold is based on lift or drag, whichever is larger
    scalar epsilonThreshold = Foam::max(epsilonLift, epsilonDrag);

    scalar epsilon = VGREAT;
    scalar epsilonMesh = VGREAT;
    const scalarField& V = mesh_.V();
    label posCellI = findCell(position_);
    if (posCellI >= 0)
    {
        // Projection width based on local cell size (from Troldborg (2008))
        epsilonMesh = 2.0*Foam::cbrt(V[posCellI]);
        epsilonMesh *= meshFactor; // Cell could have non-unity aspect ratio

        if (epsilonMesh > epsilonThreshold)
        {
            epsilon = epsilonMesh;
        }
        else
        {
            epsilon = epsilonThreshold;
        }
    }

    // Reduce epsilon over all processors
    reduce(epsilon, minOp<scalar>());

    // If epsilon is not reduced, position is not in the mesh
    if (not (epsilon < VGREAT))
    {
        // Raise fatal error since mesh size cannot be detected
        FatalErrorIn("void actuatorBernoulliLineElement::applyForceField()")
            << "Position "<< position_<<" of " << name_  << " not found in mesh"
            << abort(FatalError);
    }

    if (true)
    {
        reduce(epsilonMesh, minOp<scalar>());
        word epsilonMethod;
        if (epsilon == epsilonLift)
        {
            epsilonMethod = "lift-based";
        }
        else if (epsilon == epsilonDrag)
        {
            epsilonMethod = "drag-based";
        }
        else if (epsilon == epsilonMesh)
        {
            epsilonMethod = "mesh-based";
        }
    }

    return epsilon;
    
}


void Foam::fv::actuatorCableLineElement::applyForceField
(
    volVectorField& forceField
)
{
    // Calculate projection width
    scalar epsilon = calcProjectionEpsilon();
    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));

    // Apply force to the cells within the element's sphere of influence
    scalar sphereRadius = cableDiameter_ + projectionRadius;
    forAll(mesh_.cells(), cellI)
    {
        scalar dis = mag(mesh_.C()[cellI] - position_);
        if (dis <= sphereRadius)
        {
            scalar factor = Foam::exp(-Foam::sqr(dis/epsilon))
                          / (Foam::pow(epsilon, 3)
                          * Foam::pow(Foam::constant::mathematical::pi, 1.5));
            // forceField is opposite forceVector
            forceField[cellI] += -dragForce_*factor;
        }
    }

}

void Foam::fv::actuatorCableLineElement::addSup
(
    fvMatrix<vector>& eqn,
    volVectorField& forceField
)
{
    if (positionInMesh_)
        applyForceField(forceField);
}


void Foam::fv::actuatorCableLineElement::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    volVectorField& forceField
)
{
    if (positionInMesh_)
    {
        // Use a local scratch field for this element's contribution so that
        // multiplying by rho only affects this element's smeared force and
        // does not corrupt the contributions of elements already accumulated
        // in forceField.  This mirrors actuatorBernoulliLineElement::addSup.
        volVectorField forceFieldI
        (
            IOobject
            (
                "force." + name_,
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector
            (
                "zero",
                forceField.dimensions()/rho.dimensions(),
                vector::zero
            )
        );

        applyForceField(forceFieldI);

        // Scale forceVector_ by local density for force() bookkeeping
        multiplyForceRho(rho);

        // Multiply this element's field by density then accumulate
        forceFieldI *= rho;
        forceField  += forceFieldI;
    }
}


// ************************************************************************* //
