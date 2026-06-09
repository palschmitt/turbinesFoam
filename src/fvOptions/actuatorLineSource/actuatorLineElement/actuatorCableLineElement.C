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
    structforceVector_(vector::zero)
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

const Foam::vector& Foam::fv::actuatorCableLineElement::structforce() const
{ return structforceVector_; }

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

void Foam::fv::actuatorCableLineElement::setStructForce(vector sf)
{
    if (debug)
        Info<< "Changing structforce of " << name_ << " from " << structforceVector_ << " to " << sf << endl;
    structforceVector_ = sf;
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
    structforceVector_ = RM & structforceVector_;
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
    // 1. Hydrodynamic drag (cross-flow, Morison-style)
    // ------------------------------------------------------------------

    calculateInflowVelocity(Uin);
    relativeVelocity_ = inflowVelocity_ - velocity_;

    // Component of relative velocity normal to the span axis
    vector spanUnit = spanDirection_ / (mag(spanDirection_) + VSMALL);
    vector relNormal = relativeVelocity_
                     - spanUnit * (relativeVelocity_ & spanUnit);
    scalar magRelNormal = mag(relNormal);

    // Projected area: use cable diameter (d) * span length.
    // cableDiameter_ is the user-supplied outer diameter; chordLength_ is
    // kept as a fallback for compatibility with the base class.
    scalar diameter = (cableDiameter_ > VSMALL) ? cableDiameter_ : chordLength_;
    scalar projectedArea = diameter * spanLength_;

    // Sample local fluid density if a 'rho' field is present (compressible
    // or multiphase solvers), otherwise use the reference value.
    scalar rhoFluid = cableFluidDensity_;
    if (Uin.mesh().foundObject<volScalarField>("rho"))
    {
        const volScalarField& rhoField =
            Uin.mesh().lookupObject<volScalarField>("rho");
        // Interpolate to element position using the same cell as the
        // inflow velocity interpolation.
        label cellI = Uin.mesh().findCell(position_);
        if (cellI >= 0)
            rhoFluid = rhoField[cellI];
    }

    // Drag: F_d = 0.5 * rho_f * Cd * A_proj * |U_n|^2 * U_n/|U_n|
    vector dragForce = vector::zero;
    if (magRelNormal > VSMALL)
    {
        scalar drag = 0.5 * rhoFluid * cableDragCoeff_
                    * projectedArea * magSqr(relNormal);
        dragForce = drag * (relNormal / magRelNormal);
    }

    // ------------------------------------------------------------------
    // 2. Net buoyancy force (buoyancy - self-weight)
    //
    //   Displaced volume  V = pi/4 * d^2 * L_span
    //   Buoyancy          F_b =  rho_f * |g| * V   (acts opposite to g)
    //   Self-weight       W   =  rho_c * |g| * V   (acts along g)
    //   Net body force    F_net = (rho_f - rho_c) * g_unit * |g| * V
    //                          = (rho_f - rho_c) * g * V
    //
    // Note: gravity_ stores the vector g (e.g. (0 0 -9.81)).  The
    // buoyancy force is -rho_f*g*V (upward), self-weight is rho_c*g*V
    // (downward), so net = -(rho_f - rho_c)*g*V.
    // ------------------------------------------------------------------

    scalar pi = Foam::constant::mathematical::pi;
    scalar elemVolume = (pi / 4.0) * magSqr(diameter) * spanLength_;

    // Net upward body force on the structural node:
    //   positive when cable is lighter than fluid (floats)
    //   negative when cable is heavier than fluid (sinks)
    buoyancyForce_ = -(rhoFluid - cableDensity_) * gravity_ * elemVolume;

    // ------------------------------------------------------------------
    // 3. Combine and store
    // ------------------------------------------------------------------

    forceVector_ = dragForce + buoyancyForce_;
    cableForce_  = forceVector_;

    if (debug)
    {
        Info<< "actuatorCableLineElement " << name_ << ":" << nl
            << "  position        : " << position_ << nl
            << "  spanDirection   : " << spanDirection_ << nl
            << "  inflowVelocity  : " << inflowVelocity_ << nl
            << "  relNormal       : " << relNormal << nl
            << "  rhoFluid        : " << rhoFluid << nl
            << "  drag force      : " << dragForce << nl
            << "  buoyancy force  : " << buoyancyForce_ << nl
            << "  total cableForce: " << cableForce_ << nl
            << "  tension         : " << tension_ << endl;
    }
}


// ************************************************************************* //
