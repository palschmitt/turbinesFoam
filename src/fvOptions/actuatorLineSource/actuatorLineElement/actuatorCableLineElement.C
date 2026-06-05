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
    structforceVector_(vector::zero)
{
    // Override defaults with dictionary values if present
    if (dict.found("P1"))               dict.lookup("P1") >> P1_;
    if (dict.found("P2"))               dict.lookup("P2") >> P2_;
    if (dict.found("CableRestraints"))  dict.lookup("CableRestraints") >> cableRestraints_;
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
    // Interpolate flow velocity to element position
    calculateInflowVelocity(Uin);

    // Relative velocity (flow - element motion)
    relativeVelocity_ = inflowVelocity_ - velocity_;

    // Remove the component parallel to the span axis (cables only resist
    // cross-flow drag; axial drag is typically negligible for slender members)
    vector spanUnit = spanDirection_ / (mag(spanDirection_) + VSMALL);
    vector relNormal = relativeVelocity_
                     - spanUnit * (relativeVelocity_ & spanUnit);

    scalar magRelNormal = mag(relNormal);

    // Projected area of element normal to flow: diameter * spanLength
    // The element doesn't store a diameter directly; use chordLength_ as
    // the effective cross-sectional dimension.
    scalar projectedArea = chordLength_ * spanLength_;

    // Drag force per unit density (consistent with rest of actuator framework)
    scalar drag = 0.5 * cableDragCoeff_ * projectedArea * magSqr(relNormal);

    if (magRelNormal > VSMALL)
        forceVector_ = drag * (relNormal / magRelNormal);
    else
        forceVector_ = vector::zero;

    // Store as the structural load (used by the line source when assembling
    // the CableAnalysis load vector)
    cableForce_ = forceVector_;

    if (debug)
    {
        Info<< "actuatorCableLineElement " << name_ << ":" << nl
            << "  position        : " << position_ << nl
            << "  spanDirection   : " << spanDirection_ << nl
            << "  inflowVelocity  : " << inflowVelocity_ << nl
            << "  relativeVelocity: " << relativeVelocity_ << nl
            << "  relNormal       : " << relNormal << nl
            << "  drag force      : " << forceVector_ << nl
            << "  tension         : " << tension_ << endl;
    }
}


// ************************************************************************* //
