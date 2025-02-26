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

#include "actuatorBernoulliLineElement.H"
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
    defineTypeNameAndDebug(actuatorBernoulliLineElement, 0);
    defineRunTimeSelectionTable(actuatorBernoulliLineElement, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::actuatorBernoulliLineElement::read()
{
    // Parse dictionary
    dict_.lookup("position") >> position_;
    dict_.lookup("P1") >> P1_;
    dict_.lookup("P2") >> P2_;
    dict_.lookup("FEAmaterial") >> FEAmaterial_;
    dict_.lookup("FEAsection") >> FEAsects_;
    dict_.lookup("FEArestraints") >> FEArestraints_;

}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::fv::actuatorBernoulliLineElement::actuatorBernoulliLineElement(const word& name, const dictionary& dict, const fvMesh& mesh)
    : actuatorLineElement(name, dict, mesh)
    {
    // Initialize additional member variables specific to actuatorBernoulliLineElement
    P1_ = vector::zero;
    P2_ = vector::zero;
    FEAforce_ = vector::zero;
    FEAmoment_ = vector::zero;
    }

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorBernoulliLineElement::~actuatorBernoulliLineElement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




const Foam::vector& Foam::fv::actuatorBernoulliLineElement::P1()
{
    return P1_;
}
const Foam::vector& Foam::fv::actuatorBernoulliLineElement::P2()
{
    return P2_;
}
const Foam::vector& Foam::fv::actuatorBernoulliLineElement::FEAforce()
{
    return FEAforce_;
}
const Foam::vector& Foam::fv::actuatorBernoulliLineElement::FEAmoment()
{
    return FEAmoment_;
}
const Foam::scalar& Foam::fv::actuatorBernoulliLineElement::omega()
{
    return omega_;
}









Foam::scalar Foam::fv::actuatorBernoulliLineElement::inflowRefAngle()
{
    // Calculate inflow velocity angle in degrees (AFTAL Phi)
    scalar inflowVelAngleRad = VSMALL;
    if ((mag(chordRefDirection_)> VSMALL) and (mag(relativeVelocity_)>VSMALL) )
    {
    inflowVelAngleRad = acos
    (
        (-relativeVelocity_ & chordRefDirection_)
        / (mag(relativeVelocity_) * mag(chordRefDirection_))
    );
    }
    else
    {
        inflowVelAngleRad=0;
        }
    
    
    
    
    return radToDeg(inflowVelAngleRad);
}

void Foam::fv::actuatorBernoulliLineElement::calculateForce
(
    const volVectorField& Uin
)
{
    scalar pi = Foam::constant::mathematical::pi;

    // Calculate vector normal to chord--span plane
    planformNormal_ = -chordDirection_ ^ spanDirection_;
    planformNormal_ /= mag(planformNormal_);

    if (debug)
    {
        Info<< "Calculating force contribution from actuatorBernoulliLineElement "
            << name_ << endl;
        Info<< "    position: " << position_ << endl;
        Info<< "    chordDirection: " << chordDirection_ << endl;
        Info<< "    spanDirection: " << spanDirection_ << endl;
        Info<< "    elementVelocity: " << velocity_ << endl;
        Info<< "    planformNormal: " << planformNormal_ << endl;
    }

    // Find local flow velocity by interpolating to element location
    calculateInflowVelocity(Uin);

    // Subtract spanwise component of inflow velocity
    vector spanwiseVelocity = spanDirection_
                            * (inflowVelocity_ & spanDirection_)
                            / magSqr(spanDirection_);
    inflowVelocity_ -= spanwiseVelocity;

    // Calculate relative velocity and Reynolds number
    relativeVelocity_ = inflowVelocity_ - velocity_;
    Re_ = mag(relativeVelocity_)*chordLength_/nu_;
    
    //Segfaults for zero velocity
    if (mag(Re_)>VSMALL)
        {
        // Calculate angle of attack (radians)
        scalar angleOfAttackRad = asin((planformNormal_ & relativeVelocity_)
        / (mag(planformNormal_)
        *  mag(relativeVelocity_)));
        scalar angleOfAttackUncorrected = radToDeg(angleOfAttackRad);
        relativeVelocityGeom_ = freeStreamVelocity_ - velocity_;
        angleOfAttackGeom_ = asin((planformNormal_ & relativeVelocityGeom_)
        / (mag(planformNormal_)*mag(relativeVelocityGeom_)));
        angleOfAttackGeom_ *= 180.0/pi;
        
        // Apply flow curvature correction to angle of attack
        if (flowCurvatureActive_)
        {
        correctFlowCurvature(angleOfAttackRad);
        }
        
        // Calculate angle of attack in degrees
        angleOfAttack_ = radToDeg(angleOfAttackRad);
        
        // Update Reynolds number of profile data
        profileData_.updateRe(Re_);
        
        // Lookup lift and drag coefficients
        lookupCoefficients();
        
        if (debug)
        {
        Info<< "    inflowVelocity: " << inflowVelocity_ << endl;
        Info<< "    relativeVelocity: " << relativeVelocity_ << endl;
        Info<< "    Reynolds number: " << Re_ << endl;
        Info<< "    Geometric angle of attack (degrees): "
        << angleOfAttackGeom_ << endl;
        Info<< "    Angle of attack (uncorrected, degrees): "
        << angleOfAttackUncorrected << endl;
        Info<< "    Angle of attack (corrected, degrees): "
        << angleOfAttack_ << endl;
        }
        
        // Correct coefficients with dynamic stall model
        if (dynamicStallActive_)
        {
        dynamicStall_->correct
        (
        mag(relativeVelocity_),
        angleOfAttack_,
        liftCoefficient_,
        dragCoefficient_,
        momentCoefficient_
        );
        }
        
        // Correct for added mass effects
        if (addedMassActive_)
        {
        addedMass_.correct
        (
        liftCoefficient_,
        dragCoefficient_,
        momentCoefficient_,
        degToRad(angleOfAttack_),
        mag(chordDirection_ & relativeVelocity_),
        mag(planformNormal_ & relativeVelocity_)
        );
        }
        
        // Apply end effect correction factor to lift coefficient
        liftCoefficient_ *= endEffectFactor_;
        
        // Calculate force per unit density
        scalar area = chordLength_ * spanLength_;
        scalar magSqrU = magSqr(relativeVelocity_);
        scalar lift = 0.5*area*liftCoefficient_*magSqrU;
        scalar drag = 0.5*area*dragCoefficient_*magSqrU;
        vector liftDirection = relativeVelocity_ ^ spanDirection_;
        liftDirection /= mag(liftDirection);
        vector dragDirection = relativeVelocity_/mag(relativeVelocity_);
        forceVector_ = lift*liftDirection + drag*dragDirection;
        
        if (debug)
        {
        Info<< "    liftDirection: " << liftDirection << endl;
        Info<< "    dragDirection: " << dragDirection << endl;
        Info<< "    force (per unit density): " << forceVector_ << endl;
        }
        }
    else
    {
        relativeVelocity_=vector(VSMALL, VSMALL, VSMALL);
        Re_=VSMALL;
    }
}

void Foam::fv::actuatorBernoulliLineElement::setPosition(vector NewPosition)
{
	position_=NewPosition;
	}
void Foam::fv::actuatorBernoulliLineElement::setSpanLength(scalar SP)
			{
				spanLength_=SP;
				}
void Foam::fv::actuatorBernoulliLineElement::setSpanDirection(vector SD)
			{
				spanDirection_=SD;
				}

  
void Foam::fv::actuatorBernoulliLineElement::rotate
(
    vector rotationPoint,
    vector axis,
    scalar radians,
    bool rotateVelocity=true
)
{
    // Declare and define the rotation matrix (from SOWFA)
    //Why not use quaternions or rodrigues formula?
    tensor RM;
    scalar angle = radians;
    RM.xx() = Foam::sqr(axis.x())
            + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle);
    RM.xz() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y())
            + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z())
            + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    if (debug)
    {
        Info<< "Rotating actuatorBernoulliLineElement: " << name_ << endl;
        Info<< "Rotation point: " << rotationPoint << endl;
        Info<< "Rotation axis: " << axis << endl;
        Info<< "Rotation angle (radians): " << radians << endl;
        Info<< "Rotation matrix:" << endl << RM << endl;
        Info<< "Initial position: " << position_ << endl;
        Info<< "Initial chordDirection: " << chordDirection_ << endl;
        Info<< "Initial spanDirection: " << spanDirection_ << endl;
        Info<< "Initial velocity: " << velocity_ << endl;
        Info<< "Initial structforceVector: " << structforceVector_<< endl;
    }
        
    
    
      // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated. Ugly, write a proper function!
    vector point = position_;  
    
    point -= rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation
    // point
    point += rotationPoint;

    // Set the position of the element
    position_ = point;
    
    //Rotate P1
    point = P1_;  
    
    point -= rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation
    // point
    point += rotationPoint;

    // Set the position of the element
    P1_ = point;
    //Rotate P2
    point = P2_;  
    
    point -= rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation
    // point
    point += rotationPoint;

    // Set the position of the element
    P2_ = point;

    // Rotate the span and chord vectors of the element
    chordDirection_ = RM & chordDirection_;
    spanDirection_ = RM & spanDirection_;
    deformation_ = RM & deformation_;
    structforceVector_ = RM & structforceVector_;
    structmomentVector_= RM & structmomentVector_;    
    // Rotate the element's velocity vector if specified
    if (rotateVelocity)
    {
        velocity_ = RM & velocity_;
        chordRefDirection_ = RM & chordRefDirection_;

    }

    if (debug)
    {
        Info<< "Final position: " << position_ << endl;
        Info<< "Final chordDirection: " << chordDirection_ << endl;
        Info<< "Final chordRefDirection: " << chordRefDirection_ << endl;
        Info<< "Final spanDirection: " << spanDirection_ << endl;
        Info<< "Final velocity: " << velocity_ << endl << endl;
        Info<< "Final StructForce: " << structforceVector_ << endl << endl;
        Info<< "Final StructMoment: " << structmomentVector_ << endl << endl;
    }

    
}




void Foam::fv::actuatorBernoulliLineElement::setP1(vector P1)
{
    if (debug)
    {
        Info<< "Changing P1 of " << name_ << " from "
            << P1_ << " to " << P1 << endl << endl;
    }
    P1_ = P1;
}
void Foam::fv::actuatorBernoulliLineElement::setP2(vector P2)
{
    if (debug)
    {
        Info<< "Changing P1 of " << name_ << " from "
            << P2_ << " to " << P2 << endl << endl;
    }
    P2_ = P2;
}
void Foam::fv::actuatorBernoulliLineElement::setFEAforce(vector f)
{
    if (debug)
    {
        Info<< "Changing FEAforce of " << name_ << " from "
            << FEAforce_ << " to " << f << endl << endl;
    }
    FEAforce_ = f;
}
void Foam::fv::actuatorBernoulliLineElement::setFEAmoment(vector f)
{
    if (debug)
    {
        Info<< "Changing FEAforce of " << name_ << " from "
            << FEAmoment_<< " to " << f << endl << endl;
    }
    FEAmoment_ = f;
}
void Foam::fv::actuatorBernoulliLineElement::setStructForce(vector Structforce)
{
    if (debug)
    {
        Info<< "Changing Structforce of " << name_ << " from "
            << structforceVector_ << " to " << Structforce << endl << endl;
    }
    structforceVector_ = Structforce;
}
void Foam::fv::actuatorBernoulliLineElement::setStructMoment(vector Structmoment)
{
    if (debug)
    {
        Info<< "Changing Structmoment of " << name_ << " from "
            << structmomentVector_ << " to " << Structmoment << endl << endl;
    }
    structmomentVector_ = Structmoment;
}
void Foam::fv::actuatorBernoulliLineElement::setDeformation(vector Deformation)
{
    if (debug)
    {
        Info<< "Changing Deformation of " << name_ << " from "
            << deformation_ << " to " << Deformation << endl << endl;
    }
    deformation_ = Deformation;
}


void Foam::fv::actuatorBernoulliLineElement::setFEASects(List<scalar> s)
{
    if (debug)
    {
        Info<< "Changing Section Data of " << name_ << " from "
            << FEAsects_ << " to " << s << endl << endl;
    }
    FEAsects_ = s;
}
void Foam::fv::actuatorBernoulliLineElement::setFEAMaterial(List<scalar> s)
{
    if (debug)
    {
        Info<< "Changing Section Data of " << name_ << " from "
            << FEAmaterial_ << " to " << s << endl << endl;
    }
    FEAmaterial_ = s;
}
void Foam::fv::actuatorBernoulliLineElement::setFEARestraint(List<int> s)
{
    if (debug)
    {
        Info<< "Changing Restraint of " << name_ << " from "
            << FEArestraints_ << " to " << s << endl << endl;
    }
    FEArestraints_ = s;
}

const Foam::vector& Foam::fv::actuatorBernoulliLineElement::structforce()
{
    return structforceVector_;
}
const Foam::vector& Foam::fv::actuatorBernoulliLineElement::structmoment()
{
    return structmomentVector_;
}
const Foam::vector& Foam::fv::actuatorBernoulliLineElement::deformation()
{
    return deformation_;
}
const Foam::List<int>& Foam::fv::actuatorBernoulliLineElement::FEArestraints()
{
    return FEArestraints_;
}
const Foam::List<scalar>& Foam::fv::actuatorBernoulliLineElement::FEAsects()
{
    return FEAsects_;
}
const Foam::List<scalar>& Foam::fv::actuatorBernoulliLineElement::FEAmaterial()
{
    return FEAmaterial_;
}





const Foam::vector Foam::fv::actuatorBernoulliLineElement::pitchingMoment()
{
    //If chordmount differs from 0.25 this creates a pitching moment even if CM=0
	//return (forceVector_^(((chordMount_-0.25)*chordLength_)*chordDirection_))
      return      0.5*chordLength_*chordLength_*spanLength_
                          * momentCoefficient_*magSqr(relativeVelocity_)
                          * spanDirection_;
                          
	
	}











// ************************************************************************* //
