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

Class
    Foam::fv::actuatorBernoulliLineElement

Description
    An object, e.g. an airfoil, that imparts a force on the flow field
    calculated via lift and drag coefficients
    Should be derived from original but too lazy today for nice code

SourceFiles
    actuatorBernoulliLineElement.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorBernoulliLineElement_H
#define actuatorBernoulliLineElement_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicStallModel.H"
#include "interpolationCellPoint.H"
#include "profileData.H"
#include "addedMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class actuatorBernoulliLineElement Declaration
\*---------------------------------------------------------------------------*/

class actuatorBernoulliLineElement
{

protected:

    // Protected data

        //- Coefficients dictionary
        const dictionary dict_;

        //- Name of element
        const word name_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Mesh bounding box
        boundBox meshBoundBox_;

        //- Chord direction -- pointing from trailing to leading edge
        vector chordDirection_;

        //- Chord length
        scalar chordLength_;

        //- Span direction -- direction about which element will pitch
        //  Can be obtained by sweeping chord direction into planform normal
        // Unit Length?
        vector spanDirection_;

        //- Element span length, for calculating force
        scalar spanLength_;

        //- Element planform normal vector -- direction when the opposite of
        //  the chord direction is swept into the span direction
        vector planformNormal_;

        //- Element chord reference direction
        vector chordRefDirection_;

        //- Location of element's half chord
        vector position_;
        //- Location of element's start point
        vector P1_;
        //- Location of element's end point
        vector P2_;

        //- Velocity of the element
        vector velocity_;

        //- Free stream velocity
        vector freeStreamVelocity_;

        //- Free stream direction
        vector freeStreamDirection_;

        //- Fluid force on element -- vector that is perpendicular to
        //  spanDirection,
        vector forceVector_;
        
        //- Structural force on element -- vector that is perpendicular to
        //  spanDirection,
        vector structforceVector_;
        //- Structural moment on element -- 
        vector structmomentVector_;
        //Total applied displacement from FEA
        vector displacement_;
        //Internal (shear) force from FEA
        vector FEAforce_;
        //Internal (bending) moment from FEA
        vector FEAmoment_;
        //- Structural stiffness of element 
        List<scalar> FEAmaterial_;// E  Poisson
        List<scalar> FEAsects_;//Section data A        Iz       Iy          J        alpha
        List<int> FEArestraints_;//restrained X Y Z Rx Ry Rz 0 or 1

        //- Reference density for incompressible case
        scalar rhoRef_;

        //- Inflow velocity
        vector inflowVelocity_;

        //- Relative flow velocity
        vector relativeVelocity_;

        //- Geometric relative flow velcity
        vector relativeVelocityGeom_;

        //- Angle of attack (degrees)
        scalar angleOfAttack_;

        //- Geometric angle of attack (degrees, no flow curvature correction)
        scalar angleOfAttackGeom_;

        //- Lift coefficient
        scalar liftCoefficient_;

        //- Drag coefficient;
        scalar dragCoefficient_;

        //- 1/4 chord moment coefficient
        scalar momentCoefficient_;

        //- Profile name
        word profileName_;

        //- Profile data (encapsulates coefficients)
        profileData profileData_;

        //- Dynamic stall model
        autoPtr<dynamicStallModel> dynamicStall_;

        //- Switch for applying dynamic stall model
        bool dynamicStallActive_;

        //- Kinematic viscosity (for calculating Re)
        scalar nu_;

        //- Chord Reynolds number based on relative velocity
        scalar Re_;

        //- Angular velocity for flow curvature correction
        scalar omega_;

        //- Nondimensional chord mount distance
        scalar chordMount_;

        //- Switch for applying flow curvature correction
        bool flowCurvatureActive_;

        //- Flow curvature model name
        word flowCurvatureModelName_;

        //- Leading edge velocity vector
        vector velocityLE_;

        //- Trailing edge velocity vector
        vector velocityTE_;

        //- Switch for writing performance
        bool writePerf_;

        //- Output file stream
        OFstream* outputFile_;

        //- Nondimensional distance from actuator line root
        scalar rootDistance_;

        //- End effect correction factor [0, 1]
        scalar endEffectFactor_;

        //- Switch for added mass correction
        bool addedMassActive_;

        //- Added mass model
        addedMassModel addedMass_;

        //- Radius of circle used to sample velocities
        scalar velocitySampleRadius_;

        //- Number of elements used to sample velocities
        label nVelocitySamples_;


    // Protected Member Functions

        //- Rotate a vector
        void rotateVector
        (
            vector& vectorToRotate,
            vector rotationPoint,
            vector axis,
            scalar radians
        );

        //- Find cell containing location
        label findCell(const point& location);

        //- Lookup force coefficients
        void lookupCoefficients();

        //- Calculate projection width epsilon
        scalar calcProjectionEpsilon();

        //- Correct for flow curvatue
        void correctFlowCurvature(scalar& angleOfAttackRad);

        //- Multiply force vector by local density
        void multiplyForceRho(const volScalarField& rho);

        //- Apply force field based on force vector
        void applyForceField(volVectorField& forceField);

        //- Get inflow velocity
        void calculateInflowVelocity(const volVectorField& Uin);

        //- Create the performance output file
        virtual void createOutputFile();

        //- Write performance to CSV
        void writePerf();


public:

    //- Runtime type information
    TypeName("actuatorBernoulliLineElement");


        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            actuatorBernoulliLineElement,
            dictionary,
            (
                const dictionary& dict,
                const word& modelName
            ),
            (dict, modelName)
        );


    // Selectors

        //- Return a reference
        static autoPtr<actuatorBernoulliLineElement> New(const dictionary& dict);


    //- Constructor
    actuatorBernoulliLineElement
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~actuatorBernoulliLineElement();


    // Member functions

        // Access

            //- Return const access to the element name
            const word& name() const;

            //- Return the element chord length
            const scalar& chordLength() const;

            //- Return the element span length
            const scalar& spanLength();

            //- Return the element location
            const vector& position();
                        
            //- Return the element velocity
            const vector& velocity();
            //- Return the element P1
            const vector& P1();
            //- Return the element P2
            const vector& P2();
            //- Return the FEAforce 
            const vector& FEAforce();
            //- Return the FEAmoment 
            const vector& FEAmoment();
            
            //- Return the element angular velocity
            const scalar& omega();

            //- Return force vector from element
            const vector& force();
            
            //- Return structural force vector from element
            const vector& structforce();
            //- Return structural moment vector from element
            const vector& structmoment();
            //- Return last applied deformation
            const vector& displacement();
            //- Return stiffness vector from element
            const List<scalar>& FEAsects();
            const List<scalar>& FEAmaterial();
            const List<int>& FEArestraints();

            //- Return relative velocity
            const vector& relativeVelocity();

            //- Return geometric relative velocity
            const vector& relativeVelocityGeom();

            //- Return element angle of attack in degrees
            const scalar& angleOfAttack();

            //- Return element geometric angle of attack in degrees
            const scalar& angleOfAttackGeom();

            //- Return lift coefficient
            const scalar& liftCoefficient();

            //- Return drag coefficient
            const scalar& dragCoefficient();

            //- Return moment coefficient
            const scalar& momentCoefficient();
            
            
            const vector pitchingMoment();

            //- Return nondimensional distance from actuator line root
            const scalar& rootDistance();

            //- Return tangential reference coefficient
            scalar tangentialRefCoefficient();
            
            //- Return normal reference coefficient
            scalar normalRefCoefficient();

            //- Return tangential reference force (per unit span)
            scalar tangentialRefForce();
            
            //- Return normal reference force (per unit span)
            scalar normalRefForce();

            //- Return the inflow velocity angle w.r.t. the chord reference
            //  direction (degrees)
            //  E.g., for an axial-flow turbine this would be the angle of
            //  attack plus twist
            scalar inflowRefAngle();

			const vector&  chordDirection();
			
			const vector&  planformNormal();
			
			const vector&  spanDirection();
        // Manipulation
			void setSpanLength(scalar);
			void setSpanDirection(vector);
			void deltaChordDirection(scalar);
			
			void setPosition(vector NewPosition);

            //- Translate the element's location
            void translate(vector translationVector);

            //- Rotate the element about some axis
            void rotate
            (
                vector rotationPoint,
                vector axis,
                scalar radians,
                bool rotateVelocity
            );
            void rotate
            (
                vector rotationPoint,
                vector axis,
                scalar radians
            );

		scalar foeppl
		(
			scalar x,
			scalar pos,
			scalar exp
		);

            //- Pitch the element about some fraction of the chord
            void pitch(scalar radians, scalar chordFraction=0.25);

            //- Set velocity of element
            void setVelocity(vector velocity);
            //- Set velocity of element
            void setP1(vector P1);
            //- Set velocity of element
            void setP2(vector P2);
            //- Set FAEforce of element
            void setFEAforce(vector f);
            //- Set FAEmoment of element
            void setFEAmoment(vector m);
            
            //Set structural force
            void setStructForce(vector Structforce);
            void setStructMoment(vector Structmoment);
            
            void setDisplacement(vector Displacement);
            
            void setFEASects(List<scalar> FEAsects);            
            void setFEAMaterial(List<scalar> FEAmaterial);            
            void setFEARestraint(List<int> FEArestraint);            

            //- Set element speed (velocity magnitude)
            void setSpeed(scalar speed);

            //- Set element speed based on rotation
            void setSpeed(vector point, vector axis, scalar omega);

            //- Scale the velocity
            void scaleVelocity(scalar scale);

            //- Set dynamic stall active
            void setDynamicStallActive(bool active);

            //- Set omega for flow curvature correction
            void setOmega(scalar omega);

            //- Set end effect correction factor
            void setEndEffectFactor(scalar factor);

            //- Set velocity sample radius
            void setVelocitySampleRadius(scalar radius);

            //- Set number of velocity samples
            void setNVelocitySamples(label nSamples);


        // Evaluation

            //- Calculate forces
            void calculateForce
            (
                const volVectorField& Uin
            );

            //- Read coefficient data
            void read();

            //- Calculate and return moment about specified point
            vector moment(vector point);


        // Source term addition

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                volVectorField& force
            );

            //- Add source term to turbulence quantity
            virtual void addTurbulence(fvMatrix<scalar>& eqn, word fieldName);

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                volVectorField& force
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
