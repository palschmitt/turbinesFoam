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
#include "actuatorLineElement.H"
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

class actuatorBernoulliLineElement: public actuatorLineElement{

protected:

    // Protected data


        //- Location of element's start point
        vector P1_;
        //- Location of element's end point
        vector P2_;

       
        //- Structural force on element -- vector that is perpendicular to
        //  spanDirection,
        vector structforceVector_;
        //- Structural moment on element -- 
        vector structmomentVector_;
        //Total applied displacement from FEA
        vector deformation_;
        //Internal (shear) force from FEA
        vector FEAforce_;
        //Internal (bending) moment from FEA
        vector FEAmoment_;
        //- Structural stiffness of element 
        List<scalar> FEAmaterial_;// E  Poisson
        List<scalar> FEAsects_;//Section data A        Iz       Iy          J        alpha
        List<int> FEArestraints_;//restrained X Y Z Rx Ry Rz 0 or 1

     


    // Protected Member Functions

       


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

     
            
            //- Return structural force vector from element
            const vector& structforce();
            //- Return structural moment vector from element
            const vector& structmoment();
            //- Return last applied deformation
            const vector& deformation();
            //- Return stiffness vector from element
            const List<scalar>& FEAsects();
            const List<scalar>& FEAmaterial();
            const List<int>& FEArestraints();

          
            
            
            const vector pitchingMoment();


			const vector&  chordDirection();
			
			const vector&  planformNormal();
			
			const vector&  spanDirection();
        // Manipulation
			void setSpanLength(scalar);
			void setSpanDirection(vector);
			void deltaChordDirection(scalar);
			
			void setPosition(vector NewPosition);

            scalar inflowRefAngle();

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
            
            void setDeformation(vector Deformation);
            
            void setFEASects(List<scalar> FEAsects);            
            void setFEAMaterial(List<scalar> FEAmaterial);            
            void setFEARestraint(List<int> FEArestraint);            

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
           

        // Evaluation

        //- Read coefficient data
        void read();
            
        void calculateForce
            (
                const volVectorField& Uin
            );

       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
