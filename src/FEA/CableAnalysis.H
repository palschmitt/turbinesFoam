/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::CableAnalysis

Description
    Finite-element analysis of cable/rope networks using a tension-only
    formulation.  Each node carries 3 translational DOF (no rotations).
    The element stiffness combines:
      - Linear elastic axial stiffness  k_e = EA/L
      - Geometric (stress) stiffness    k_g = T/L * (I - e*e^T)
    where T is the current element tension and e is the unit direction vector.

    Because cable structures are geometrically nonlinear, the solver uses a
    Newton-Raphson iteration.  Tension-only behaviour is enforced by zeroing
    the stiffness contribution of any element whose tension becomes negative
    (slack cable).

    Input conventions (matching FrameAnalysis where possible):
      nodes       – (nnodes x 3)  nodal coordinates [x y z]
      elems       – (nelems x 2)  node connectivity  [n1 n2]
      restraints  – (nnodes x 3)  1 = fixed DOF, 0 = free
      mats        – (nelems x 2)  [E  A]      (Young's modulus, cross-section area)
      loads       – (nnodes x 3)  nodal forces [Fx Fy Fz]
      prescribed  – (nnodes x 3)  prescribed displacement (applied at restrained DOF)
      pretension  – (nelems x 1)  initial element pretension T0 (may be zero)

    Outputs:
      nodedisp()   – (nnodes x 3) nodal displacements
      tensions()   – (nelems x 1) element tensions in deformed state

SourceFiles
    CableAnalysis.C

\*---------------------------------------------------------------------------*/

#ifndef CableAnalysis_H
#define CableAnalysis_H

#include <armadillo>

#include "List.H"
#include "dictionary.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of stream operators
class Istream;
class Ostream;
class CableAnalysis;
Istream& operator>>(Istream&, CableAnalysis&);
Ostream& operator<<(Ostream&, const CableAnalysis&);


/*---------------------------------------------------------------------------*\
                       Class CableAnalysis Declaration
\*---------------------------------------------------------------------------*/

class CableAnalysis
{
    // -----------------------------------------------------------------------
    // Private Data
    // -----------------------------------------------------------------------

    // Problem dimensions
    int nelems_;          //- Number of elements
    int nnodes_;          //- Number of nodes
    int nnode_;           //- Nodes per element (always 2)
    int ndof_;            //- DOF per node     (always 3 – translations only)
    int totdof_;          //- Total system DOF = nnodes_ * ndof_
    int nfree_;           //- Number of free (unconstrained) DOF

    // Mesh
    arma::Mat<double> nodes_;        //- Reference nodal coordinates (nnodes x 3)
    arma::Mat<int>    elems_;        //- Element connectivity       (nelems  x 2)
    arma::Mat<int>    restraints_;   //- Restraint flags            (nnodes  x 3)

    // Material / section
    arma::Mat<double> mats_;         //- [E  A] per element         (nelems  x 2)
    arma::Mat<double> pretension_;   //- Initial pretension T0      (nelems  x 1)

    // Loading
    arma::Mat<double> loads_;        //- Nodal forces               (nnodes  x 3)
    arma::Mat<double> prescribed_;   //- Prescribed displacements   (nnodes  x 3)

    // Current element geometry (recomputed each iteration)
    double E_;            //- Young's modulus
    double A_;            //- Cross-sectional area
    double L0_;           //- Unstressed (reference) element length
    double L_;            //- Deformed element length
    double T_;            //- Current element tension
    double Cx_, Cy_, Cz_; //- Direction cosines of deformed element

    int ielem_;           //- Current element index (loop variable)
    arma::Mat<int> elnodes_; //- Node indices of current element (2 x 1)

    // System matrices (assembled)
    arma::Mat<double> K_;            //- Tangent stiffness (full system, totdof x totdof)
    arma::Mat<double> Ke_;           //- Element stiffness (6 x 6 in global frame)
    arma::Mat<double> Kff_, Kfr_;
    arma::Mat<double> Krf_, Krr_;

    // Vectors
    arma::Mat<double> Fext_;         //- External force vector   (totdof x 1)
    arma::Mat<double> Fint_;         //- Internal force vector   (totdof x 1)
    arma::Mat<double> Ff_, Fr_;
    arma::Mat<double> defr_;         //- Prescribed displacement values (restrained DOF)
    arma::Mat<double> deff_;         //- Free-DOF displacement solution

    // DOF reordering
    arma::Mat<int> order1_;          //- Maps reordered index -> original DOF index
    arma::Mat<int> order2_;          //- Maps original DOF index -> reordered index

    // Results
    arma::Mat<double> nodedisp_;     //- Nodal displacements (nnodes x 3)
    arma::Mat<double> tensions_;     //- Element tensions    (nelems  x 1)

    // Solver parameters
    int    maxIter_;      //- Maximum Newton-Raphson iterations
    double tol_;          //- Convergence tolerance (relative residual norm)
    double minL0_;        //- Minimum reference element length (for step limiting)

    // Warm-start initial displacement (totdof x 1), default zero.
    // When supplied this is used as the starting point for Newton-Raphson
    // instead of the zero vector, so the solver begins near the previous
    // converged solution and avoids large transient overshoots.
    arma::Mat<double> u0_;


    // -----------------------------------------------------------------------
    // Private Member Functions
    // -----------------------------------------------------------------------

        //- Convert nested OpenFOAM List to Armadillo double matrix
        arma::Mat<double> List2Mat(const List<List<scalar>>& lst);

        //- Convert nested OpenFOAM int List to Armadillo int matrix
        arma::Mat<int>    List2intMat(const List<List<int>>& lst);

        //- Convert Armadillo double matrix to nested OpenFOAM List
        Foam::List<Foam::List<Foam::scalar>>
            Mat2List(const arma::Mat<double>& mat);

        //- Update current-element geometry from deformed nodal positions
        //  (nodes_ + current displacement u)
        void updateElemGeom(const arma::Mat<double>& u);

        //- Build 6x6 element tangent stiffness (elastic + geometric) in global
        //  frame for the current element.  Slack elements (T_ < 0) contribute
        //  only a small stabilisation stiffness to avoid singularity.
        void elemStiffness();

        //- Assemble element stiffness Ke_ into global K_
        void assemble();

        //- Build DOF reordering arrays (order1_, order2_) from restraints_
        void neworder();

        //- Partition reordered K_ into Kff_, Kfr_, Krf_, Krr_
        void reorder();

        //- Assemble external load vector Ff_ (free DOF) and Fr_ (restrained DOF)
        void assembleLoads();

        //- Set prescribed displacement vector defr_
        void applyPrescribed();

        //- Compute internal force vector Fint_ from current tensions
        void internalForces(const arma::Mat<double>& u);

        //- Main Newton-Raphson solve loop; populates nodedisp_ and tensions_
        void solve();


public:

    // Constructors

        //- Construct null
        CableAnalysis();

        //- Construct from components
        CableAnalysis
        (
            const List<List<scalar>>& nodes,
            const List<List<int>>&    elems,
            const List<List<int>>&    restraints,
            const List<List<scalar>>& mats,
            const List<List<scalar>>& loads,
            const List<List<scalar>>& prescribed,
            const List<List<scalar>>& pretension,
            int    maxIter = 200,
            double tol     = 1.0e-8,
            const List<List<scalar>>& u0 = List<List<scalar>>()
        );

        //- Construct from Istream
        CableAnalysis(Istream&);

        //- Copy constructor
        CableAnalysis(const CableAnalysis&);


    // Selectors

        //- Select null-constructed
        static autoPtr<CableAnalysis> New();


    //- Destructor
    ~CableAnalysis();


    // Access (results)

        //- Nodal displacement matrix (nnodes x 3)
        const arma::Mat<double>& nodedisp() const;

        //- Nodal displacement as nested OpenFOAM List
        Foam::List<Foam::List<Foam::scalar>> nodedispList();

        //- Element tensions (nelems x 1)
        const arma::Mat<double>& tensions() const;

        //- Element tensions as flat OpenFOAM List<scalar>
        Foam::List<Foam::scalar> tensionsList();


    // IOstream Operators

        friend Istream& operator>>(Istream&, CableAnalysis&);
        friend Ostream& operator<<(Ostream&, const CableAnalysis&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
