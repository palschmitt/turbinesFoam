/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.  (GPL v3 – see header in .H file)

\*---------------------------------------------------------------------------*/

#include "CableAnalysis.H"

// * * * * * * * * * * * * * Private Helper Functions * * * * * * * * * * * //

arma::Mat<double>
Foam::CableAnalysis::List2Mat(const List<List<scalar>>& lst)
{
    arma::Mat<double> mat(lst.size(), lst[0].size());
    forAll(lst, i)
        forAll(lst[i], j)
            mat(i, j) = lst[i][j];
    return mat;
}

arma::Mat<int>
Foam::CableAnalysis::List2intMat(const List<List<int>>& lst)
{
    arma::Mat<int> mat(lst.size(), lst[0].size());
    forAll(lst, i)
        forAll(lst[i], j)
            mat(i, j) = lst[i][j];
    return mat;
}

Foam::List<Foam::List<Foam::scalar>>
Foam::CableAnalysis::Mat2List(const arma::Mat<double>& mat)
{
    List<List<scalar>> out(mat.n_rows);
    for (uint i = 0; i < mat.n_rows; i++)
    {
        List<scalar> row(mat.n_cols);
        for (uint j = 0; j < mat.n_cols; j++)
            row[j] = mat(i, j);
        out[i] = row;
    }
    return out;
}


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::CableAnalysis::updateElemGeom(const arma::Mat<double>& u)
{
    // Node indices for current element
    int n1 = elnodes_(0);
    int n2 = elnodes_(1);

    // Deformed nodal positions
    arma::vec x1 = {nodes_(n1,0) + u(n1*ndof_+0),
                    nodes_(n1,1) + u(n1*ndof_+1),
                    nodes_(n1,2) + u(n1*ndof_+2)};
    arma::vec x2 = {nodes_(n2,0) + u(n2*ndof_+0),
                    nodes_(n2,1) + u(n2*ndof_+1),
                    nodes_(n2,2) + u(n2*ndof_+2)};

    arma::vec dv = x2 - x1;
    L_ = arma::norm(dv);

    // Unstressed length from reference configuration
    arma::vec dv0 = {nodes_(n2,0) - nodes_(n1,0),
                     nodes_(n2,1) - nodes_(n1,1),
                     nodes_(n2,2) - nodes_(n1,2)};
    L0_ = arma::norm(dv0);

    // Unit direction cosines (deformed)
    Cx_ = dv(0) / L_;
    Cy_ = dv(1) / L_;
    Cz_ = dv(2) / L_;

    // Current tension: elastic stretch + pretension
    double T0 = pretension_(ielem_, 0);
    double strain = (L_ - L0_) / L0_;
    T_ = E_ * A_ * strain + T0;
}


void Foam::CableAnalysis::elemStiffness()
{
    // --- Direction cosine vector e = [Cx Cy Cz]
    arma::vec e = {Cx_, Cy_, Cz_};

    // --- Elastic stiffness contribution (axial only)
    //  k_e = (EA/L) * [e*e^T , -e*e^T ; -e*e^T , e*e^T]
    double kAxial = E_ * A_ / L_;

    // --- Geometric (stress) stiffness contribution
    //  k_g = (T/L) * [ I - e*e^T , -(I - e*e^T) ; -(I-e*e^T) , (I-e*e^T) ]
    // For a tension-only cable we zero k_g if cable is slack (T_ < 0).
    // A small stabilisation stiffness kSlack is retained to avoid singularity.
    double kGeo = 0.0;
    if (T_ > 0.0)
        kGeo = T_ / L_;
    else
    {
        // Slack cable: remove elastic contribution too (no compression)
        kAxial = 0.0;
        // Tiny stabiliser keeps the system non-singular
        kGeo = 1.0e-10 * E_ * A_ / L0_;
    }

    arma::mat I3 = arma::eye(3, 3);
    arma::mat eeT = e * e.t();

    // 3x3 sub-block contributions
    arma::mat Kaa =  kAxial * eeT + kGeo * (I3 - eeT);
    arma::mat Kab = -kAxial * eeT - kGeo * (I3 - eeT);
    arma::mat Kba =  Kab;
    arma::mat Kbb =  Kaa;

    // Assemble 6x6 element stiffness in global frame
    Ke_.zeros(6, 6);
    Ke_.submat(0,0,2,2) = Kaa;
    Ke_.submat(0,3,2,5) = Kab;
    Ke_.submat(3,0,5,2) = Kba;
    Ke_.submat(3,3,5,5) = Kbb;
}


void Foam::CableAnalysis::assemble()
{
    int n1 = elnodes_(0);
    int n2 = elnodes_(1);
    int ids[2] = {n1, n2};

    for (int i = 0; i < nnode_; i++)
    {
        int gi = ids[i] * ndof_;
        int li = i      * ndof_;
        for (int j = 0; j < nnode_; j++)
        {
            int gj = ids[j] * ndof_;
            int lj = j      * ndof_;
            for (int p = 0; p < ndof_; p++)
                for (int q = 0; q < ndof_; q++)
                    K_(gi+p, gj+q) += Ke_(li+p, lj+q);
        }
    }
}


void Foam::CableAnalysis::neworder()
{
    totdof_ = nnodes_ * ndof_;
    order1_.zeros(1, totdof_);
    order2_.zeros(1, totdof_);

    arma::Mat<int> restraintList;
    restraintList.zeros(1, totdof_);
    int cnt = 0;
    int nRestrained = 0;

    for (int in = 0; in < nnodes_; in++)
        for (int id = 0; id < ndof_; id++)
        {
            restraintList(cnt) = restraints_(in, id);
            if (restraints_(in, id) == 1) nRestrained++;
            cnt++;
        }

    nfree_ = totdof_ - nRestrained;

    int iFree = 0, iRes = 0;
    for (int id = 0; id < totdof_; id++)
    {
        if (restraintList(id) == 1)
        {
            order1_(nfree_ + iRes) = id;
            order2_(id)            = nfree_ + iRes;
            iRes++;
        }
        else
        {
            order1_(iFree) = id;
            order2_(id)    = iFree;
            iFree++;
        }
    }
}


void Foam::CableAnalysis::reorder()
{
    Kff_.zeros(nfree_, nfree_);
    Kfr_.zeros(nfree_, totdof_ - nfree_);
    Krf_.zeros(totdof_ - nfree_, nfree_);
    Krr_.zeros(totdof_ - nfree_, totdof_ - nfree_);

    for (int i = 0; i < totdof_; i++)
    {
        int pi = order2_(i);
        for (int j = 0; j < totdof_; j++)
        {
            int pj = order2_(j);
            if      (pi   < nfree_ && pj   < nfree_) Kff_(pi,        pj)        += K_(i,j);
            else if (pi   < nfree_ && pj  >= nfree_) Kfr_(pi,        pj-nfree_) += K_(i,j);
            else if (pi  >= nfree_ && pj   < nfree_) Krf_(pi-nfree_, pj)        += K_(i,j);
            else                                      Krr_(pi-nfree_, pj-nfree_) += K_(i,j);
        }
    }
}


void Foam::CableAnalysis::assembleLoads()
{
    Fext_.zeros(totdof_, 1);
    int cnt = 0;
    for (int in = 0; in < nnodes_; in++)
        for (int id = 0; id < ndof_; id++)
            Fext_(cnt++) = loads_(in, id);

    Ff_.zeros(nfree_, 1);
    Fr_.zeros(totdof_ - nfree_, 1);
    for (int id = 0; id < totdof_; id++)
    {
        int pos = order2_(id);
        if (pos < nfree_)
            Ff_(pos) = Fext_(id);
        else
            Fr_(pos - nfree_) = Fext_(id);
    }
}


void Foam::CableAnalysis::applyPrescribed()
{
    defr_.zeros(totdof_ - nfree_, 1);
    arma::Mat<double> dv;
    dv.zeros(totdof_, 1);
    int cnt = 0;
    for (int in = 0; in < nnodes_; in++)
        for (int id = 0; id < ndof_; id++)
            dv(cnt++) = prescribed_(in, id);

    for (int id = 0; id < totdof_; id++)
    {
        int pos = order2_(id);
        if (pos >= nfree_)
            defr_(pos - nfree_) = dv(id);
    }
}


void Foam::CableAnalysis::internalForces(const arma::Mat<double>& u)
{
    // Assemble internal force vector from current element tensions
    Fint_.zeros(totdof_, 1);

    for (ielem_ = 0; ielem_ < nelems_; ielem_++)
    {
        E_  = mats_(ielem_, 0);
        A_  = mats_(ielem_, 1);
        elnodes_.zeros(2, 1);
        elnodes_(0) = elems_(ielem_, 0);
        elnodes_(1) = elems_(ielem_, 1);

        updateElemGeom(u);

        arma::vec e = {Cx_, Cy_, Cz_};

        // Internal force contribution: T * e at node 2, -T * e at node 1
        // (only non-slack elements contribute tension)
        double Teff = (T_ > 0.0) ? T_ : 0.0;

        int n1 = elnodes_(0);
        int n2 = elnodes_(1);
        for (int d = 0; d < ndof_; d++)
        {
            Fint_(n1*ndof_ + d) -= Teff * e(d);
            Fint_(n2*ndof_ + d) += Teff * e(d);
        }
    }
}


void Foam::CableAnalysis::solve()
{
    // -----------------------------------------------------------------------
    // Newton-Raphson iteration for geometrically nonlinear cable analysis
    // -----------------------------------------------------------------------

    // Current displacement (global, all DOF)
    arma::Mat<double> u;
    u.zeros(totdof_, 1);

    // Apply prescribed displacements immediately
    for (int id = 0; id < totdof_; id++)
    {
        int pos = order2_(id);
        if (pos >= nfree_)
            u(id) = defr_(pos - nfree_);
    }

    for (int iter = 0; iter < maxIter_; iter++)
    {
        // -- 1. Assemble tangent stiffness from deformed geometry
        K_.zeros(totdof_, totdof_);
        for (ielem_ = 0; ielem_ < nelems_; ielem_++)
        {
            E_  = mats_(ielem_, 0);
            A_  = mats_(ielem_, 1);
            elnodes_.zeros(2, 1);
            elnodes_(0) = elems_(ielem_, 0);
            elnodes_(1) = elems_(ielem_, 1);
            updateElemGeom(u);
            elemStiffness();
            assemble();
        }

        // -- 2. Partition stiffness
        reorder();

        // -- 3. Internal force vector
        internalForces(u);

        // -- 4. Residual  R = Fext - Fint  (free DOF only)
        arma::Mat<double> Rfull = Fext_ - Fint_;

        arma::Mat<double> Rf(nfree_, 1);
        for (int id = 0; id < totdof_; id++)
        {
            int pos = order2_(id);
            if (pos < nfree_)
                Rf(pos) = Rfull(id);
        }
        // Subtract coupling from prescribed DOF
        Rf -= Kfr_ * defr_;

        // -- 5. Check convergence
        double resNorm = arma::norm(Rf);
        double fNorm   = arma::norm(Ff_);
        if (fNorm < 1.0e-30) fNorm = 1.0;   // guard against zero-load case
        if (resNorm / fNorm < tol_ && iter > 0)
            break;

        // -- 6. Solve for displacement increment
        arma::Mat<double> du = arma::solve(
            Kff_, Rf, arma::solve_opts::no_approx);

        // -- 7. Update free DOF displacements
        for (int id = 0; id < totdof_; id++)
        {
            int pos = order2_(id);
            if (pos < nfree_)
                u(id) += du(pos);
        }
    }

    // -----------------------------------------------------------------------
    // Extract results
    // -----------------------------------------------------------------------

    // Nodal displacements reshaped to (nnodes x 3)
    nodedisp_.set_size(nnodes_, ndof_);
    for (int in = 0; in < nnodes_; in++)
        for (int id = 0; id < ndof_; id++)
            nodedisp_(in, id) = u(in*ndof_ + id);

    // Element tensions in deformed state
    tensions_.set_size(nelems_, 1);
    for (ielem_ = 0; ielem_ < nelems_; ielem_++)
    {
        E_  = mats_(ielem_, 0);
        A_  = mats_(ielem_, 1);
        elnodes_.zeros(2, 1);
        elnodes_(0) = elems_(ielem_, 0);
        elnodes_(1) = elems_(ielem_, 1);
        updateElemGeom(u);
        tensions_(ielem_) = (T_ > 0.0) ? T_ : 0.0;  // slack = zero tension
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CableAnalysis::CableAnalysis()
:
    nelems_(0), nnodes_(0), nnode_(2), ndof_(3),
    totdof_(0), nfree_(0), maxIter_(200), tol_(1.0e-8)
{}


Foam::CableAnalysis::CableAnalysis
(
    const List<List<scalar>>& nodes,
    const List<List<int>>&    elems,
    const List<List<int>>&    restraints,
    const List<List<scalar>>& mats,
    const List<List<scalar>>& loads,
    const List<List<scalar>>& prescribed,
    const List<List<scalar>>& pretension,
    int    maxIter,
    double tol
)
:
    nnode_(2),
    ndof_(3),
    maxIter_(maxIter),
    tol_(tol)
{
    nelems_      = elems.size();
    nnodes_      = nodes.size();
    totdof_      = nnodes_ * ndof_;

    nodes_       = List2Mat(nodes);
    elems_       = List2intMat(elems);
    restraints_  = List2intMat(restraints);
    mats_        = List2Mat(mats);
    loads_       = List2Mat(loads);
    prescribed_  = List2Mat(prescribed);
    pretension_  = List2Mat(pretension);

    // Build DOF ordering from restraints
    neworder();

    // Assemble load vectors
    assembleLoads();
    applyPrescribed();

    // Run nonlinear solver
    solve();
}


Foam::CableAnalysis::CableAnalysis(const CableAnalysis&)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::CableAnalysis>
Foam::CableAnalysis::New()
{
    return autoPtr<CableAnalysis>(new CableAnalysis);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::CableAnalysis::~CableAnalysis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const arma::Mat<double>& Foam::CableAnalysis::nodedisp() const
{
    return nodedisp_;
}

Foam::List<Foam::List<Foam::scalar>>
Foam::CableAnalysis::nodedispList()
{
    return Mat2List(nodedisp_);
}

const arma::Mat<double>& Foam::CableAnalysis::tensions() const
{
    return tensions_;
}

Foam::List<Foam::scalar>
Foam::CableAnalysis::tensionsList()
{
    List<scalar> out(nelems_);
    for (int i = 0; i < nelems_; i++)
        out[i] = tensions_(i);
    return out;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// (Implement operator>> / operator<< as needed for dictionary I/O)

// ************************************************************************* //
