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
#include "error.H"
#include <limits>
#include <cmath>

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

    // Guard: zero reference length means coincident nodes in the input mesh –
    // this is a setup error.  Set L0_ to a small positive value so we do not
    // divide by zero, then flag excessive strain below.
    if (L0_ < 1.0e-14)
    {
        L0_ = 1.0e-14;
        WarningIn("CableAnalysis::updateElemGeom()")
            << "Element " << ielem_ << " has zero reference length (L0=0). "
            << "Check that no two consecutive nodes in the cable mesh are "
               "coincident.  Axial strain will be unreliable for this element."
            << endl;
    }

    // Guard: deformed length collapsing to zero (e.g. Newton overshoot).
    // Clamp L_ to 1 % of L0_ so direction cosines remain finite.
    // The resulting large compressive strain will make T_ strongly negative,
    // the element will be treated as slack, and the geometric stiffness will
    // resist further collapse on the next iteration.
    if (L_ < 1.0e-2 * L0_)
        L_ = 1.0e-2 * L0_;

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

    // Guard: L_ should never be zero here (updateElemGeom clamps it), but
    // protect the divisions anyway.
    double Lsafe = std::max(L_, 1.0e-14);

    // --- Elastic stiffness contribution (axial only)
    //  k_e = (EA/L) * [e*e^T , -e*e^T ; -e*e^T , e*e^T]
    double kAxial = E_ * A_ / Lsafe;

    // --- Geometric (stress) stiffness contribution
    //  k_g = (T/L) * [ I - e*e^T , -(I - e*e^T) ; -(I-e*e^T) , (I-e*e^T) ]
    // For a tension-only cable we zero k_g if cable is slack (T_ < 0).
    // A small stabilisation stiffness kSlack is retained to avoid singularity.
    double kGeo = 0.0;
    if (T_ > 0.0)
        kGeo = T_ / Lsafe;
    else
    {
        // Slack cable: remove elastic contribution too (no compression)
        kAxial = 0.0;
        // Tiny stabiliser keeps the system non-singular
        kGeo = 1.0e-10 * E_ * A_ / std::max(L0_, 1.0e-14);
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

    // Current displacement (global, all DOF).
    // Initialise from warm-start vector u0_ (supplied by caller from the
    // previous converged step).  This prevents large first-step overshoots
    // when gravity / buoyancy creates an O(T0) residual from u=0.
    arma::Mat<double> u = u0_;

    // Apply prescribed displacements at restrained DOF, overriding whatever
    // u0_ has at those positions.
    for (int id = 0; id < totdof_; id++)
    {
        int pos = order2_(id);
        if (pos >= nfree_)
            u(id) = defr_(pos - nfree_);
    }

    // Reference load norm for relative convergence check.
    // We need a force scale representative of the problem.  On incremental
    // calls from actuatorCableLineSource the external load vector contains
    // only the *change* in fluid force since the last step, which is zero
    // when the flow is steady.  In that case norm(Fext_)=0 and the old
    // guard (fNorm=1) gave an absolute tolerance of tol_=1e-8 N -- far
    // tighter than needed and impossible to reach when pretension creates
    // O(T0) internal forces at the free nodes.
    //
    // Correct approach: use the larger of
    //   (a) norm(Fext_)           -- external load scale
    //   (b) norm(Fint at u=0)     -- internal prestress scale
    // so the relative tolerance is always meaningful.
    //
    // Compute Fint at u=0 (reference internal forces from pretension).
    {
        arma::Mat<double> u0;
        u0.zeros(totdof_, 1);
        internalForces(u0);
    }
    double fNormExt  = arma::norm(Fext_);
    double fNormInt  = arma::norm(Fint_);
    double fNorm     = std::max(fNormExt, fNormInt);
    if (fNorm < 1.0e-30) fNorm = 1.0;

    bool converged = false;

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

        // -- 5. Check convergence (every iteration, including iter == 0)
        double resNorm = arma::norm(Rf);
        if (resNorm / fNorm < tol_)
        {
            converged = true;
            break;
        }

        // -- 6. Solve for displacement increment.
        //    Primary attempt: direct solve.
        //    If Kff_ is singular (e.g. all-slack cable at t=0) armadillo
        //    returns false rather than throwing when we use the bool-return
        //    overload.  On failure, apply Tikhonov regularisation scaled to
        //    the diagonal magnitude and retry.  A second failure is a genuine
        //    problem and reported via OpenFOAM FatalError so that all MPI
        //    ranks receive a clean abort message rather than an uncaught C++
        //    exception that kills ranks independently.
        arma::Mat<double> du(nfree_, 1);
        bool solveOk = arma::solve(du, Kff_, Rf,
                                   arma::solve_opts::no_approx +
                                   arma::solve_opts::equilibrate);
        if (!solveOk)
        {
            double alpha = 1.0e-6 * arma::abs(Kff_.diag()).max();
            if (alpha < 1.0e-30) alpha = 1.0;
            arma::Mat<double> Kreg = Kff_;
            Kreg.diag() += alpha;
            if (!arma::solve(du, Kreg, Rf))
            {
                FatalErrorIn("CableAnalysis::solve()")
                    << "Tangent stiffness matrix is singular and Tikhonov "
                       "regularisation failed at Newton-Raphson iteration "
                    << iter << ".\n"
                    << "Check that at least one node per free-DOF direction "
                       "is restrained and that EA > 0 for all elements.\n"
                    << "nfree=" << nfree_
                    << "  nelems=" << nelems_
                    << "  nnodes=" << nnodes_
                    << abort(FatalError);
            }
        }

        // -- 7. Update free DOF displacements.
        //    Limit the step so no free DOF moves more than half the shortest
        //    reference element length in a single iteration.  This prevents
        //    Newton overshoot from collapsing elements to zero length on the
        //    first iteration when pretension dominates the residual.
        double duMax = arma::abs(du).max();
        if (duMax > 0.5 * minL0_)
            du *= (0.5 * minL0_) / duMax;

        for (int id = 0; id < totdof_; id++)
        {
            int pos = order2_(id);
            if (pos < nfree_)
                u(id) += du(pos);
        }
    }

    if (!converged)
    {
        WarningIn("CableAnalysis::solve()")
            << "Newton-Raphson did not converge in " << maxIter_
            << " iterations (tol=" << tol_ << ").\n"
            << "Continuing with last iterate.  Consider increasing "
               "CableMaxIter or relaxing CableTolerance, or verify that "
               "boundary conditions fully constrain rigid-body motion."
            << endl;
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
    totdof_(0), nfree_(0), maxIter_(200), tol_(1.0e-8), minL0_(1.0)
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
    double tol,
    const List<List<scalar>>& u0
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

    // Warm-start displacement.  If u0 is supplied and has the right size,
    // store it; otherwise initialise to zero.
    u0_.zeros(totdof_, 1);
    if (u0.size() == nnodes_)
    {
        for (int in = 0; in < nnodes_; in++)
            if (u0[in].size() == ndof_)
                for (int id = 0; id < ndof_; id++)
                    u0_(in*ndof_ + id) = u0[in][id];
    }

    // Build DOF ordering from restraints
    neworder();

    // Compute minimum reference element length (used for step limiting in solve)
    minL0_ = std::numeric_limits<double>::max();
    for (int ie = 0; ie < nelems_; ie++)
    {
        int n1 = elems_(ie, 0);
        int n2 = elems_(ie, 1);
        double dx = nodes_(n2,0) - nodes_(n1,0);
        double dy = nodes_(n2,1) - nodes_(n1,1);
        double dz = nodes_(n2,2) - nodes_(n1,2);
        double L0 = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (L0 > 1.0e-14 && L0 < minL0_) minL0_ = L0;
    }
    if (minL0_ > 1.0e29) minL0_ = 1.0;  // fallback if all lengths are zero

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
