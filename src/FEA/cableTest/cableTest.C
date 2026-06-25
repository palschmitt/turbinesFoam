
#include "argList.H"
#include "List.H"
#include "vector.H"
#include "IOstreams.H"
#include "CableAnalysis.H"
#include <iomanip>
using namespace Foam;
#include <cmath>
#include <algorithm>
void printDeformedAscii
(
    const Foam::List<Foam::List<Foam::scalar>>& nodes,
    const Foam::List<Foam::List<Foam::scalar>>& disp
)
{
    using namespace Foam;

    const label n = nodes.size();

    List<scalar> z(n), x(n);

    scalar zMin = GREAT, zMax = -GREAT;
    scalar xMin = GREAT, xMax = -GREAT;

    // Build coordinates
    for (label i = 0; i < n; ++i)
    {
        z[i] = nodes[i][2] + disp[i][2];
        x[i] = nodes[i][0] + disp[i][0];

        zMin = min(zMin, z[i]);
        zMax = max(zMax, z[i]);

        xMin = min(xMin, x[i]);
        xMax = max(xMax, x[i]);
    }

    // Find max displacement in x
    scalar maxAbsX = 0;
    for (label i = 0; i < n; ++i)
    {
        maxAbsX = max(maxAbsX, mag(x[i]));
    }

    // Force symmetric axis for clarity
    xMin = -maxAbsX;
    xMax =  maxAbsX;

    const label height = 25;
    const label width  = 60;

    List<string> canvas(height, string(width, ' '));

    // Mapping
    auto zToRow = [&](scalar zz)
    {
        scalar t = (zz - zMin)/(zMax - zMin);
        return label((1.0 - t)*(height - 1));
    };

    auto xToCol = [&](scalar xx)
    {
        scalar t = (xx - xMin)/(xMax - xMin);
        return label(t*(width - 1));
    };

    // Draw interpolated curve
    for (label i = 0; i < n-1; ++i)
    {
        for (label s = 0; s <= 20; ++s)
        {
            scalar t = scalar(s)/20.0;

            scalar zi = (1.0 - t)*z[i] + t*z[i+1];
            scalar xi = (1.0 - t)*x[i] + t*x[i+1];

            label r = zToRow(zi);
            label c = xToCol(xi);

            if (r >= 0 && r < height && c >= 0 && c < width)
                canvas[r][c] = '.';
        }
    }

    // Draw nodes
    for (label i = 0; i < n; ++i)
    {
        label r = zToRow(z[i]);
        label c = xToCol(x[i]);

        if (r >= 0 && r < height && c >= 0 && c < width)
            canvas[r][c] = 'o';
    }

    // Print plot
    Info<< "\nDeformed cable (x-z projection)\n\n";

    for (label r = 0; r < height; ++r)
    {
        scalar zHere = zMax - scalar(r)/(height-1)*(zMax - zMin);

        std::ostringstream zss;
        zss << std::fixed << std::setw(8) << std::setprecision(3) << zHere;

        Info<< "z=" << zss.str() << " |" << canvas[r] << "|\n";
    }

    // === X AXIS ===

    string axis(width, '-');

    // mark x = 0
    label zeroCol = xToCol(0.0);
    if (zeroCol >=0 && zeroCol < width)
        axis[zeroCol] = '+';

    // mark max displacement
    label maxCol = xToCol(maxAbsX);
    if (maxCol >=0 && maxCol < width)
        axis[maxCol] = '|';

    Info<< "          +" << string(width, '-') << "+\n";
    Info<< "          |" << axis << "|\n";

    // labels
    std::ostringstream left, mid, right;
    left  << std::fixed << std::setprecision(3) << xMin;
    mid   << "0.000";
    right << std::fixed << std::setprecision(3) << xMax;

    Info<< "          " 
        << left.str()
        << std::string(width/2 - left.str().size(), ' ')
        << mid.str()
        << std::string(width/2 - mid.str().size(), ' ')
        << right.str()
        << nl << nl;
}


int main(int argc, char *argv[])
{
    argList args(argc, argv);

    const label nElems = 20; // <-- change this only
    const label nNodes = nElems + 1;

    List<List<scalar>> nodes(nNodes);
    for (label i = 0; i < nNodes; ++i)
    {
        nodes[i].setSize(3);
        nodes[i][0] = 0.0;
        nodes[i][1] = 0.0;
        nodes[i][2] = scalar(i) * 0.1;
    }

    List<List<int>> elems(nElems);
    for (label e = 0; e < nElems; ++e)
    {
        elems[e].setSize(2);
        elems[e][0] = e;
        elems[e][1] = e + 1;
    }

    List<List<int>> restraints(nNodes);
    for (label i = 0; i < nNodes; ++i)
    {
        restraints[i].setSize(3);
        restraints[i][0] = 0;
        restraints[i][1] = 0;
        restraints[i][2] = 0;
    }
    restraints[0][0] = 1;
    restraints[0][1] = 1;
    restraints[0][2] = 1;

    List<List<scalar>> mats(nElems);
    for (label e = 0; e < nElems; ++e)
    {
        mats[e].setSize(2);
        mats[e][0] = 5e2;
        mats[e][1] = 1;
    }

    List<List<scalar>> loads(nNodes);
    for (label i = 0; i < nNodes; ++i)
    {
        loads[i].setSize(3);
        loads[i][0] = 10.0;
        loads[i][1] = 0.0;
        loads[i][2] = 10.0;
    }

    List<List<scalar>> prescribed(nNodes);
    for (label i = 0; i < nNodes; ++i)
    {
        prescribed[i].setSize(3);
        prescribed[i][0] = 0.0;
        prescribed[i][1] = 0.0;
        prescribed[i][2] = 0.0;
    }

    List<List<scalar>> pretension(nElems);
    for (label e = 0; e < nElems; ++e)
    {
        pretension[e].setSize(1);
        pretension[e][0] = 0.1;
    }

    List<List<scalar>> u0(nNodes);
    for (label i = 0; i < nNodes; ++i)
    {
        u0[i].setSize(3);
        u0[i][0] = 0.0;
        u0[i][1] = 0.0;
        u0[i][2] = 0.0;
    }

    const label maxIter = 200;
    const scalar tol = 1.0e-8;

    CableAnalysis cable(
        nodes,
        elems,
        restraints,
        mats,
        loads,
        prescribed,
        pretension,
        maxIter,
        tol,
        u0);

    Info<< "Run complete with " << nElems << " elements." << endl;
    List<List<scalar>> disp = cable.nodedispList();
	List<scalar> tension = cable.tensionsList();
    printDeformedAscii(nodes, disp);
    
    Info<< "\nNodal displacements [m]:\n";

	for (label i = 0; i < disp.size(); ++i)
	{
    Info<< "node " << i << ": ("
        << disp[i][0] << ", "
        << disp[i][1] << ", "
        << disp[i][2] << ")\n";
	}

    
    
    return 0;
    
    
    
}
