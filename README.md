ALFEA
Actuator Line modelling  extended with Finite Element Analysis
adapted from turbinesFoam
============

![OpenFOAM 8](https://img.shields.io/badge/OpenFOAM-8-brightgreen.svg)
![OpenFOAM 7](https://img.shields.io/badge/OpenFOAM-7-brightgreen.svg)

ALFEA adds Finite Element Analysis capabilities to turbinesFoam, a library for simulating wind and marine hydrokinetic turbines
in OpenFOAM using the actuator line method.

[![](https://cloud.githubusercontent.com/assets/4604869/10141523/f2e3ad9a-65da-11e5-971c-b736abd30c3b.png)](https://www.youtube.com/watch?v=THZvV4R1vow)

Be sure to check out the
[development snapshot videos on YouTube](https://www.youtube.com/playlist?list=PLOlLyh5gytG8n8D3V1lDeZ3e9fJf9ux-e).


Contributing
------------

Pull requests are very welcome!
See the [issue tracker](https://github.com/petebachant/turbinesFoam/issues)
for more details.


Features
--------

`fvOptions` classes for adding actuator lines and turbines constructed from
actuator lines to any compatible solver or turbulence model, e.g.,
`simpleFoam`, `pimpleFoam`, `interFoam`, etc.


Installation
------------
Tested in OF-8 and OF-7.

Requires the Armadillo library:
http://arma.sourceforge.net/

Before installing Armadillo, first install OpenBLAS and LAPACK, along with the corresponding development/header files
     
Recommended packages to install before installing Armadillo:
        Fedora & Red Hat: cmake, openblas-devel, lapack-devel, arpack-devel, SuperLU-devel
        Ubuntu & Debian: cmake, libopenblas-dev, liblapack-dev, libarpack2-dev, libsuperlu-dev

Pre-built Armadillo packages are provided by many Linux-based operating systems: Fedora, Debian, Ubuntu, openSUSE, Arch
the pre-built packages may not be the latest version; if you're encountering problems, use the official stable version provided here


```bash
cd $WM_PROJECT_USER_DIR
git clone https://github.com/turbinesFoam/turbinesFoam.git
cd turbinesFoam
./Allwmake
```


Usage
-----

There are tutorials located in `turbinesFoam/tutorials`.
Each case contains an Allrun and Allclean script. 
Results can be plotted against analytical solutions using the GNUoctave toolbox and the ComparAnalyticalsolution.m script files.

Test cases for ALFEA:

tutorials/actuatorLine/bendingtestFEA
 demonstrates bending deformation of a flexible wing 

tutorials/actuatorLine/torsiontestFEA
 demonstrates torsional deformation of a beam supporting a wing

tutorials/actuatorLine/changeAoA
 demonstrates change in angle of attack by the bending deformation of the support structure of a wing
 fixedAoA runs the same case but for a "stiff" support structure for comprison



Publications
------------

Bachant, P., Goude, A., and Wosnik, M. (2016) [_Actuator line modeling of vertical-axis turbines_](https://arxiv.org/abs/1605.01449). arXiv preprint 1605.01449.


How to cite
-----------

The latest release of turbinesFoam can be cited via DOI thanks to Zenodo: [![DOI](https://zenodo.org/badge/4234/turbinesFoam/turbinesFoam.svg)](https://zenodo.org/badge/latestdoi/4234/turbinesFoam/turbinesFoam)


Acknowledgements
----------------
The Bryden Centre project funded Pa\'l Schmitt and is supported by the European 	Union’s INTERREG VA Programme, managed by the Special EU Programmes Body (SEUPB).
Computations were performed on the Kelvin-2 HPC system based at QUB, funded by EPSRC grant EP/T022175/1.
%  
\section*{DISCLAIMER}
The views and opinions expressed in this paper do not necessarily reflect those of the European Commission or the Special EU Programmes Body (SEUPB).



This work was funded through a National Science Foundation CAREER award,
principal investigator Martin Wosnik ([NSF CBET
1150797](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1150797), Energy for
Sustainability, original program manager Geoffrey A. Prentice, current program
manager Gregory L. Rorrer).

OpenFOAM is free, open source software for computational fluid dynamics (CFD),
developed primarily by [CFD Direct](http://cfd.direct), on behalf of the
[OpenFOAM](http://openfoam.org) Foundation.

Interpolation, Gaussian projection, and vector rotation functions adapted from
NREL's [SOWFA](https://github.com/NREL/SOWFA).
