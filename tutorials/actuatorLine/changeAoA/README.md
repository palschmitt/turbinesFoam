OpenFOAM/turbinesFoam changing AoA case
==============================================

This case  simulates a suspended actuator line experiencing  bending deformation of its support structure and as a consequence changes in AoA of the main actuator line sectionas shown in Figure 7 in the accompanying paper.
The corresponding stiff case can be found under ../fixedAoA.
Both cases are setup and run by executing their respective 

./Allrun

scripts.

Figure 8 showing analytical solutions against the simulation results for deformation of the support structure  can be generated using the GNUoctave toolbox and the ComparAnalyticalsolution.m script file with the command

octave ComparAnalyticalsolution.m

Figure 9 showing the change in AoA and coefficnet of lift over time can be generated using the  gnuplot command file plotAoAandCL.plt, 

gnuplot plotAoAandCL.plt



