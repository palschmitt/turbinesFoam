OpenFOAM/turbinesFoam actuator line bending simulation
======================================================

This case  simulates a single vertical actuator line experiencing  bending deformation from lift forces as shown in Figure 3 in the accompanying paper.
The case is setup and run by executing 

./Allrun

Figure 4, showing analytical solutions against the simulation results, can be generated using the GNUoctave toolbox and the ComparAnalyticalsolution.m script file with the command
 
octave ComparAnalyticalsolution.m


