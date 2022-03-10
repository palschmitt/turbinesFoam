OpenFOAM/turbinesFoam actuator line torsion simulation
======================================================

This case  simulates a  vertical actuator line experiencing torsional  deformation from lift forces acting on a second acturator line at a 90deg angle as shown in Figure 5 in the accompanying paper.
The case is setup and run by executing 

./Allrun

Figure 6, showing analytical solutions against the simulation results, can be generated using the GNUoctave toolbox and the ComparAnalyticalsolution.m script file with the command

octave ComparAnalyticalsolution.m



