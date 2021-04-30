%Assuems elliptic shape to provide moments..
chordlength=0.07
thicknessratio=0.12
Ix=pi*0.5*chordlength*(0.5*chordlength*thicknessratio)^3/4
Iy=pi*(0.5*chordlength)^3*0.5*chordlength*thicknessratio/4
It=pi/16*chordlength^3*(thicknessratio*chordlength)^3/(chordlength^2+(thicknessratio*chordlength)^2)
A=pi/4*chordlength^2*thicknessratio

disp([A Ix Iy It ])