
set datafile separator ','

set y2tics
set ylabel 'AoA [deg]'
set xlabel 'Time [s]'
set y2label 'C_L'
set yrange [8.2:10.1]
set xrange [0:2]
set ytics nomirror
#set y2range [0.75:1.25]
#Beware to set element index to middle or half the max elements
set key left  at 1.4, 9.8
plot  '../fixedAoA/postProcessing/actuatorBernoulliLineElements/0/leftblade.element30.csv' us 1:9 linetype 1 w l tit 'AoA_{Case 2}',\
'../fixedAoA/postProcessing/actuatorBernoulliLineElements/0/leftblade.element30.csv' us 1:10 axis x1y2 linetype 1 dashtype 2 w l tit 'C_{L,Case 2}',\
'postProcessing/actuatorBernoulliLineElements/0/leftblade.element30.csv' us 1:9 w l linetype 2  tit 'AoA_{Case 1}',\
'postProcessing/actuatorBernoulliLineElements/0/leftblade.element30.csv' us 1:10   w l axis x1y2 linetype 2 dashtype 2 tit 'C_{L,Case 1}'
set term pdf
set out 'AoAandCl.pdf'
replot

reset
set term pdf
set out 'Wake.pdf'
set datafile separator ','
set xlabel 'Velocity [m/s]'
set ylabel 'Depth [m]'
plot 'Wakedatastiff.csv' us 1:5 w l linetype 1 tit 'Case 2',\
'Wakedatasoft.csv' us 1:5 w l  linetype 2 tit 'Case 1'

replot
