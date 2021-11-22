
set datafile separator ','

set y2tics
set ylabel 'AoA [deg]'
set xlabel 'Time [s]'
set y2label '$C_L$'
set yrange [8.2:10.1]
#set y2range [0.75:1.25]
#Beware to set element index to middle or half the max elements
plot  '../fixedAoA/postProcessing/actuatorBernoulliLineElements/0/leftblade.element30.csv' us 1:9 w l tit '$AoA_{f}$',\
'../fixedAoA/postProcessing/actuatorBernoulliLineElements/0/leftblade.element30.csv' us 1:10 axis x1y2  w l tit '$C_{L,f}$',\
'postProcessing/actuatorBernoulliLineElements/0/leftblade.element30.csv' us 1:9  w l tit '$AoA$',\
'postProcessing/actuatorBernoulliLineElements/0/leftblade.element30.csv' us 1:10  w l axis x1y2 tit '$C_L$'



