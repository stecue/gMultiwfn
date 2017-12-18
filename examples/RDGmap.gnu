set terminal postscript landscape enhanced color 'Helvetica' 20
set encoding iso_8859_1
set output 'RDGmap.ps' 
set key 
set ylabel 'RDG (a.u)' font "Helvetica, 20" 
set xlabel 'sign({/Symbol l}_2){/Symbol r} (a.u.)' font "Helvetica, 20"
set pm3d map
set palette defined (-0.035 "blue",-0.0075 "green", 0.02 "red")
set format y "%.2f"
set format x "%.2f"
set format cb "%.3f"
set border lw 2
set xtic  -0.05,0.01,0.05 nomirror rotate font "Helvetica"
set ytic   0.0,0.2,2.0 nomirror font "Helvetica"
set cbtic  -0.035,0.005,0.02 nomirror font "Helvetica"
set xrange [-0.05:0.05]  
set yrange [0.0:2.0] 
set cbrange [-0.035:0.02]
plot 'output.txt' u 1:2:1 with points pointtype 31 pointsize 0.3 palette t ''  



