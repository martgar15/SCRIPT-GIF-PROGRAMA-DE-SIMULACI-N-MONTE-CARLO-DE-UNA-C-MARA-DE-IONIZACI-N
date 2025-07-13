set term gif animate font 'Arial,18'
set output 'animacionpp.gift'

set yra[-0.65:0.75]

set xlabel 'Eje Y'
set ylabel 'Eje Z'
#posicion leyenda(key)
set key top right

stats 'posicionesneg' using 2 name "p" nooutput
x_max=p_max

set xra[-0.15:x_max+0.15]
set xtics 0.1

set arrow from 0,-0.5 to 0,0.5 nohead lc rgb "red" lw 5 
set arrow from x_max,-0.5 to x_max,0.5 nohead lc rgb "blue" lw 5 
do for[a=0:380:1]{plot 'posicionespos' i a using 2:1 t 'Iones positivos (+)' lc rgb 'red' ps 0.7,'posicionesneg' i a using 2:1 t 'Iones negativos (-)' lc 'blue' ps 0.7}
