set term gif animate
set output 'animacioncil.gift'
#set xra[-0.5:0.5]
#set yra[-0.35:0.35]


set table 'absx.tmp'
    plot 'posicionesneg' using (abs($1)) with table
unset table

stats 'absx.tmp' using 1 name "p" nooutput

print "Máximo abs(x):", p_max
print "Mínimo abs(x):", p_min

x_max=p_max
x_min=p_min*1000

set xra[-x_max-0.2:x_max+0.2]
set yra[-x_max-0.05:x_max+0.05]


 set parametric
do for[a=0:1200:5]{plot 'posicionesneg' i a using 2:1 t 'Electrones' lc rgb 'blue','posicionespos' i a using 2:1 t 'Cationes' lc 'red',  x_max*cos(t),x_max*sin(t) lc rgb 'blue' lw 3 t 'Cátodo(-)',  x_min*cos(t),x_min*sin(t) lc rgb 'red' lw 3 t 'Ánodo(+)'}
