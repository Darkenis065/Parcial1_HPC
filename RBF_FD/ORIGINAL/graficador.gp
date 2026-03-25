#!/usr/bin/env gnuplot
# ==============================================================================
# SCRIPT DE GNUPLOT CORREGIDO
# ==============================================================================

set terminal pngcairo size 1200,900 enhanced font 'Verdana,10'
set output 'visualizacion_rbf_ns.png'

# --- Coordenadas del Obstáculo ---
obs_x1 = 0.5
obs_x2 = 1.0
obs_y  = 0.5

set multiplot layout 2,1 title "Navier-Stokes RBF-FD: Resultados de Simulación" font 'Verdana,14 bold'
set lmargin 10

# Dibujar el obstáculo en ambos paneles
set object 1 rect from obs_x1, 0.0 to obs_x2, obs_y front fillstyle solid 1.0 noborder fillcolor rgb "#222222"

# --- Configuración de Interpolación compatible ---
# Cambiamos 'RBF' por 'gauss' con una norma de 0.05 (ajustable según h)
set dgrid3d 100,100 gauss 0.05, 0.05

# ==============================================================================
# PANEL 1: Psi (Streamfunction)
# ==============================================================================
set title 'Función de Corriente ({/Symbol Y})' font 'Verdana,12'
set xlabel 'X'; set ylabel 'Y'
set xrange [0.0:2.0]; set yrange [-0.05:1.05]
set view map
set contour base
set cntrparam levels 25
unset key

set palette defined (0 "blue", 1 "green", 2 "yellow", 3 "red")
set cblabel "{/Symbol Y}"

splot 'FLUJO' using 1:2:3 with pm3d, \
      'FLUJO' using 1:2:3 with lines lc rgb "#555555" lw 0.5

# ==============================================================================
# PANEL 2: Omega (Vorticity)
# ==============================================================================
set title 'Vorticidad ({/Symbol w})' font 'Verdana,12'
set xlabel 'X'; set ylabel 'Y'
set xrange [0.0:2.0]; set yrange [-0.05:1.05]
set view map
set contour base
set cntrparam levels 30
unset clabel

# Paleta para vorticidad (Cyan-White-Red para mostrar rotación +/-)
set palette rgbformulae 33,13,10
set cblabel "{/Symbol w}"

splot 'VORT' using 1:2:3 with pm3d, \
      'VORT' using 1:2:3 with lines lc rgb "black" lw 0.5

unset multiplot
