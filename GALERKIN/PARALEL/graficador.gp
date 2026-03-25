# ==============================================================================
# SCRIPT DE VISUALIZACIÓN CFD - NAVIER-STOKES (PSI-XI)
# ==============================================================================

set terminal pngcairo size 1200,900 enhanced font 'Verdana,11'
set output 'resultados_simulacion_final.png'

set multiplot layout 2,1 title "Simulación Navier-Stokes 2D (Galerkin P1) - Re = 2" font ",14"

set view map
set size ratio -1          
set xrange [0:5]
set yrange [0:1]
set palette rgbformulae 22,13,-31  

# Dibujar el obstáculo
set object 1 rect from 1.0, 0.0 to 1.5, 0.5 fc rgb "black" fillstyle solid 1.0 front

# ==============================================================================
# GRÁFICA 1: FUNCIÓN DE CORRIENTE (CON LÍNEAS DE FLUJO)
# ==============================================================================
set title "Función de Corriente (Psi) y Líneas de Flujo"
set cblabel "Valor de Psi"

# Radio de Gauss ajustado a 0.05 para que alcance a "ver" a los vecinos
set dgrid3d 150,30 gauss 0.05
set pm3d map
set contour base
set cntrparam levels 25    

splot "FLUJO" u 2:3:4 with pm3d notitle, \
      "FLUJO" u 2:3:4 nosurface with lines lc rgb "white" lw 1.2 notitle

# ==============================================================================
# GRÁFICA 2: VORTICIDAD (MAPA DE CALOR CONTINUO)
# ==============================================================================
unset contour

set title "Distribución de Vorticidad (Xi)"
set cblabel "Valor de Xi"

# Radio ajustado a 0.05
set dgrid3d 150,30 gauss 0.05
splot "VORT" u 2:3:4 with pm3d notitle

unset multiplot

