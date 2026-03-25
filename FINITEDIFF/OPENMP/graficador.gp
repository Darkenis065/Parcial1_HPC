# Configuración de salida
set terminal pngcairo size 800,700 enhanced font 'Verdana,12'
set xlabel "X (I)"
set ylabel "Y (J)"

# Paleta de colores cálida del libro
set palette rgbformulae 22,13,10

# DEFINICIÓN DEL OBSTÁCULO (Basado en I1, I2 y J1 del código)
# Se dibuja un rectángulo sólido para representar el cuerpo físico
set object 1 rectangle from 46,1 to 80,46 fc rgb "black" fillstyle solid 1.0 front

# Configuración de visualización
set view map
set pm3d at b interpolate 2,2
set contour base
set cntrparam level 25
unset key

# --- Graficar el FLUJO ---
set output 'flujo_con_obstaculo.png'
set title "Función de Corriente con Obstáculo Sólido"
splot "FLUJO" using 1:2:3 with pm3d, \
      "FLUJO" using 1:2:3 with lines nosurface lc rgb "black" lw 0.5

# --- Graficar la VORTICIDAD ---
set output 'vorticidad_con_bloque.png'
set title "Distribución de Vorticidad"
set cntrparam level 15
splot "VORT" using 1:2:3 with pm3d, \
      "VORT" using 1:2:3 with lines nosurface lc rgb "black" lw 0.5
