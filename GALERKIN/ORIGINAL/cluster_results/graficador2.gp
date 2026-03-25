# ==============================================================================
# SCRIPT DE VISUALIZACIÓN AL ESTILO CLÁSICO CFD
# ==============================================================================

# Configuración de salida
set terminal pngcairo size 900,400 enhanced font 'Verdana,12'
set xlabel "X"
set ylabel "Y"

# Paleta de colores cálida (Jet), clásica en mecánica de fluidos
set palette rgbformulae 22,13,-31

# Definición del obstáculo sólido (Coordenadas físicas reales)
set object 1 rect from 1.0, 0.0 to 1.5, 0.5 fc rgb "black" fillstyle solid 1.0 front

# Configuración de visualización
set view map
# Interpolación para suavizar los colores y dgrid3d para leer la malla
set pm3d map interpolate 2,2
set dgrid3d 150,30 gauss 0.05
unset key

# ==============================================================================
# --- GRÁFICA 1: FUNCIÓN DE CORRIENTE ---
# ==============================================================================
set output 'flujo_con_obstaculo.png'
set title "Función de Corriente con Obstáculo Sólido"

# Activar contornos con líneas finas negras
set contour base
set cntrparam levels 25

# Graficamos el mapa de calor (pm3d) y encima las líneas negras (lines)
splot "FLUJO" using 2:3:4 with pm3d, \
      "FLUJO" using 2:3:4 nosurface with lines lc rgb "black" lw 0.6

# ==============================================================================
# --- GRÁFICA 2: VORTICIDAD ---
# ==============================================================================
set output 'vorticidad_con_bloque.png'
set title "Distribución de Vorticidad"

# Para la vorticidad suele ser mejor menos líneas (o ninguna). 
set cntrparam levels 15

splot "VORT" using 2:3:4 with pm3d, \
      "VORT" using 2:3:4 nosurface with lines lc rgb "black" lw 0.4

