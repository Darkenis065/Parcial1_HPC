import numpy as np

# --- PARÁMETROS DEL DOMINIO ---
L, H = 5.0, 1.0        # Largo y ancho del canal
obs_x, obs_y = 1.0, 0.0 # Esquina inferior izquierda del obstáculo
obs_w, obs_h = 0.5, 0.5 # Dimensiones del obstáculo (rectángulo)
nx, ny = 151, 31         # Resolución de la malla

dx = L / (nx - 1)
dy = H / (ny - 1)

nodes = []
# --- GENERACIÓN DE NODOS Y CLASIFICACIÓN ---
node_id = 1
for j in range(ny):
    for i in range(nx):
        curr_x = i * dx
        curr_y = j * dy
        
        # Clasificación por defecto: Interno
        btype = 0
        
        # 1. Identificar si el nodo está DENTRO del obstáculo (se omite o se marca)
        # Para este ejemplo, los marcaremos para que el solver sepa ignorarlos o usarlos como pared.
        is_inside = (obs_x < curr_x < obs_x + obs_w) and (obs_y < curr_y < obs_y + obs_h)
        
        if is_inside:
            continue # No generamos nodos dentro del sólido para mallas no estructuradas simples
            
        # 2. Fronteras del Canal
        if j == 0: btype = 1                     # Pared Inferior
        elif j == ny-1: btype = 2                # Pared Superior
        elif i == 0: btype = 3                   # Entrada
        elif i == nx-1: btype = 4                # Salida
        
        # 3. Fronteras del Obstáculo (Paredes internas)
        # Margen de tolerancia para detectar la superficie del obstáculo
        tol = 1e-7
        if abs(curr_x - obs_x) < tol and obs_y <= curr_y <= obs_y + obs_h: btype = 5 # Obs Izq
        elif abs(curr_x - (obs_x + obs_w)) < tol and obs_y <= curr_y <= obs_y + obs_h: btype = 6 # Obs Der
        elif abs(curr_y - obs_y) < tol and obs_x <= curr_x <= obs_x + obs_w: btype = 7 # Obs Inf
        elif abs(curr_y - (obs_y + obs_h)) < tol and obs_x <= curr_x <= obs_x + obs_w: btype = 8 # Obs Sup

        nodes.append((node_id, curr_x, curr_y, btype))
        node_id += 1

# --- ESCRITURA DEL ARCHIVO ---
with open("nodos.txt", "w") as f:
    f.write(f"{len(nodes)}\n") # Cabecera con el total de nodos
    for n in nodes:
        f.write(f"{n[0]} {n[1]:.6f} {n[2]:.6f} {n[3]}\n")

print(f"Archivo nodos.txt generado con {len(nodes)} nodos.")
