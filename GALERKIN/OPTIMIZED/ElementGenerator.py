import numpy as np
from scipy.spatial import Delaunay

# --- MISMOS PARÁMETROS DEL SCRIPT DE NODOS ---
L, H = 5.0, 1.0
obs_x, obs_y = 1.0, 0.0
obs_w, obs_h = 0.5, 0.5
tol = 1e-5

# 1. Cargar los nodos generados para asegurar consistencia
nodes = []
with open("nodos.txt", "r") as f:
    n_nodes = int(f.readline())
    for line in f:
        nodes.append(list(map(float, line.split())))

nodes = np.array(nodes) # [ID, x, y, type]
points = nodes[:, 1:3]  # Solo coordenadas X, Y

# 2. Generar Triangulación de Delaunay
tri = Delaunay(points)

# 3. Filtrar elementos que están "dentro" del obstáculo
# Calculamos el baricentro de cada triángulo para decidir si se mantiene
final_elements = []
for i, t in enumerate(tri.simplices):
    # t contiene los índices de los 3 nodos (0-based)
    p1, p2, p3 = points[t]
    
    # Baricentro (centro del triángulo)
    cx = (p1[0] + p2[0] + p3[0]) / 3.0
    cy = (p1[1] + p2[1] + p3[1]) / 3.0
    
    # Condición: El centro no debe estar dentro del rectángulo del obstáculo
    # Usamos la tolerancia para no borrar triángulos que tocan la pared
    is_inside_obs = (obs_x + tol < cx < obs_x + obs_w - tol) and \
                    (obs_y + tol < cy < obs_y + obs_h - tol)
    
    if not is_inside_obs:
        # Guardamos: [ID_Elemento, Nodo1, Nodo2, Nodo3]
        # Sumamos +1 a los índices porque Fortran suele usar base 1
        n1, n2, n3 = int(nodes[t[0], 0]), int(nodes[t[1], 0]), int(nodes[t[2], 0])
        final_elements.append((len(final_elements) + 1, n1, n2, n3))

# 4. Escritura del archivo elementos.txt
with open("elementos.txt", "w") as f:
    f.write(f"{len(final_elements)}\n")
    for e in final_elements:
        f.write(f"{e[0]} {e[1]} {e[2]} {e[3]}\n")

print(f"Archivo elementos.txt generado con {len(final_elements)} triángulos.")
