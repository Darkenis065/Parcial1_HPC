import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 5:
    print("Usage: python plot_cfd.py <FLUJO_file> <VORT_file> <Output_file.png> <Title>")
    sys.exit(1)

flujo_file = sys.argv[1]
vort_file = sys.argv[2]
output_file = sys.argv[3]
plot_title = sys.argv[4]

def load_data(filename):
    try:
        data = np.loadtxt(filename)
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        sys.exit(1)

    if len(data) == 0:
        return np.array([]), np.array([]), np.array([]), False

    if data.shape[1] >= 4:
        # Format might be ID X Y Value
        x = data[:, 1]
        y = data[:, 2]
        z = data[:, 3]
    else:
        # Standard format X Y Value
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]

    # Find unique grid points
    ux = np.unique(x)
    uy = np.unique(y)

    # If not a perfect grid, use tricontourf
    is_grid = (len(ux) * len(uy) == len(z)) and len(x) > 0

    if is_grid:
        X, Y = np.meshgrid(ux, uy)
        Z = np.zeros_like(X)
        for i in range(len(x)):
            idx_x = np.where(ux == x[i])[0][0]
            idx_y = np.where(uy == y[i])[0][0]
            Z[idx_y, idx_x] = z[i]
        return X, Y, Z, is_grid
    else:
        return x, y, z, is_grid

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
fig.suptitle(plot_title, fontsize=16)

# Try to determine obstacle bounds based on plot_title
if "Diferencias Finitas" in plot_title or "FINITEDIFF" in plot_title:
    obs_x_start = 46; obs_x_end = 80
    obs_y_start = 1; obs_y_end = 46
elif "RBF-FD" in plot_title or "RBF_FD" in plot_title:
    obs_x_start = 0.5; obs_x_end = 1.0
    obs_y_start = 0.0; obs_y_end = 0.5
else:
    obs_x_start = 1.0; obs_x_end = 1.5
    obs_y_start = 0.0; obs_y_end = 0.5

# Load Flujo
x1, y1, z1, is_grid1 = load_data(flujo_file)
if len(x1) > 0:
    if is_grid1:
        cf1 = ax1.contourf(x1, y1, z1, 25, cmap='magma')
        ax1.contour(x1, y1, z1, 25, colors='black', linewidths=0.5, alpha=0.5)
    else:
        cf1 = ax1.tricontourf(x1, y1, z1, 25, cmap='magma')
        ax1.tricontour(x1, y1, z1, 25, colors='black', linewidths=0.5, alpha=0.5)

    rect1 = plt.Rectangle((obs_x_start, obs_y_start), obs_x_end - obs_x_start, obs_y_end - obs_y_start,
                        fc='black', ec='black', zorder=10)
    ax1.add_patch(rect1)
    fig.colorbar(cf1, ax=ax1)

ax1.set_title('Flujo (Streamfunction)')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

# Load Vort
x2, y2, z2, is_grid2 = load_data(vort_file)
if len(x2) > 0:
    if is_grid2:
        cf2 = ax2.contourf(x2, y2, z2, 25, cmap='magma')
    else:
        cf2 = ax2.tricontourf(x2, y2, z2, 25, cmap='magma')

    rect2 = plt.Rectangle((obs_x_start, obs_y_start), obs_x_end - obs_x_start, obs_y_end - obs_y_start,
                        fc='black', ec='black', zorder=10)
    ax2.add_patch(rect2)
    fig.colorbar(cf2, ax=ax2)

ax2.set_title('Vorticidad (Vorticity)')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')

plt.tight_layout()
plt.savefig(output_file, dpi=300)
print(f"Plot saved to {output_file}")
