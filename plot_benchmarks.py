import matplotlib.pyplot as plt
import os
import sys

def parse_comparison_file(filepath):
    implementations = []
    times = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    imp, t = parts[0], parts[1]
                    try:
                        t = float(t)
                        implementations.append(imp)
                        times.append(t)
                    except ValueError:
                        pass
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return implementations, times

def parse_scaling_file(filepath):
    threads = []
    times = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0] != "Threads":
                    try:
                        th = int(parts[0])
                        t = float(parts[1])
                        threads.append(th)
                        times.append(t)
                    except ValueError:
                        pass
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return threads, times

os.makedirs('bilder', exist_ok=True)
methods = ['Finite_Difference', 'Galerkin', 'RBF-FD']
method_titles = ['Diferencias Finitas', 'Galerkin', 'RBF-FD']

for method, title in zip(methods, method_titles):
    comp_file = f'data/{method}_comparison.dat'
    scal_file = f'data/{method}_parallel_scaling.dat'

    # 1. Comparison Bar Chart
    imps, comp_times = parse_comparison_file(comp_file)
    if imps:
        plt.figure(figsize=(8, 6))
        bars = plt.bar(imps, comp_times, color=['#440154', '#21918c', '#fde725'])
        plt.title(f'{title} - Comparación de Tiempos de Ejecución (Total)')
        plt.ylabel('Tiempo (segundos)')
        plt.xlabel('Implementación')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, yval, round(yval, 3), ha='center', va='bottom')
        plt.tight_layout()
        plt.savefig(f'bilder/Comparacion_{method}_Total.png', dpi=300)
        plt.close()

    # 2. Parallel Scaling Line Chart
    threads, scal_times = parse_scaling_file(scal_file)
    if threads:
        plt.figure(figsize=(8, 6))
        plt.plot(threads, scal_times, marker='o', linestyle='-', color='#e16462', linewidth=2)
        plt.title(f'{title} - Escalabilidad Paralela (Parallel)')
        plt.ylabel('Tiempo (segundos)')
        plt.xlabel('Número de Hilos (OpenMP)')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.xticks(threads)
        for i, t in enumerate(scal_times):
             plt.text(threads[i], t, f'{t:.2f}', ha='center', va='bottom', fontsize=9)
        plt.tight_layout()
        plt.savefig(f'bilder/Comparacion_{method}_Parallel.png', dpi=300)
        plt.close()

print("Benchmark plots generated successfully in 'bilder' folder.")
