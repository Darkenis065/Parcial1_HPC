#!/bin/bash
set -e

mkdir -p bilder

echo "--- Generating Simulation Visualizations ---"

# FINITE DIFFERENCE
echo "Plotting Finite Difference (Original)..."
python3 FINITEDIFF/ORIGINAL/plot_cfd.py FINITEDIFF/ORIGINAL/FLUJO FINITEDIFF/ORIGINAL/VORT bilder/FINITEDIFF_Original_Simulation.png "Diferencias Finitas (Original)" || echo "No data for FINITEDIFF Original"

echo "Plotting Finite Difference (Optimized)..."
python3 FINITEDIFF/OPTIMIZED/plot_cfd.py FINITEDIFF/OPTIMIZED/FLUJO FINITEDIFF/OPTIMIZED/VORT bilder/FINITEDIFF_Optimized_Simulation.png "Diferencias Finitas (Optimized)" || echo "No data for FINITEDIFF Optimized"

echo "Plotting Finite Difference (Parallel)..."
python3 FINITEDIFF/PARALEL/plot_cfd.py FINITEDIFF/PARALEL/FLUJO FINITEDIFF/PARALEL/VORT bilder/FINITEDIFF_Parallel_Simulation.png "Diferencias Finitas (Parallel)" || echo "No data for FINITEDIFF Parallel"


# GALERKIN
echo "Plotting Galerkin (Original)..."
python3 GALERKIN/ORIGINAL/plot_cfd.py GALERKIN/ORIGINAL/FLUJO GALERKIN/ORIGINAL/VORT bilder/Galerkin_Original_Simulation.png "Galerkin (Original)" || echo "No data for Galerkin Original"

echo "Plotting Galerkin (Optimized)..."
python3 GALERKIN/OPTIMIZED/plot_cfd.py GALERKIN/OPTIMIZED/FLUJO GALERKIN/OPTIMIZED/VORT bilder/Galerkin_Optimized_Simulation.png "Galerkin (Optimized)" || echo "No data for Galerkin Optimized"

echo "Plotting Galerkin (Parallel)..."
python3 GALERKIN/PARALEL/plot_cfd.py GALERKIN/PARALEL/FLUJO GALERKIN/PARALEL/VORT bilder/Galerkin_Parallel_Simulation.png "Galerkin (Parallel)" || echo "No data for Galerkin Parallel"


# RBF-FD
echo "Plotting RBF-FD (Original)..."
python3 RBF_FD/ORIGINAL/plot_cfd.py RBF_FD/ORIGINAL/FLUJO RBF_FD/ORIGINAL/VORT bilder/RBF-FD_Original_Simulation.png "RBF-FD (Original)" || echo "No data for RBF-FD Original"

echo "Plotting RBF-FD (Optimized)..."
python3 RBF_FD/OPTIMIZED/plot_cfd.py RBF_FD/OPTIMIZED/FLUJO RBF_FD/OPTIMIZED/VORT bilder/RBF-FD_Optimized_Simulation.png "RBF-FD (Optimized)" || echo "No data for RBF-FD Optimized"

echo "Plotting RBF-FD (Parallel)..."
python3 RBF_FD/PARALEL/plot_cfd.py RBF_FD/PARALEL/FLUJO RBF_FD/PARALEL/VORT bilder/RBF-FD_Parallel_Simulation.png "RBF-FD (Parallel)" || echo "No data for RBF-FD Parallel"

echo "--- Generating Performance Benchmarks ---"
python3 plot_benchmarks.py

echo "All plots generated successfully in 'bilder' folder."
