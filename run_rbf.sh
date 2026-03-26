#!/bin/bash
set -e

mkdir -p data

run_rbf() {
    echo "--- Benchmarking RBF-FD ---"

    # Original
    cd RBF_FD/ORIGINAL
    gfortran -O3 types_mod.mod sparse_mod.mod kdtree_mod.mod rbf_weights_mod.mod solver_ns_mod.mod NavierStokes.f90 -llapack -lblas -o sim || gfortran -O3 NavierStokes.f90 -llapack -lblas -o sim || echo "Failed compilation RBF_FD Original"
    if [ -f sim ]; then T_ORIG=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T_ORIG="0"; fi
    cd ../..

    # Optimized
    cd RBF_FD/OPTIMIZED
    gfortran -O3 -march=native -ffast-math -funroll-loops types_mod.mod sparse_mod.mod kdtree_mod.mod rbf_weights_mod.mod solver_ns_mod.mod NavierStokes.f90 -llapack -lblas -o sim || gfortran -O3 -march=native -ffast-math -funroll-loops NavierStokes.f90 -llapack -lblas -o sim || echo "Failed compilation RBF_FD Optimized"
    if [ -f sim ]; then T_OPT=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T_OPT="0"; fi
    cd ../..

    # Parallel
    cd RBF_FD/PARALEL
    gfortran -O3 -march=native -fopenmp -ffast-math -funroll-loops types_mod.mod sparse_mod.mod kdtree_mod.mod rbf_weights_mod.mod solver_ns_mod.mod NavierStokes.f90 -llapack -lblas -o sim || gfortran -O3 -march=native -fopenmp -ffast-math -funroll-loops NavierStokes.f90 -llapack -lblas -o sim || echo "Failed compilation RBF_FD Parallel"
    export OMP_NUM_THREADS=20
    if [ -f sim ]; then T_PAR=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T_PAR="0"; fi
    cd ../..

    echo "Original: $T_ORIG, Optimized: $T_OPT, Parallel(20): $T_PAR"
    echo -e "Original\t$T_ORIG\nOptimized\t$T_OPT\nParallel\t$T_PAR" > data/RBF-FD_comparison.dat

    # Parallel scaling
    echo -e "Threads\tTime" > data/RBF-FD_parallel_scaling.dat
    for t in {5..40..5}; do
        export OMP_NUM_THREADS=$t
        cd RBF_FD/PARALEL
        if [ -f sim ]; then T=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T="0"; fi
        cd ../..
        echo -e "$t\t$T" >> data/RBF-FD_parallel_scaling.dat
    done
}

run_rbf

echo "All benchmarks completed."
