#!/bin/bash
set -e

mkdir -p data

run_finitediff() {
    echo "--- Benchmarking Finite Difference ---"

    # Original
    cd FINITEDIFF/ORIGINAL
    gfortran -O3 Navier_Stokes.f -o sim
    TIMEFORMAT=%R
    T_ORIG=$( { time ./sim > /dev/null 2>&1; } 2>&1 )
    cd ../..

    # Optimized
    cd FINITEDIFF/OPTIMIZED
    gfortran -O3 -march=native -ffast-math -funroll-loops NavierStokes.f -o sim
    T_OPT=$( { time ./sim > /dev/null 2>&1; } 2>&1 )
    cd ../..

    # Parallel
    cd FINITEDIFF/PARALEL
    gfortran -O3 -march=native -fopenmp -ffast-math -funroll-loops NavierStokes.f -o sim
    export OMP_NUM_THREADS=20
    T_PAR=$( { time ./sim > /dev/null 2>&1; } 2>&1 )
    cd ../..

    echo "Original: $T_ORIG, Optimized: $T_OPT, Parallel(20): $T_PAR"
    echo -e "Original\t$T_ORIG\nOptimized\t$T_OPT\nParallel\t$T_PAR" > data/Finite_Difference_comparison.dat

    # Parallel scaling
    echo -e "Threads\tTime" > data/Finite_Difference_parallel_scaling.dat
    for t in {5..40..5}; do
        export OMP_NUM_THREADS=$t
        cd FINITEDIFF/PARALEL
        T=$( { time ./sim > /dev/null 2>&1; } 2>&1 )
        cd ../..
        echo -e "$t\t$T" >> data/Finite_Difference_parallel_scaling.dat
    done
}

run_rbf() {
    echo "--- Benchmarking RBF-FD ---"

    # Original
    cd RBF_FD/ORIGINAL
    gfortran -O3 NavierStokes.f90 -llapack -lblas -o sim || echo "Failed compilation RBF_FD Original"
    if [ -f sim ]; then T_ORIG=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T_ORIG="0"; fi
    cd ../..

    # Optimized
    cd RBF_FD/OPTIMIZED
    gfortran -O3 -march=native -ffast-math -funroll-loops NavierStokes.f90 -llapack -lblas -o sim || echo "Failed compilation RBF_FD Optimized"
    if [ -f sim ]; then T_OPT=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T_OPT="0"; fi
    cd ../..

    # Parallel
    cd RBF_FD/PARALEL
    gfortran -O3 -march=native -fopenmp -ffast-math -funroll-loops NavierStokes.f90 -llapack -lblas -o sim || echo "Failed compilation RBF_FD Parallel"
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

run_galerkin() {
    echo "--- Benchmarking Galerkin ---"

    # Need to generate mesh first for Galerkin!
    cd GALERKIN/ORIGINAL
    python3 MeshGenerator.py > /dev/null 2>&1 || true
    python3 ElementGenerator.py > /dev/null 2>&1 || true
    gfortran -O3 NavierStokes.f90 -o sim || echo "Failed compilation Galerkin Original"
    if [ -f sim ]; then T_ORIG=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T_ORIG="0"; fi
    echo "Galerkin Original finished"
    cd ../..

    cd GALERKIN/OPTIMIZED
    python3 MeshGenerator.py > /dev/null 2>&1 || true
    python3 ElementGenerator.py > /dev/null 2>&1 || true
    gfortran -O3 -march=native -ffast-math -funroll-loops NavierStokes.f90 -o sim || echo "Failed compilation Galerkin Optimized"
    if [ -f sim ]; then T_OPT=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T_OPT="0"; fi
    echo "Galerkin Optimized finished"
    cd ../..

    cd GALERKIN/PARALEL
    cp ../ORIGINAL/nodos.txt . || true
    cp ../ORIGINAL/elementos.txt . || true
    gfortran -O3 -march=native -fopenmp -ffast-math -funroll-loops NavierStokes.f90 -o sim || echo "Failed compilation Galerkin Parallel"
    export OMP_NUM_THREADS=20
    if [ -f sim ]; then T_PAR=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T_PAR="0"; fi
    cd ../..

    echo "Original: $T_ORIG, Optimized: $T_OPT, Parallel(20): $T_PAR"
    echo -e "Original\t$T_ORIG\nOptimized\t$T_OPT\nParallel\t$T_PAR" > data/Galerkin_comparison.dat

    # Parallel scaling
    echo -e "Threads\tTime" > data/Galerkin_parallel_scaling.dat
    for t in {5..40..5}; do
        export OMP_NUM_THREADS=$t
        cd GALERKIN/PARALEL
        if [ -f sim ]; then T=$( { time ./sim > /dev/null 2>&1; } 2>&1 ); else T="0"; fi
        cd ../..
        echo -e "$t\t$T" >> data/Galerkin_parallel_scaling.dat
        echo "Galerkin Parallel Thread $t finished"
    done
}

run_finitediff
run_galerkin
run_rbf

echo "All benchmarks completed."
