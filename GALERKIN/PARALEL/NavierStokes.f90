! ==============================================================================
! SOLVER NAVIER-STOKES 2D - GALERKIN P1 (Ψ-ξ)
! ARQUITECTURA HPC MULTI-HILO (OPENMP 4.5+ / NUMA AWARE / MULTI-COLORING)
! ==============================================================================

module GlobalData
    use omp_lib
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    ! --- PANEL DE CONTROL FÍSICO Y NUMÉRICO ---
    real(dp) :: V0 = 1.0_dp
    real(dp) :: nu = 0.25_dp
    real(dp) :: omega_thom = 0.002_dp 
    real(dp) :: tol = 1.0d-6
    integer  :: max_picard = 25000
    integer  :: max_gs     = 15000

    ! --- ESTRUCTURAS DE DATOS TOPOLÓGICAS ---
    integer :: n_nodes, n_elements
    real(dp), allocatable :: x(:), y(:)
    
    ! [OPTIMIZACIÓN L1/TLB]: Formato Column-Major Estricto.
    integer, allocatable :: conn(:,:) 
    integer, allocatable :: node_type(:)

    ! --- MATRIZ DISPERSA (ELLPACK ALINEADA) ---
    integer, parameter :: max_adj = 50 
    real(dp), allocatable :: A_val(:,:)
    integer, allocatable :: A_col(:,:)
    integer, allocatable :: num_adj(:)
    real(dp), allocatable :: b_vec(:)

    ! --- VARIABLES DE ESTADO DEL FLUIDO ---
    real(dp), allocatable :: Psi(:), Xi(:), Psi_old(:)

    ! --- PRE-CÓMPUTO DE INVARIANTES TOPOLÓGICOS ---
    integer, allocatable :: closest_int_node(:)
    real(dp), allocatable :: inv_h_dist_sq(:) 

    ! [OpenMP] ESTRUCTURAS DE COLORACIÓN DE GRAFOS PARA GAUSS-SEIDEL
    integer, allocatable :: node_color(:)
    integer :: n_colors
end module GlobalData


module MeshOps
    use GlobalData
    implicit none
contains

    subroutine read_mesh()
        integer :: i, id_tmp
        
        open(10, file='nodos.txt', status='old')
        read(10, *) n_nodes
        allocate(x(n_nodes), y(n_nodes), node_type(n_nodes))
        do i = 1, n_nodes
            read(10, *) id_tmp, x(i), y(i), node_type(i)
        end do
        close(10)

        open(20, file='elementos.txt', status='old')
        read(20, *) n_elements
        ! Asignación Column-Major para el arreglo bidimensional
        allocate(conn(3, n_elements))
        do i = 1, n_elements
            read(20, *) id_tmp, conn(1,i), conn(2,i), conn(3,i)
        end do
        close(20)
    end subroutine read_mesh

    subroutine precompute_topology()
        integer :: i, j, c_int
        real(dp) :: min_dist_sq, dx, dy, d_sq

        allocate(closest_int_node(n_nodes), inv_h_dist_sq(n_nodes))
        closest_int_node = -1
        inv_h_dist_sq = 0.0_dp

        ! [INDEPENDIZACIÓN CÍCLICA]: Evaluación geométrica ejecutada antes del bucle transitorio.
        do i = 1, n_nodes
            if (node_type(i) == 1 .or. node_type(i) >= 5) then
                min_dist_sq = 1.0d30
                c_int = -1
                do j = 1, n_nodes
                    if (node_type(j) == 0) then
                        dx = x(i) - x(j)
                        dy = y(i) - y(j)
                        d_sq = dx*dx + dy*dy 
                        if (d_sq < min_dist_sq) then
                            min_dist_sq = d_sq
                            c_int = j
                        end if
                    end if
                end do
                closest_int_node(i) = c_int
                inv_h_dist_sq(i) = 1.0_dp / max(min_dist_sq, 1.0d-16)
            end if
        end do
    end subroutine precompute_topology

    ! [OpenMP] PASO 3: GRAPH COLORING (ALGORITMO AVARO)
    ! Se ejecuta serialmente. Asegura que ningún nodo adyacente posea el mismo color.
    subroutine compute_graph_coloring()
        integer :: i, k, j, c
        logical :: available_colors(max_adj)
        
        allocate(node_color(n_nodes))
        node_color = 0
        n_colors = 0
        
        print*, "[HPC Profiler] Generando Grafo de Colores para Gauss-Seidel Concurrente..."
        
        do i = 1, n_nodes
            available_colors = .true.
            
            ! Bloquear colores usados por los vecinos
            do k = 1, num_adj(i)
                j = A_col(k, i)
                if (j > 0 .and. j /= i .and. node_color(j) > 0) then
                    available_colors(node_color(j)) = .false.
                end if
            end do
            
            ! Asignar el primer color libre
            do c = 1, max_adj
                if (available_colors(c)) then
                    node_color(i) = c
                    n_colors = max(n_colors, c)
                    exit
                end if
            end do
        end do
        
        print*, "[HPC Profiler] Grafo completado. Colores independientes extraídos: ", n_colors
    end subroutine compute_graph_coloring

end module MeshOps


module SolverOps
    use GlobalData
    implicit none
contains

    ! [OpenMP] PASO 2: FASE SIMBÓLICA (Topología de la Matriz)
    ! Reemplaza el init_sparse original. Pre-calcula la conectividad serialmente
    ! para evitar reestructuraciones de memoria durante el multithreading.
    subroutine symbolic_assembly()
        integer :: e, i, j, n(3), k
        logical :: found
        
        num_adj = 0
        A_col = 0
        
        do e = 1, n_elements
            n(1)=conn(1,e); n(2)=conn(2,e); n(3)=conn(3,e)
            do i = 1, 3
                do j = 1, 3
                    found = .false.
                    do k = 1, num_adj(n(i))
                        if (A_col(k, n(i)) == n(j)) then
                            found = .true.
                            exit
                        end if
                    end do
                    if (.not. found) then
                        num_adj(n(i)) = num_adj(n(i)) + 1
                        A_col(num_adj(n(i)), n(i)) = n(j)
                    end if
                end do
            end do
        end do
    end subroutine symbolic_assembly

    ! [OpenMP] Reinicio numérico First-Touch para la fase de ensamblaje iterativo
    subroutine init_sparse_numeric()
        integer :: i, k
        !$omp parallel do schedule(static) private(i, k)
        do i = 1, n_nodes
            b_vec(i) = 0.0_dp
            do k = 1, num_adj(i)
                A_val(k, i) = 0.0_dp
            end do
        end do
        !$omp end parallel do
    end subroutine init_sparse_numeric

    ! [OpenMP] PASO 2: ATOMIC UPDATES
    ! Las adiciones superpuestas en la matriz desde múltiples elementos 
    ! se protegen en el hardware usando instrucciones atómicas ultrarrápidas.
    subroutine add_to_sparse_atomic(row, col, val)
        integer, intent(in) :: row, col
        real(dp), intent(in) :: val
        integer :: k
        do k = 1, num_adj(row)
            if (A_col(k, row) == col) then
                !$omp atomic update
                A_val(k, row) = A_val(k, row) + val
                return
            end if
        end do
    end subroutine add_to_sparse_atomic

    ! [OpenMP] PASO 1: PARALELIZACIÓN DE CONDICIÓN DE THOM
    subroutine update_thom_bc()
        integer :: i, c_int
        real(dp) :: Xi_new

        !$omp parallel do default(shared) private(i, c_int, Xi_new) schedule(static)
        do i = 1, n_nodes
            if (node_type(i) == 1 .or. node_type(i) >= 5) then
                c_int = closest_int_node(i)
                if (c_int /= -1) then
                    Xi_new = -2.0_dp * (Psi(c_int) - Psi(i)) * inv_h_dist_sq(i)
                    Xi(i) = (1.0_dp - omega_thom) * Xi(i) + omega_thom * Xi_new
                end if
            else if (node_type(i) == 2 .or. node_type(i) == 3) then
                Xi(i) = 0.0_dp
            end if
        end do
        !$omp end parallel do
    end subroutine update_thom_bc

    ! [OpenMP] PASO 2: ENSAMBLAJE NUMÉRICO CONCURRENTE (TRANSPORTE)
    subroutine assemble_vorticity_transport()
        integer :: e, tmp, i, j
        real(dp) :: x1, y1, x2, y2, x3, y3, detJ, Area
        real(dp) :: b(3), c(3), u, v, C_e, K_e, val
        real(dp) :: inv_2Area, inv_4Area, coeff_diff
        integer :: n(3)

        call init_sparse_numeric()

        !$omp parallel do default(shared) private(e, tmp, i, j, x1, y1, x2, y2, x3, y3, detJ, Area, b, c, u, v, C_e, K_e, val, inv_2Area, inv_4Area, coeff_diff, n) schedule(guided)
        do e = 1, n_elements
            n(1)=conn(1,e); n(2)=conn(2,e); n(3)=conn(3,e)
            x1=x(n(1)); y1=y(n(1)); x2=x(n(2)); y2=y(n(2)); x3=x(n(3)); y3=y(n(3))

            detJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            if (detJ < 0.0_dp) then
                tmp=n(2); n(2)=n(3); n(3)=tmp
                x2=x(n(2)); y2=y(n(2)); x3=x(n(3)); y3=y(n(3))
                detJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            end if
            
            Area = max(0.5_dp * detJ, 1.0d-12)

            b(1)=y2-y3; b(2)=y3-y1; b(3)=y1-y2
            c(1)=x3-x2; c(2)=x1-x3; c(3)=x2-x1

            inv_2Area = 0.5_dp / Area
            inv_4Area = 0.25_dp / Area
            coeff_diff = nu * inv_4Area

            u =  (Psi(n(1))*c(1) + Psi(n(2))*c(2) + Psi(n(3))*c(3)) * inv_2Area
            v = -(Psi(n(1))*b(1) + Psi(n(2))*b(2) + Psi(n(3))*b(3)) * inv_2Area

            do i = 1, 3
                do j = 1, 3
                    C_e = 0.1666666666666667_dp * (u*b(j) + v*c(j))
                    K_e = coeff_diff * (b(i)*b(j) + c(i)*c(j))
                    val = C_e + K_e
                    call add_to_sparse_atomic(n(i), n(j), val)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine assemble_vorticity_transport

    ! [OpenMP] PASO 2: ENSAMBLAJE NUMÉRICO CONCURRENTE (POISSON)
    subroutine assemble_poisson()
        integer :: e, tmp, i, j
        real(dp) :: x1, y1, x2, y2, x3, y3, detJ, Area
        real(dp) :: b(3), c(3), K_e, f_e(3)
        real(dp) :: inv_4Area, area_div_12
        integer :: n(3)

        call init_sparse_numeric()

        !$omp parallel do default(shared) private(e, tmp, i, j, x1, y1, x2, y2, x3, y3, detJ, Area, b, c, K_e, f_e, inv_4Area, area_div_12, n) schedule(guided)
        do e = 1, n_elements
            n(1)=conn(1,e); n(2)=conn(2,e); n(3)=conn(3,e)
            x1=x(n(1)); y1=y(n(1)); x2=x(n(2)); y2=y(n(2)); x3=x(n(3)); y3=y(n(3))

            detJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            if (detJ < 0.0_dp) then
                tmp=n(2); n(2)=n(3); n(3)=tmp
                x2=x(n(2)); y2=y(n(2)); x3=x(n(3)); y3=y(n(3))
                detJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            end if
            
            Area = max(0.5_dp * detJ, 1.0d-12)

            b(1)=y2-y3; b(2)=y3-y1; b(3)=y1-y2
            c(1)=x3-x2; c(2)=x1-x3; c(3)=x2-x1

            inv_4Area = 0.25_dp / Area
            area_div_12 = Area * 0.0833333333333333_dp

            f_e(1) = area_div_12 * (2.0_dp*Xi(n(1)) + Xi(n(2)) + Xi(n(3)))
            f_e(2) = area_div_12 * (Xi(n(1)) + 2.0_dp*Xi(n(2)) + Xi(n(3)))
            f_e(3) = area_div_12 * (Xi(n(1)) + Xi(n(2)) + 2.0_dp*Xi(n(3)))

            do i = 1, 3
                ! [OpenMP] Modificación del vector fuerza con política atómica
                !$omp atomic update
                b_vec(n(i)) = b_vec(n(i)) + f_e(i)
                do j = 1, 3
                    K_e = inv_4Area * (b(i)*b(j) + c(i)*c(j))
                    call add_to_sparse_atomic(n(i), n(j), K_e)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine assemble_poisson

    ! [OpenMP] PASO 3: GAUSS-SEIDEL BASADO EN MULTI-COLORING
    subroutine solve_gs(var)
        real(dp), intent(inout) :: var(:)
        integer :: it, c, i, k, j
        real(dp) :: sum_j, diag_val
        
        do it = 1, max_gs
            ! Iteración secuencial exterior (barrido de colores)
            do c = 1, n_colors
                
                ! Iteración totalmente paralela y segura. Los nodos del mismo color NO son vecinos.
                !$omp parallel do default(shared) private(i, k, j, sum_j, diag_val) schedule(static)
                do i = 1, n_nodes
                    if (node_color(i) == c) then
                        if (node_type(i) == 0 .or. node_type(i) == 4) then
                            sum_j = 0.0_dp
                            diag_val = 1.0_dp
                            
                            do k = 1, num_adj(i)
                                j = A_col(k, i)
                                if (j == i) then
                                    diag_val = A_val(k, i)
                                else if (j > 0) then
                                    sum_j = sum_j + A_val(k, i) * var(j)
                                end if
                            end do
                            
                            if (abs(diag_val) > 1.0d-12) then
                                var(i) = (b_vec(i) - sum_j) / diag_val
                            end if
                        end if
                    end if
                end do
                !$omp end parallel do

            end do
        end do
    end subroutine solve_gs

end module SolverOps


program NavierStokes_Galerkin
    use GlobalData
    use MeshOps
    use SolverOps
    implicit none
    
    integer :: iter, i
    real(dp) :: error_sum, diff, error_global

    call read_mesh()
    call precompute_topology() 
    
    allocate(Psi(n_nodes), Xi(n_nodes), Psi_old(n_nodes))
    allocate(A_val(max_adj, n_nodes), A_col(max_adj, n_nodes))
    allocate(num_adj(n_nodes), b_vec(n_nodes))
    
    ! [OpenMP] PASO 1: POLÍTICA FIRST-TOUCH NUMA
    ! Cada núcleo asignará físicamente las páginas RAM que luego computará.
    !$omp parallel do schedule(static) private(i)
    do i = 1, n_nodes
        Psi(i) = 0.0_dp
        Xi(i) = 0.0_dp
        Psi_old(i) = 0.0_dp
        
        if (node_type(i) == 3) Psi(i) = V0 * y(i)      
        if (node_type(i) == 2) Psi(i) = V0 * 1.0_dp    
        if (node_type(i) == 1 .or. node_type(i) >= 5) Psi(i) = 0.0_dp 
    end do
    !$omp end parallel do

    ! Fases Pre-Transitorias y Estructurales
    call symbolic_assembly()
    call compute_graph_coloring()

    print*, "=============================================================="
    print*, " INICIANDO N-S HPC (OPENMP MULTI-HILO / NUMA AWARE)"
    print*, " Nodos:", n_nodes, " | Hilos Lógicos Asignados: ", omp_get_max_threads()
    print*, "=============================================================="

    do iter = 1, max_picard
        !$omp parallel do schedule(static) private(i)
        do i = 1, n_nodes
            Psi_old(i) = Psi(i)
        end do
        !$omp end parallel do

        call update_thom_bc()
        call assemble_vorticity_transport()
        call solve_gs(Xi)
        call assemble_poisson()
        call solve_gs(Psi)

        ! [OpenMP] REDUCTION SEGURA PARA EL CÁLCULO DEL ERROR GLOBAL
        error_sum = 0.0_dp
        !$omp parallel do default(shared) private(i, diff) reduction(+:error_sum) schedule(static)
        do i = 1, n_nodes
            diff = Psi(i) - Psi_old(i)
            error_sum = error_sum + (diff * diff)
        end do
        !$omp end parallel do
        
        error_global = sqrt(error_sum / n_nodes)
        
        if (mod(iter, 50) == 0) then
            print*, "Ciclo Integrativo:", iter, "| Residual Global:", error_global
        end if
        
        if (error_global > 1.0d3 .or. error_global /= error_global) then
            print*, "[ALERTA CRÍTICA] Divergencia geométrica detectada en ciclo", iter
            exit
        end if
        
        if (error_global < tol) then
            print*, " >> CONVERGENCIA ANALÍTICA ALCANZADA EN CICLO:", iter
            exit
        end if
    end do

    open(30, file='FLUJO')
    open(40, file='VORT')
    do i = 1, n_nodes
        write(30,*) i, x(i), y(i), Psi(i)
        write(40,*) i, x(i), y(i), Xi(i)
    end do
    close(30)
    close(40)
    
    print*, "Descarga de campo tensorial finalizada (FLUJO, VORT)."

end program NavierStokes_Galerkin
