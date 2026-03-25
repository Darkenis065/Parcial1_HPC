! ==============================================================================
! SOLUCIONADOR RBF-FD NAVIER-STOKES 2D (FLUJO-VORTICIDAD)
! ESTABILIZACION: RK4, BiCGSTAB Matrix-Free, Hiperviscosidad, Ghost Nodes Nativos
! ==============================================================================
MODULE types_mod
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
END MODULE types_mod

! ==============================================================================
! MODULO 1: KD-TREE (TOPOLOGIA NODAL)
! ==============================================================================
MODULE kdtree_mod
    USE types_mod
    IMPLICIT NONE
    TYPE :: kd_node
        INTEGER :: id
        REAL(DP), DIMENSION(2) :: coord
        TYPE(kd_node), POINTER :: left => NULL(), right => NULL()
    END TYPE kd_node

    CONTAINS

    RECURSIVE SUBROUTINE build_kdtree(ids, coords, n, depth, root)
        INTEGER, DIMENSION(:), INTENT(INOUT) :: ids
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: coords
        INTEGER, INTENT(IN) :: n, depth
        TYPE(kd_node), POINTER, INTENT(OUT) :: root
        INTEGER :: axis, median_idx, i, j, temp_id
        REAL(DP) :: temp_coord(2)

        IF (n <= 0) THEN
            root => NULL(); RETURN
        END IF

        ALLOCATE(root)
        axis = MOD(depth, 2) + 1

        DO i = 1, n - 1
            DO j = i + 1, n
                IF (coords(axis, j) < coords(axis, i)) THEN
                    temp_coord = coords(:, i); coords(:, i) = coords(:, j); coords(:, j) = temp_coord
                    temp_id = ids(i); ids(i) = ids(j); ids(j) = temp_id
                END IF
            END DO
        END DO
        median_idx = n / 2 + 1
        root%id = ids(median_idx)
        root%coord = coords(:, median_idx)
        CALL build_kdtree(ids(1:median_idx-1), coords(:, 1:median_idx-1), median_idx-1, depth+1, root%left)
        CALL build_kdtree(ids(median_idx+1:n), coords(:, median_idx+1:n), n-median_idx, depth+1, root%right)
    END SUBROUTINE build_kdtree

    RECURSIVE SUBROUTINE search_knn(node, target, k, k_ids, k_dists, depth)
        TYPE(kd_node), POINTER, INTENT(IN) :: node
        REAL(DP), DIMENSION(2), INTENT(IN) :: target
        INTEGER, INTENT(IN) :: k
        INTEGER, DIMENSION(k), INTENT(INOUT) :: k_ids
        REAL(DP), DIMENSION(k), INTENT(INOUT) :: k_dists
        INTEGER, INTENT(IN) :: depth
        REAL(DP) :: dist, dx, dy
        INTEGER :: axis, i, pos
        TYPE(kd_node), POINTER :: next_branch, opposite_branch

        IF (.NOT. ASSOCIATED(node)) RETURN
        dx = node%coord(1) - target(1); dy = node%coord(2) - target(2)
        dist = SQRT(dx**2 + dy**2)

        IF (dist < k_dists(k)) THEN
            pos = k
            DO i = 1, k
                IF (dist < k_dists(i)) THEN
                    pos = i; EXIT
                END IF
            END DO
            DO i = k, pos + 1, -1
                k_dists(i) = k_dists(i-1); k_ids(i) = k_ids(i-1)
            END DO
            k_dists(pos) = dist; k_ids(pos) = node%id
        END IF
        axis = MOD(depth, 2) + 1
        IF (target(axis) < node%coord(axis)) THEN
            next_branch => node%left; opposite_branch => node%right
        ELSE
            next_branch => node%right; opposite_branch => node%left
        END IF
        CALL search_knn(next_branch, target, k, k_ids, k_dists, depth+1)
        IF (ABS(target(axis) - node%coord(axis)) < k_dists(k)) THEN
            CALL search_knn(opposite_branch, target, k, k_ids, k_dists, depth+1)
        END IF
    END SUBROUTINE search_knn

    SUBROUTINE find_knn(root, target, k, neighbors_ids, distances)
        TYPE(kd_node), POINTER, INTENT(IN) :: root
        REAL(DP), DIMENSION(2), INTENT(IN) :: target
        INTEGER, INTENT(IN) :: k
        INTEGER, DIMENSION(k), INTENT(OUT) :: neighbors_ids
        REAL(DP), DIMENSION(k), INTENT(OUT) :: distances
        distances = HUGE(1.0_DP); neighbors_ids = -1
        CALL search_knn(root, target, k, neighbors_ids, distances, 0)
    END SUBROUTINE find_knn
END MODULE kdtree_mod

! ==============================================================================
! MODULO 2: ESTRUCTURAS DISPERSAS CSR
! ==============================================================================
MODULE sparse_mod
    USE types_mod
    IMPLICIT NONE
    TYPE :: csr_matrix
        INTEGER :: n_rows, nnz
        REAL(DP), ALLOCATABLE :: val(:)
        INTEGER, ALLOCATABLE :: col_ind(:), row_ptr(:)
    END TYPE csr_matrix
    CONTAINS
    SUBROUTINE spmv(A, x, y)
        TYPE(csr_matrix), INTENT(IN) :: A
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(OUT) :: y
        INTEGER :: i, j
        y = 0.0_DP
        DO i = 1, A%n_rows
            DO j = A%row_ptr(i), A%row_ptr(i+1) - 1
                y(i) = y(i) + A%val(j) * x(A%col_ind(j))
            END DO
        END DO
    END SUBROUTINE spmv
END MODULE sparse_mod

! ==============================================================================
! MODULO 3: PESOS RBF-FD (NORMALIZACION Y HIPERVISCOSIDAD)
! ==============================================================================
MODULE rbf_weights_mod
    USE types_mod
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE DGETRF(M, N, A, LDA, IPIV, INFO)
            IMPORT :: DP
            INTEGER :: M, N, LDA, INFO
            REAL(DP) :: A(LDA, *)
            INTEGER :: IPIV(*)
        END SUBROUTINE DGETRF
        SUBROUTINE DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            IMPORT :: DP
            CHARACTER(LEN=1) :: TRANS
            INTEGER :: N, NRHS, LDA, LDB, INFO
            REAL(DP) :: A(LDA, *), B(LDB, *)
            INTEGER :: IPIV(*)
        END SUBROUTINE DGETRS
    END INTERFACE

    CONTAINS
    SUBROUTINE compute_weights(stencil_coords, target_coord, w_lap, w_dx, w_dy, w_hyp)
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: stencil_coords
        REAL(DP), DIMENSION(2), INTENT(IN) :: target_coord
        REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: w_lap, w_dx, w_dy, w_hyp
        
        INTEGER :: ns, np, n_tot, i, j, info
        REAL(DP) :: dx, dy, r_norm, R_scale
        REAL(DP), ALLOCATABLE :: A(:,:), B(:,:)
        INTEGER, ALLOCATABLE :: ipiv(:)

        ns = SIZE(stencil_coords, 2); np = 10
        n_tot = ns + np
        ALLOCATE(A(n_tot, n_tot), B(n_tot, 4), ipiv(n_tot))
        ALLOCATE(w_lap(ns), w_dx(ns), w_dy(ns), w_hyp(ns))
        A = 0.0_DP; B = 0.0_DP

        R_scale = 0.0_DP
        DO i = 1, ns
            dx = stencil_coords(1, i) - target_coord(1)
            dy = stencil_coords(2, i) - target_coord(2)
            IF (SQRT(dx**2 + dy**2) > R_scale) R_scale = SQRT(dx**2 + dy**2)
        END DO
        IF (R_scale < 1.0D-14) R_scale = 1.0_DP

        DO i = 1, ns
            DO j = 1, ns
                dx = (stencil_coords(1, i) - stencil_coords(1, j)) / R_scale
                dy = (stencil_coords(2, i) - stencil_coords(2, j)) / R_scale
                A(i, j) = (SQRT(dx**2 + dy**2))**7
            END DO
            dx = (stencil_coords(1, i) - target_coord(1)) / R_scale
            dy = (stencil_coords(2, i) - target_coord(2)) / R_scale
            
            A(i, ns+1) = 1.0_DP;     A(i, ns+2) = dx;         A(i, ns+3) = dy
            A(i, ns+4) = dx**2;      A(i, ns+5) = dx*dy;      A(i, ns+6) = dy**2
            A(i, ns+7) = dx**3;      A(i, ns+8) = dx**2 * dy; A(i, ns+9) = dx * dy**2; A(i, ns+10)= dy**3
            A(ns+1:n_tot, i) = A(i, ns+1:n_tot)
        END DO

        DO i = 1, ns
            dx = (stencil_coords(1, i) - target_coord(1)) / R_scale
            dy = (stencil_coords(2, i) - target_coord(2)) / R_scale
            r_norm = SQRT(dx**2 + dy**2)
            
            B(i, 1) = 49.0_DP * (r_norm**5) / (R_scale**2)             
            B(i, 2) = -7.0_DP * dx * (r_norm**5) / R_scale             
            B(i, 3) = -7.0_DP * dy * (r_norm**5) / R_scale             
            B(i, 4) = 11025.0_DP * r_norm / (R_scale**6)               
        END DO

        B(ns+4, 1) = 2.0_DP / (R_scale**2); B(ns+6, 1) = 2.0_DP / (R_scale**2)
        B(ns+2, 2) = 1.0_DP / R_scale;      B(ns+3, 3) = 1.0_DP / R_scale
        
        CALL DGETRF(n_tot, n_tot, A, n_tot, ipiv, info)
        CALL DGETRS('N', n_tot, 4, A, n_tot, ipiv, B, n_tot, info)

        w_lap = B(1:ns, 1); w_dx = B(1:ns, 2)
        w_dy  = B(1:ns, 3); w_hyp = B(1:ns, 4)
        DEALLOCATE(A, B, ipiv)
    END SUBROUTINE compute_weights
END MODULE rbf_weights_mod

! ==============================================================================
! MODULO 4: SOLVER R-K4 Y BiCGSTAB (DINAMICA DE FLUIDOS AVANZADA)
! ==============================================================================
MODULE solver_ns_mod
    USE types_mod
    USE sparse_mod
    IMPLICIT NONE

    CONTAINS
    
    ! --- Operador libre de Matrices para BiCGSTAB ---
    SUBROUTINE apply_A(x_in, y_out, L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y)
        REAL(DP), DIMENSION(:), INTENT(IN) :: x_in
        REAL(DP), DIMENSION(:), INTENT(OUT) :: y_out
        TYPE(csr_matrix), INTENT(IN) :: L_lap, L_dx, L_dy
        INTEGER, INTENT(IN) :: num_nodes
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type, out_mirror, g_wall_idx
        REAL(DP), DIMENSION(:), INTENT(IN) :: g_norm_x, g_norm_y
        
        REAL(DP), ALLOCATABLE :: lap_x(:), dx_x(:), dy_x(:)
        INTEGER :: i, w_idx
        
        ALLOCATE(lap_x(num_nodes), dx_x(num_nodes), dy_x(num_nodes))
        CALL spmv(L_lap, x_in, lap_x); CALL spmv(L_dx, x_in, dx_x); CALL spmv(L_dy, x_in, dy_x)
        
        DO i = 1, num_nodes
            IF (node_type(i) == 0) THEN 
                y_out(i) = lap_x(i)
            ELSE IF (node_type(i) == 1 .OR. node_type(i) == 3 .OR. node_type(i) == 5) THEN
                y_out(i) = x_in(i) ! Fronteras de Dirichlet fuertes
            ELSE IF (node_type(i) == 2 .AND. out_mirror(i) > 0) THEN
                y_out(i) = x_in(i) - x_in(out_mirror(i)) ! Gradiente nulo (Neumann d/dx)
            ELSE IF (node_type(i) == 4) THEN
                ! Derivada direccional de pared aplicada al nodo fantasma
                w_idx = g_wall_idx(i)
                y_out(i) = g_norm_x(i) * dx_x(w_idx) + g_norm_y(i) * dy_x(w_idx)
            ELSE
                y_out(i) = x_in(i)
            END IF
        END DO
        DEALLOCATE(lap_x, dx_x, dy_x)
    END SUBROUTINE apply_A

    ! --- Inversor Elíptico BiCGSTAB ---
    SUBROUTINE solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omega)
        TYPE(csr_matrix), INTENT(IN) :: L_lap, L_dx, L_dy
        INTEGER, INTENT(IN) :: num_nodes
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type, out_mirror, g_wall_idx
        REAL(DP), DIMENSION(:), INTENT(IN) :: g_norm_x, g_norm_y
        REAL(DP), INTENT(IN) :: V0, Ly
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: coords
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: psi, omega
        
        REAL(DP), ALLOCATABLE :: b(:), r(:), r0_hat(:), p(:), v(:), s(:), t(:)
        REAL(DP) :: rho, rho_prev, alpha, omega_c, beta, err
        INTEGER :: i, iter
        
        ALLOCATE(b(num_nodes), r(num_nodes), r0_hat(num_nodes), p(num_nodes), v(num_nodes), s(num_nodes), t(num_nodes))
        
        ! 1. Vector derecho estático (Condiciones de Frontera + Término Fuente)
        b = 0.0_DP
        DO i = 1, num_nodes
            IF (node_type(i) == 0) b(i) = -omega(i)
            IF (node_type(i) == 1) b(i) = V0 * coords(2,i)
            IF (node_type(i) == 5) b(i) = V0 * Ly
            ! Wall, Outlet, y Ghost mantienen 0.0 en RHS analitico.
        END DO
        
        ! 2. Ciclo principal BiCGSTAB
        CALL apply_A(psi, v, L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y)
        r = b - v
        r0_hat = r
        rho_prev = 1.0_DP; alpha = 1.0_DP; omega_c = 1.0_DP
        v = 0.0_DP; p = 0.0_DP
        
        DO iter = 1, 2000
            rho = DOT_PRODUCT(r0_hat, r)
            IF (ABS(rho) < 1D-14) EXIT
            
            IF (iter == 1) THEN
                p = r
            ELSE
                beta = (rho / rho_prev) * (alpha / omega_c)
                p = r + beta * (p - omega_c * v)
            END IF
            
            CALL apply_A(p, v, L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y)
            alpha = rho / DOT_PRODUCT(r0_hat, v)
            s = r - alpha * v
            
            err = SQRT(DOT_PRODUCT(s, s)) / SQRT(DOT_PRODUCT(b,b) + 1D-14)
            IF (err < 1D-6) THEN
                psi = psi + alpha * p; EXIT
            END IF
            
            CALL apply_A(s, t, L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y)
            omega_c = DOT_PRODUCT(t, s) / DOT_PRODUCT(t, t)
            
            psi = psi + alpha * p + omega_c * s
            r = s - omega_c * t
            
            err = SQRT(DOT_PRODUCT(r, r)) / SQRT(DOT_PRODUCT(b,b) + 1D-14)
            IF (err < 1D-6) EXIT
            rho_prev = rho
        END DO
        DEALLOCATE(b, r, r0_hat, p, v, s, t)
    END SUBROUTINE solve_poisson_bicgstab

    SUBROUTINE update_boundary_vorticity(L_lap, num_nodes, node_type, psi, omega)
        TYPE(csr_matrix), INTENT(IN) :: L_lap
        INTEGER, INTENT(IN) :: num_nodes
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type
        REAL(DP), DIMENSION(:), INTENT(IN) :: psi
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: omega
        REAL(DP), ALLOCATABLE :: lap_psi(:)
        INTEGER :: i
        ALLOCATE(lap_psi(num_nodes))
        CALL spmv(L_lap, psi, lap_psi)
        DO i = 1, num_nodes
            IF (node_type(i) == 3) omega(i) = -lap_psi(i)
        END DO
        DEALLOCATE(lap_psi)
    END SUBROUTINE update_boundary_vorticity

    SUBROUTINE compute_rhs(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omega, rhs_out)
        REAL(DP), INTENT(IN) :: Re, gamma_hyp
        INTEGER, INTENT(IN) :: num_nodes
        TYPE(csr_matrix), INTENT(IN) :: L_lap, L_dx, L_dy, L_hyp
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type
        REAL(DP), DIMENSION(:), INTENT(IN) :: psi, omega
        REAL(DP), DIMENSION(:), INTENT(OUT) :: rhs_out
        
        REAL(DP), ALLOCATABLE :: dpsi_dx(:), dpsi_dy(:), domg_dx(:), domg_dy(:), lap_omg(:), hyp_omg(:)
        INTEGER :: i
        ALLOCATE(dpsi_dx(num_nodes), dpsi_dy(num_nodes), domg_dx(num_nodes), domg_dy(num_nodes))
        ALLOCATE(lap_omg(num_nodes), hyp_omg(num_nodes))
        
        CALL spmv(L_dx, psi, dpsi_dx); CALL spmv(L_dy, psi, dpsi_dy)
        CALL spmv(L_dx, omega, domg_dx); CALL spmv(L_dy, omega, domg_dy)
        CALL spmv(L_lap, omega, lap_omg); CALL spmv(L_hyp, omega, hyp_omg)
        
        rhs_out = 0.0_DP
        DO i = 1, num_nodes
            IF (node_type(i) == 0) THEN 
                rhs_out(i) = (1.0_DP/Re)*lap_omg(i) - (dpsi_dy(i)*domg_dx(i) - dpsi_dx(i)*domg_dy(i)) + gamma_hyp*hyp_omg(i)
            END IF
        END DO
        DEALLOCATE(dpsi_dx, dpsi_dy, domg_dx, domg_dy, lap_omg, hyp_omg)
    END SUBROUTINE compute_rhs

    SUBROUTINE solve_ns(Re, V0, Ly, dt, num_nodes, coords, L_lap, L_dx, L_dy, L_hyp, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, psi, omega)
        REAL(DP), INTENT(IN) :: Re, V0, Ly, dt
        INTEGER, INTENT(IN) :: num_nodes
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: coords
        TYPE(csr_matrix), INTENT(IN) :: L_lap, L_dx, L_dy, L_hyp
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type, out_mirror, g_wall_idx
        REAL(DP), DIMENSION(:), INTENT(IN) :: g_norm_x, g_norm_y
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: psi, omega

        INTEGER :: iter, max_iter, i
        REAL(DP), ALLOCATABLE :: omg_n(:), k1(:), k2(:), k3(:), k4(:), omg_temp(:)
        REAL(DP) :: err_omega, gamma_hyp

        gamma_hyp = 0.01_DP * (0.05_DP**6) 
        max_iter = 10000

        ALLOCATE(omg_n(num_nodes), omg_temp(num_nodes))
        ALLOCATE(k1(num_nodes), k2(num_nodes), k3(num_nodes), k4(num_nodes))

        DO iter = 1, max_iter
            omg_n = omega
            
            ! --- ETAPA RK 1 ---
            CALL solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omega)
            CALL update_boundary_vorticity(L_lap, num_nodes, node_type, psi, omega)
            CALL compute_rhs(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omega, k1)
            
            ! --- ETAPA RK 2 ---
            omg_temp = omg_n + 0.5_DP * dt * k1
            CALL solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omg_temp)
            CALL update_boundary_vorticity(L_lap, num_nodes, node_type, psi, omg_temp)
            CALL compute_rhs(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omg_temp, k2)

            ! --- ETAPA RK 3 ---
            omg_temp = omg_n + 0.5_DP * dt * k2
            CALL solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omg_temp)
            CALL update_boundary_vorticity(L_lap, num_nodes, node_type, psi, omg_temp)
            CALL compute_rhs(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omg_temp, k3)

            ! --- ETAPA RK 4 ---
            omg_temp = omg_n + dt * k3
            CALL solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omg_temp)
            CALL update_boundary_vorticity(L_lap, num_nodes, node_type, psi, omg_temp)
            CALL compute_rhs(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omg_temp, k4)

            DO i = 1, num_nodes
                IF (node_type(i) == 0) THEN
                    omega(i) = omg_n(i) + (dt / 6.0_DP) * (k1(i) + 2.0_DP*k2(i) + 2.0_DP*k3(i) + k4(i))
                ELSE IF (node_type(i) == 2 .AND. out_mirror(i) > 0) THEN
                    omega(i) = omega(out_mirror(i))
                END IF
            END DO

            err_omega = MAXVAL(ABS(omega - omg_n))
            IF (MOD(iter, 10) == 0) THEN
                PRINT '(A,I6, A,ES10.3)', 'RK4 Step: ', iter, ' Err_Om: ', err_omega
                IF (err_omega > 1.0D5 .OR. err_omega /= err_omega) STOP 'DIVERGENCIA ESPECTRAL DETECTADA'
            END IF
            IF (err_omega < 1D-7) EXIT
        END DO
    END SUBROUTINE solve_ns
END MODULE solver_ns_mod

! ==============================================================================
! PROGRAMA PRINCIPAL
! ==============================================================================
PROGRAM main_navier_stokes
    USE types_mod
    USE kdtree_mod
    USE sparse_mod
    USE rbf_weights_mod
    USE solver_ns_mod
    IMPLICIT NONE

    REAL(DP) :: Re = 2.0_DP, V0 = 1.0_DP, h = 0.05_DP, dt_rk = 0.00005_DP
    REAL(DP) :: Lx = 2.0_DP, Ly = 1.0_DP, obs_x1 = 0.5_DP, obs_x2 = 1.0_DP, obs_y = 0.5_DP
    INTEGER  :: stencil_size = 31, max_nodes = 10000, n_points
    INTEGER  :: i, j, nx, ny, idx, u_vort, u_flujo
    REAL(DP) :: x, y
    
    REAL(DP), ALLOCATABLE :: coords(:,:), psi(:), omega(:)
    INTEGER, ALLOCATABLE :: node_type(:), out_mirror(:), id_map(:,:)
    INTEGER, ALLOCATABLE :: g_wall_idx(:)
    REAL(DP), ALLOCATABLE :: g_norm_x(:), g_norm_y(:)
    
    TYPE(kd_node), POINTER :: root => NULL()
    INTEGER, ALLOCATABLE :: ids_build(:), knn_ids(:)
    REAL(DP), ALLOCATABLE :: coords_build(:,:), knn_dists(:), w_l(:), w_x(:), w_y(:), w_h(:)
    TYPE(csr_matrix) :: M_lap, M_dx, M_dy, M_hyp

    ALLOCATE(coords(2, max_nodes), psi(max_nodes), omega(max_nodes))
    ALLOCATE(node_type(max_nodes), out_mirror(max_nodes))
    ALLOCATE(g_wall_idx(max_nodes), g_norm_x(max_nodes), g_norm_y(max_nodes))
    
    g_wall_idx = 1; g_norm_x = 0.0_DP; g_norm_y = 0.0_DP
    n_points = 0
    nx = NINT(Lx/h) + 1; ny = NINT(Ly/h) + 1
    ALLOCATE(id_map(0:nx+1, -1:ny+1)); id_map = -1

    DO i = 0, nx
        DO j = -1, ny
            x = i * h; y = j * h
            IF (x > obs_x1+1e-4 .AND. x < obs_x2-1e-4 .AND. y >= 0.0 .AND. y < obs_y-1e-4) CYCLE
            n_points = n_points + 1
            id_map(i, j) = n_points
            coords(1, n_points) = x; coords(2, n_points) = y
            psi(n_points) = 0.0_DP; omega(n_points) = 0.0_DP
            out_mirror(n_points) = -1
        END DO
    END DO

    DO i = 0, nx
        DO j = -1, ny
            idx = id_map(i, j); IF (idx < 1) CYCLE
            x = coords(1,idx); y = coords(2,idx)

            IF (i == 0 .AND. j >= 0 .AND. j <= ny) THEN
                node_type(idx) = 1; psi(idx) = V0 * y
            ELSE IF (i == nx .AND. j >= 0 .AND. j <= ny) THEN
                node_type(idx) = 2; IF (id_map(i-1, j) > 0) out_mirror(idx) = id_map(i-1, j)
            ELSE IF (j == ny) THEN
                node_type(idx) = 5; psi(idx) = V0 * Ly
            ELSE IF (j == -1 .AND. (x < obs_x1 .OR. x > obs_x2)) THEN
                node_type(idx) = 4; 
                g_wall_idx(idx) = id_map(i, 0)
                g_norm_x(idx) = 0.0_DP; g_norm_y(idx) = -1.0_DP
            ELSE IF (j == 0 .AND. (x <= obs_x1 .OR. x >= obs_x2)) THEN
                node_type(idx) = 3 
            ELSE IF (x >= obs_x1-1e-4 .AND. x <= obs_x2+1e-4 .AND. y >= -h-1e-4 .AND. y <= obs_y+1e-4) THEN
                IF (ABS(y - obs_y) < 1e-4 .OR. ABS(x - obs_x1) < 1e-4 .OR. ABS(x - obs_x2) < 1e-4) THEN
                    IF (j >= 0) THEN
                        node_type(idx) = 3 
                    ELSE
                        node_type(idx) = 4; g_wall_idx(idx) = id_map(i, 0)
                        g_norm_x(idx) = 0.0_DP; g_norm_y(idx) = -1.0_DP
                    END IF
                ELSE
                    node_type(idx) = 4
                    IF (ABS(y - (obs_y - h)) < 1e-4 .AND. id_map(i, j+1) > 0) THEN
                        g_wall_idx(idx) = id_map(i, j+1)
                        g_norm_x(idx) = 0.0_DP; g_norm_y(idx) = -1.0_DP
                    ELSE IF (ABS(x - (obs_x1 + h)) < 1e-4 .AND. id_map(i-1, j) > 0) THEN
                        g_wall_idx(idx) = id_map(i-1, j)
                        g_norm_x(idx) = 1.0_DP; g_norm_y(idx) = 0.0_DP
                    ELSE IF (ABS(x - (obs_x2 - h)) < 1e-4 .AND. id_map(i+1, j) > 0) THEN
                        g_wall_idx(idx) = id_map(i+1, j)
                        g_norm_x(idx) = -1.0_DP; g_norm_y(idx) = 0.0_DP
                    END IF
                END IF
            ELSE
                node_type(idx) = 0 
            END IF
        END DO
    END DO

    ALLOCATE(ids_build(n_points), coords_build(2, n_points))
    DO i = 1, n_points; ids_build(i) = i; coords_build(:, i) = coords(:, i); END DO
    CALL build_kdtree(ids_build, coords_build, n_points, 0, root)

    M_lap%n_rows = n_points; M_lap%nnz = n_points * stencil_size; M_dx = M_lap; M_dy = M_lap; M_hyp = M_lap
    ALLOCATE(M_lap%row_ptr(n_points+1), M_lap%col_ind(M_lap%nnz), M_lap%val(M_lap%nnz))
    ALLOCATE(M_dx%row_ptr(n_points+1),  M_dx%col_ind(M_dx%nnz),   M_dx%val(M_dx%nnz))
    ALLOCATE(M_dy%row_ptr(n_points+1),  M_dy%col_ind(M_dy%nnz),   M_dy%val(M_dy%nnz))
    ALLOCATE(M_hyp%row_ptr(n_points+1), M_hyp%col_ind(M_hyp%nnz), M_hyp%val(M_hyp%nnz))
    
    ALLOCATE(knn_ids(stencil_size), knn_dists(stencil_size))
    M_lap%row_ptr(1) = 1; M_dx%row_ptr(1) = 1; M_dy%row_ptr(1) = 1; M_hyp%row_ptr(1) = 1

    DO i = 1, n_points
        CALL find_knn(root, coords(:, i), stencil_size, knn_ids, knn_dists)
        DO j = 1, stencil_size; coords_build(:, j) = coords(:, knn_ids(j)); END DO
        CALL compute_weights(coords_build(:, 1:stencil_size), coords(:, i), w_l, w_x, w_y, w_h)
        DO j = 1, stencil_size
            idx = (i-1)*stencil_size + j
            M_lap%col_ind(idx) = knn_ids(j); M_lap%val(idx) = w_l(j)
            M_dx%col_ind(idx)  = knn_ids(j); M_dx%val(idx)  = w_x(j)
            M_dy%col_ind(idx)  = knn_ids(j); M_dy%val(idx)  = w_y(j)
            M_hyp%col_ind(idx) = knn_ids(j); M_hyp%val(idx) = w_h(j)
        END DO
        M_lap%row_ptr(i+1) = i*stencil_size + 1; M_dx%row_ptr(i+1) = i*stencil_size + 1
        M_dy%row_ptr(i+1)  = i*stencil_size + 1; M_hyp%row_ptr(i+1) = i*stencil_size + 1
    END DO

    PRINT *, 'Iniciando simulacion transitoria (RK4 + BiCGSTAB Acoplado)...'
    CALL solve_ns(Re, V0, Ly, dt_rk, n_points, coords, M_lap, M_dx, M_dy, M_hyp, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, psi, omega)

    OPEN(NEWUNIT=u_vort, FILE='VORT', STATUS='REPLACE')
    OPEN(NEWUNIT=u_flujo, FILE='FLUJO', STATUS='REPLACE')
    WRITE(u_vort,  '(A)') '# X            Y            Omega'
    WRITE(u_flujo, '(A)') '# X            Y            Psi'
    DO i = 1, n_points
        IF (node_type(i) /= 4) THEN 
            WRITE(u_vort,  '(3(ES15.6, 2X))') coords(1,i), coords(2,i), omega(i)
            WRITE(u_flujo, '(3(ES15.6, 2X))') coords(1,i), coords(2,i), psi(i)
        END IF
    END DO
    CLOSE(u_vort); CLOSE(u_flujo)
    PRINT *, 'Simulación RK4 completada. Dominio estabilizado exitosamente.'
END PROGRAM main_navier_stokes

