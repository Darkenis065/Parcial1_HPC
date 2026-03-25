! ==============================================================================
! SOLUCIONADOR RBF-FD NAVIER-STOKES 2D (FLUJO-VORTICIDAD)
! ARQUITECTURA HPC: Single-Core Optimizado, Fusión CSR, KD-Tree Max-Heap, BiCGSTAB
! ==============================================================================
MODULE types_mod
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
END MODULE types_mod

! ==============================================================================
! MODULO 1: KD-TREE Y TOPOLOGÍA DE VECINDAD
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

    SUBROUTINE partition(ids, coords, left, right, axis, pivot_idx)
        INTEGER, DIMENSION(:), INTENT(INOUT) :: ids
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: coords
        INTEGER, INTENT(IN) :: left, right, axis
        INTEGER, INTENT(OUT) :: pivot_idx
        REAL(DP) :: pivot_val, temp_c(2)
        INTEGER :: i, j, temp_id

        pivot_val = coords(axis, right)
        i = left - 1
        DO j = left, right - 1
            IF (coords(axis, j) <= pivot_val) THEN
                i = i + 1
                temp_c = coords(:, i); coords(:, i) = coords(:, j); coords(:, j) = temp_c
                temp_id = ids(i); ids(i) = ids(j); ids(j) = temp_id
            END IF
        END DO
        i = i + 1
        temp_c = coords(:, i); coords(:, i) = coords(:, right); coords(:, right) = temp_c
        temp_id = ids(i); ids(i) = ids(right); ids(right) = temp_id
        pivot_idx = i
    END SUBROUTINE partition

    SUBROUTINE quickselect(ids, coords, n, k_target, axis)
        INTEGER, DIMENSION(:), INTENT(INOUT) :: ids
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: coords
        INTEGER, INTENT(IN) :: n, k_target, axis
        INTEGER :: left, right, pivot_idx
        left = 1; right = n
        DO WHILE (left < right)
            CALL partition(ids, coords, left, right, axis, pivot_idx)
            IF (pivot_idx == k_target) THEN
                EXIT
            ELSE IF (k_target < pivot_idx) THEN
                right = pivot_idx - 1
            ELSE
                left = pivot_idx + 1
            END IF
        END DO
    END SUBROUTINE quickselect

    RECURSIVE SUBROUTINE build_kdtree(ids, coords, n, depth, root)
        INTEGER, DIMENSION(:), INTENT(INOUT) :: ids
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: coords
        INTEGER, INTENT(IN) :: n, depth
        TYPE(kd_node), POINTER, INTENT(OUT) :: root
        INTEGER :: axis, median_idx

        IF (n <= 0) THEN
            root => NULL(); RETURN
        END IF

        ALLOCATE(root)
        axis = MOD(depth, 2) + 1
        median_idx = n / 2 + 1
        
        CALL quickselect(ids, coords, n, median_idx, axis)

        root%id = ids(median_idx)
        root%coord = coords(:, median_idx)
        CALL build_kdtree(ids(1:median_idx-1), coords(:, 1:median_idx-1), median_idx-1, depth+1, root%left)
        CALL build_kdtree(ids(median_idx+1:n), coords(:, median_idx+1:n), n-median_idx, depth+1, root%right)
    END SUBROUTINE build_kdtree

    RECURSIVE SUBROUTINE heapify_down(k_dists, k_ids, k, idx)
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: k_dists
        INTEGER, DIMENSION(:), INTENT(INOUT) :: k_ids
        INTEGER, INTENT(IN) :: k, idx
        INTEGER :: largest, l, r, temp_id
        REAL(DP) :: temp_dist
        
        largest = idx; l = 2 * idx; r = 2 * idx + 1
        IF (l <= k) THEN; IF (k_dists(l) > k_dists(largest)) largest = l; END IF
        IF (r <= k) THEN; IF (k_dists(r) > k_dists(largest)) largest = r; END IF
        
        IF (largest /= idx) THEN
            temp_dist = k_dists(idx); k_dists(idx) = k_dists(largest); k_dists(largest) = temp_dist
            temp_id = k_ids(idx); k_ids(idx) = k_ids(largest); k_ids(largest) = temp_id
            CALL heapify_down(k_dists, k_ids, k, largest)
        END IF
    END SUBROUTINE heapify_down

    RECURSIVE SUBROUTINE search_knn(node, target, k, k_ids, k_dists, depth)
        TYPE(kd_node), POINTER, INTENT(IN) :: node
        REAL(DP), DIMENSION(2), INTENT(IN) :: target
        INTEGER, INTENT(IN) :: k, depth
        INTEGER, DIMENSION(k), INTENT(INOUT) :: k_ids
        REAL(DP), DIMENSION(k), INTENT(INOUT) :: k_dists
        REAL(DP) :: dist_sq, dx, dy
        INTEGER :: axis
        TYPE(kd_node), POINTER :: next_branch, opposite_branch

        IF (.NOT. ASSOCIATED(node)) RETURN
        dx = node%coord(1) - target(1); dy = node%coord(2) - target(2)
        dist_sq = dx*dx + dy*dy

        IF (dist_sq < k_dists(1)) THEN
            k_dists(1) = dist_sq; k_ids(1) = node%id
            CALL heapify_down(k_dists, k_ids, k, 1)
        END IF

        axis = MOD(depth, 2) + 1
        dx = target(axis) - node%coord(axis)
        IF (dx < 0.0_DP) THEN
            next_branch => node%left; opposite_branch => node%right
        ELSE
            next_branch => node%right; opposite_branch => node%left
        END IF
        
        CALL search_knn(next_branch, target, k, k_ids, k_dists, depth+1)
        IF ((dx*dx) < k_dists(1)) THEN
            CALL search_knn(opposite_branch, target, k, k_ids, k_dists, depth+1)
        END IF
    END SUBROUTINE search_knn

    SUBROUTINE find_knn(root, target, k, neighbors_ids, distances)
        TYPE(kd_node), POINTER, INTENT(IN) :: root
        REAL(DP), DIMENSION(2), INTENT(IN) :: target
        INTEGER, INTENT(IN) :: k
        INTEGER, DIMENSION(k), INTENT(OUT) :: neighbors_ids
        REAL(DP), DIMENSION(k), INTENT(OUT) :: distances
        INTEGER :: i
        distances = HUGE(1.0_DP); neighbors_ids = -1
        CALL search_knn(root, target, k, neighbors_ids, distances, 0)
        DO i = 1, k
            distances(i) = SQRT(distances(i))
        END DO
    END SUBROUTINE find_knn
END MODULE kdtree_mod

! ==============================================================================
! MODULO 2: ARQUITECTURA CSR
! ==============================================================================
MODULE sparse_mod
    USE types_mod
    IMPLICIT NONE
    TYPE :: csr_matrix
        INTEGER :: n_rows, nnz
        REAL(DP), ALLOCATABLE :: val(:)
        INTEGER, ALLOCATABLE :: col_ind(:), row_ptr(:)
    END TYPE csr_matrix
END MODULE sparse_mod

! ==============================================================================
! MODULO 3: PESOS RBF-FD
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
    SUBROUTINE compute_weights(stencil_coords, target_coord, w_lap, w_dx, w_dy, w_hyp, A, B, ipiv)
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: stencil_coords
        REAL(DP), DIMENSION(2), INTENT(IN) :: target_coord
        REAL(DP), DIMENSION(:), INTENT(OUT) :: w_lap, w_dx, w_dy, w_hyp
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: A, B
        INTEGER, DIMENSION(:), INTENT(INOUT) :: ipiv
        
        INTEGER :: ns, np, n_tot, i, j, info
        REAL(DP) :: dx, dy, dist_sq, R_scale_sq, R_scale, inv_R, inv_R2, inv_R6
        REAL(DP) :: r_norm_sq, r_norm_5, r_norm_7

        ns = SIZE(stencil_coords, 2); np = 10; n_tot = ns + np
        A = 0.0_DP; B = 0.0_DP

        R_scale_sq = 0.0_DP
        DO i = 1, ns
            dx = stencil_coords(1, i) - target_coord(1)
            dy = stencil_coords(2, i) - target_coord(2)
            dist_sq = dx*dx + dy*dy
            IF (dist_sq > R_scale_sq) R_scale_sq = dist_sq
        END DO
        IF (R_scale_sq < 1.0D-28) THEN
            R_scale = 1.0_DP; inv_R = 1.0_DP
        ELSE
            R_scale = SQRT(R_scale_sq); inv_R = 1.0_DP / R_scale
        END IF
        inv_R2 = inv_R * inv_R; inv_R6 = inv_R2 * inv_R2 * inv_R2

        DO i = 1, ns
            DO j = 1, ns
                dx = (stencil_coords(1, i) - stencil_coords(1, j)) * inv_R
                dy = (stencil_coords(2, i) - stencil_coords(2, j)) * inv_R
                r_norm_sq = dx*dx + dy*dy
                r_norm_7 = r_norm_sq * r_norm_sq * r_norm_sq * SQRT(r_norm_sq)
                A(i, j) = r_norm_7
            END DO
            dx = (stencil_coords(1, i) - target_coord(1)) * inv_R
            dy = (stencil_coords(2, i) - target_coord(2)) * inv_R
            
            A(i, ns+1) = 1.0_DP;     A(i, ns+2) = dx;         A(i, ns+3) = dy
            A(i, ns+4) = dx*dx;      A(i, ns+5) = dx*dy;      A(i, ns+6) = dy*dy
            A(i, ns+7) = dx*dx*dx;   A(i, ns+8) = dx*dx*dy;   A(i, ns+9) = dx*dy*dy; A(i, ns+10)= dy*dy*dy
            A(ns+1:n_tot, i) = A(i, ns+1:n_tot)
        END DO

        DO i = 1, ns
            dx = (stencil_coords(1, i) - target_coord(1)) * inv_R
            dy = (stencil_coords(2, i) - target_coord(2)) * inv_R
            r_norm_sq = dx*dx + dy*dy
            r_norm_5 = r_norm_sq * r_norm_sq * SQRT(r_norm_sq)
            
            B(i, 1) = 49.0_DP * r_norm_5 * inv_R2             
            B(i, 2) = -7.0_DP * dx * r_norm_5 * inv_R             
            B(i, 3) = -7.0_DP * dy * r_norm_5 * inv_R             
            B(i, 4) = 11025.0_DP * SQRT(r_norm_sq) * inv_R6               
        END DO

        B(ns+4, 1) = 2.0_DP * inv_R2; B(ns+6, 1) = 2.0_DP * inv_R2
        B(ns+2, 2) = 1.0_DP * inv_R;  B(ns+3, 3) = 1.0_DP * inv_R
        
        CALL DGETRF(n_tot, n_tot, A, n_tot, ipiv, info)
        CALL DGETRS('N', n_tot, 4, A, n_tot, ipiv, B, n_tot, info)

        w_lap = B(1:ns, 1); w_dx = B(1:ns, 2)
        w_dy  = B(1:ns, 3); w_hyp = B(1:ns, 4)
    END SUBROUTINE compute_weights
END MODULE rbf_weights_mod

! ==============================================================================
! MODULO 4: SOLVER R-K4 Y BiCGSTAB
! ==============================================================================
MODULE solver_ns_mod
    USE types_mod
    USE sparse_mod
    IMPLICIT NONE

    CONTAINS
    
    SUBROUTINE apply_A_fused(x_in, y_out, lap_x, dx_x, dy_x, L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y)
        REAL(DP), DIMENSION(:), INTENT(IN) :: x_in
        REAL(DP), DIMENSION(:), INTENT(OUT) :: y_out, lap_x, dx_x, dy_x
        TYPE(csr_matrix), INTENT(IN) :: L_lap, L_dx, L_dy
        INTEGER, INTENT(IN) :: num_nodes
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type, out_mirror, g_wall_idx
        REAL(DP), DIMENSION(:), INTENT(IN) :: g_norm_x, g_norm_y
        INTEGER :: i, j, c, w_idx
        REAL(DP) :: x_val, s_lap, s_dx, s_dy
        
        DO i = 1, num_nodes
            s_lap = 0.0_DP; s_dx = 0.0_DP; s_dy = 0.0_DP
            DO j = L_lap%row_ptr(i), L_lap%row_ptr(i+1)-1
                c = L_lap%col_ind(j); x_val = x_in(c)
                s_lap = s_lap + L_lap%val(j) * x_val
                s_dx  = s_dx  + L_dx%val(j)  * x_val
                s_dy  = s_dy  + L_dy%val(j)  * x_val
            END DO
            lap_x(i) = s_lap; dx_x(i) = s_dx; dy_x(i) = s_dy
        END DO
        
        DO i = 1, num_nodes
            IF (node_type(i) == 0) THEN 
                y_out(i) = lap_x(i)
            ELSE IF (node_type(i) == 1 .OR. node_type(i) == 3 .OR. node_type(i) == 5) THEN
                y_out(i) = x_in(i)
            ELSE IF (node_type(i) == 2 .AND. out_mirror(i) > 0) THEN
                y_out(i) = x_in(i) - x_in(out_mirror(i))
            ELSE IF (node_type(i) == 4) THEN
                w_idx = g_wall_idx(i)
                y_out(i) = g_norm_x(i) * dx_x(w_idx) + g_norm_y(i) * dy_x(w_idx)
            ELSE
                y_out(i) = x_in(i)
            END IF
        END DO
    END SUBROUTINE apply_A_fused

    SUBROUTINE solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omega, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x)
        TYPE(csr_matrix), INTENT(IN) :: L_lap, L_dx, L_dy
        INTEGER, INTENT(IN) :: num_nodes
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type, out_mirror, g_wall_idx
        REAL(DP), DIMENSION(:), INTENT(IN) :: g_norm_x, g_norm_y, omega
        REAL(DP), INTENT(IN) :: V0, Ly
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: coords
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: psi, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x
        REAL(DP) :: rho, rho_prev, alpha, omega_c, beta, err, dot_r0_r, dot_r0_v, dot_s_s, dot_t_s, dot_t_t, dot_r_r, dot_b_b
        INTEGER :: i, iter
        
        DO i = 1, num_nodes
            IF (node_type(i) == 0) THEN; b(i) = -omega(i)
            ELSE IF (node_type(i) == 1) THEN; b(i) = V0 * coords(2,i)
            ELSE IF (node_type(i) == 5) THEN; b(i) = V0 * Ly
            ELSE; b(i) = 0.0_DP; END IF
        END DO
        
        CALL apply_A_fused(psi, v, lap_x, dx_x, dy_x, L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y)
        
        dot_b_b = 0.0_DP
        DO i = 1, num_nodes
            r(i) = b(i) - v(i); r0_hat(i) = r(i)
            v(i) = 0.0_DP; p(i) = 0.0_DP
            dot_b_b = dot_b_b + b(i)*b(i)
        END DO
        rho_prev = 1.0_DP; alpha = 1.0_DP; omega_c = 1.0_DP
        
        DO iter = 1, 2000
            dot_r0_r = 0.0_DP
            DO i = 1, num_nodes; dot_r0_r = dot_r0_r + r0_hat(i) * r(i); END DO
            rho = dot_r0_r
            IF (ABS(rho) < 1D-14) EXIT
            
            IF (iter == 1) THEN
                DO i = 1, num_nodes; p(i) = r(i); END DO
            ELSE
                beta = (rho / rho_prev) * (alpha / omega_c)
                DO i = 1, num_nodes
                    p(i) = r(i) + beta * (p(i) - omega_c * v(i))
                END DO
            END IF
            
            CALL apply_A_fused(p, v, lap_x, dx_x, dy_x, L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y)
            
            dot_r0_v = 0.0_DP
            DO i = 1, num_nodes; dot_r0_v = dot_r0_v + r0_hat(i) * v(i); END DO
            IF (ABS(dot_r0_v) < 1D-14) dot_r0_v = SIGN(1.0D-14, dot_r0_v) ! Blindaje singular
            alpha = rho / dot_r0_v
            
            dot_s_s = 0.0_DP
            DO i = 1, num_nodes
                s(i) = r(i) - alpha * v(i)
                dot_s_s = dot_s_s + s(i)*s(i)
            END DO
            
            err = SQRT(dot_s_s) / (SQRT(dot_b_b) + 1D-14)
            IF (err < 1D-6) THEN
                DO i = 1, num_nodes; psi(i) = psi(i) + alpha * p(i); END DO
                EXIT
            END IF
            
            CALL apply_A_fused(s, t, lap_x, dx_x, dy_x, L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y)
            
            dot_t_s = 0.0_DP; dot_t_t = 0.0_DP
            DO i = 1, num_nodes
                dot_t_s = dot_t_s + t(i)*s(i); dot_t_t = dot_t_t + t(i)*t(i)
            END DO
            IF (ABS(dot_t_t) < 1D-14) dot_t_t = SIGN(1.0D-14, dot_t_t) ! Blindaje singular
            omega_c = dot_t_s / dot_t_t
            
            dot_r_r = 0.0_DP
            DO i = 1, num_nodes
                psi(i) = psi(i) + alpha * p(i) + omega_c * s(i)
                r(i) = s(i) - omega_c * t(i)
                dot_r_r = dot_r_r + r(i)*r(i)
            END DO
            
            err = SQRT(dot_r_r) / (SQRT(dot_b_b) + 1D-14)
            IF (err < 1D-6) EXIT
            rho_prev = rho
        END DO
    END SUBROUTINE solve_poisson_bicgstab

    ! --- CORRECCION 2: Extrapolación de Vorticidad para estabilizar pared ---
    SUBROUTINE update_boundary_vorticity(L_lap, num_nodes, node_type, g_wall_idx, psi, omega)
        TYPE(csr_matrix), INTENT(IN) :: L_lap
        INTEGER, INTENT(IN) :: num_nodes
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type, g_wall_idx
        REAL(DP), DIMENSION(:), INTENT(IN) :: psi
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: omega
        INTEGER :: i, j
        REAL(DP) :: lap_val
        
        DO i = 1, num_nodes
            IF (node_type(i) == 3) THEN
                lap_val = 0.0_DP
                DO j = L_lap%row_ptr(i), L_lap%row_ptr(i+1) - 1
                    lap_val = lap_val + L_lap%val(j) * psi(L_lap%col_ind(j))
                END DO
                omega(i) = -lap_val
            END IF
        END DO
        ! Prevención del acantilado limitrofe propagando la pared hacia el fantasma
        DO i = 1, num_nodes
            IF (node_type(i) == 4) omega(i) = omega(g_wall_idx(i))
        END DO
    END SUBROUTINE update_boundary_vorticity

    ! --- CORRECCION 1: Eliminación del Bucle Infinito Matemático ---
    SUBROUTINE compute_rhs_fused(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omega, rhs_out)
        REAL(DP), INTENT(IN) :: Re, gamma_hyp
        INTEGER, INTENT(IN) :: num_nodes
        TYPE(csr_matrix), INTENT(IN) :: L_lap, L_dx, L_dy, L_hyp
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type
        REAL(DP), DIMENSION(:), INTENT(IN) :: psi, omega
        REAL(DP), DIMENSION(:), INTENT(OUT) :: rhs_out
        INTEGER :: i, j, c
        REAL(DP) :: sl_omg, sdx_psi, sdy_psi, sdx_omg, sdy_omg, shyp_omg, p_val, o_val
        
        DO i = 1, num_nodes
            sl_omg = 0.0_DP; sdx_psi = 0.0_DP; sdy_psi = 0.0_DP; sdx_omg = 0.0_DP; sdy_omg = 0.0_DP; shyp_omg = 0.0_DP
            DO j = L_lap%row_ptr(i), L_lap%row_ptr(i+1)-1
                c = L_lap%col_ind(j); p_val = psi(c); o_val = omega(c)
                sl_omg = sl_omg + L_lap%val(j) * o_val
                sdx_psi = sdx_psi + L_dx%val(j) * p_val; sdy_psi = sdy_psi + L_dy%val(j) * p_val
                sdx_omg = sdx_omg + L_dx%val(j) * o_val; sdy_omg = sdy_omg + L_dy%val(j) * o_val
                shyp_omg = shyp_omg + L_hyp%val(j) * o_val
            END DO
            IF (node_type(i) == 0) THEN 
                rhs_out(i) = (1.0_DP/Re)*sl_omg - (sdy_psi*sdx_omg - sdx_psi*sdy_omg) + gamma_hyp*shyp_omg
            ELSE
                ! La frontera NO se retroalimenta en RHS, ya que se actualiza en update_boundary_vorticity
                rhs_out(i) = 0.0_DP 
            END IF
        END DO
    END SUBROUTINE compute_rhs_fused

    SUBROUTINE solve_ns(Re, V0, Ly, dt, num_nodes, coords, L_lap, L_dx, L_dy, L_hyp, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, psi, omega, &
                        omg_n, k1, k2, k3, k4, omg_temp, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x)
        REAL(DP), INTENT(IN) :: Re, V0, Ly, dt
        INTEGER, INTENT(IN) :: num_nodes
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: coords
        TYPE(csr_matrix), INTENT(IN) :: L_lap, L_dx, L_dy, L_hyp
        INTEGER, DIMENSION(:), INTENT(IN) :: node_type, out_mirror, g_wall_idx
        REAL(DP), DIMENSION(:), INTENT(IN) :: g_norm_x, g_norm_y
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: psi, omega, omg_n, k1, k2, k3, k4, omg_temp, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x

        INTEGER :: iter, max_iter, i
        REAL(DP) :: err_omega, gamma_hyp

        gamma_hyp = 0.01_DP * (0.05_DP**6); max_iter = 250

        DO iter = 1, max_iter
            DO i = 1, num_nodes; omg_n(i) = omega(i); END DO
            
            ! --- ETAPA RK 1 ---
            CALL solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omega, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x)
            CALL update_boundary_vorticity(L_lap, num_nodes, node_type, g_wall_idx, psi, omega)
            CALL compute_rhs_fused(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omega, k1)
            
            ! --- ETAPA RK 2 ---
            DO i = 1, num_nodes; omg_temp(i) = omg_n(i) + 0.5_DP * dt * k1(i); END DO
            CALL solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omg_temp, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x)
            CALL update_boundary_vorticity(L_lap, num_nodes, node_type, g_wall_idx, psi, omg_temp)
            CALL compute_rhs_fused(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omg_temp, k2)

            ! --- ETAPA RK 3 ---
            DO i = 1, num_nodes; omg_temp(i) = omg_n(i) + 0.5_DP * dt * k2(i); END DO
            CALL solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omg_temp, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x)
            CALL update_boundary_vorticity(L_lap, num_nodes, node_type, g_wall_idx, psi, omg_temp)
            CALL compute_rhs_fused(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omg_temp, k3)

            ! --- ETAPA RK 4 ---
            DO i = 1, num_nodes; omg_temp(i) = omg_n(i) + dt * k3(i); END DO
            CALL solve_poisson_bicgstab(L_lap, L_dx, L_dy, num_nodes, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, V0, Ly, coords, psi, omg_temp, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x)
            CALL update_boundary_vorticity(L_lap, num_nodes, node_type, g_wall_idx, psi, omg_temp)
            CALL compute_rhs_fused(Re, gamma_hyp, num_nodes, L_lap, L_dx, L_dy, L_hyp, node_type, psi, omg_temp, k4)

            err_omega = 0.0_DP
            DO i = 1, num_nodes
                IF (node_type(i) == 0) THEN
                    omega(i) = omg_n(i) + (dt / 6.0_DP) * (k1(i) + 2.0_DP*k2(i) + 2.0_DP*k3(i) + k4(i))
                ELSE IF (node_type(i) == 2 .AND. out_mirror(i) > 0) THEN
                    omega(i) = omega(out_mirror(i))
                END IF
                IF (ABS(omega(i) - omg_n(i)) > err_omega) err_omega = ABS(omega(i) - omg_n(i))
            END DO

            IF (MOD(iter, 50) == 0) THEN
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
    
    REAL(DP), ALLOCATABLE :: coords(:,:), psi(:), omega(:), omg_n(:), k1(:), k2(:), k3(:), k4(:), omg_temp(:)
    REAL(DP), ALLOCATABLE :: b(:), r(:), r0_hat(:), p(:), v(:), s(:), t(:), lap_x(:), dx_x(:), dy_x(:)
    INTEGER, ALLOCATABLE :: node_type(:), out_mirror(:), id_map(:,:), g_wall_idx(:)
    REAL(DP), ALLOCATABLE :: g_norm_x(:), g_norm_y(:)
    
    TYPE(kd_node), POINTER :: root => NULL()
    INTEGER, ALLOCATABLE :: ids_build(:), knn_ids(:), ipiv(:)
    REAL(DP), ALLOCATABLE :: coords_build(:,:), knn_dists(:), w_l(:), w_x(:), w_y(:), w_h(:)
    REAL(DP), ALLOCATABLE :: A_mat(:,:), B_mat(:,:)
    TYPE(csr_matrix) :: M_lap, M_dx, M_dy, M_hyp

    ALLOCATE(coords(2, max_nodes), psi(max_nodes), omega(max_nodes), omg_n(max_nodes), k1(max_nodes), k2(max_nodes), k3(max_nodes), k4(max_nodes), omg_temp(max_nodes))
    ALLOCATE(b(max_nodes), r(max_nodes), r0_hat(max_nodes), p(max_nodes), v(max_nodes), s(max_nodes), t(max_nodes), lap_x(max_nodes), dx_x(max_nodes), dy_x(max_nodes))
    ALLOCATE(node_type(max_nodes), out_mirror(max_nodes), g_wall_idx(max_nodes), g_norm_x(max_nodes), g_norm_y(max_nodes))
    ALLOCATE(ipiv(stencil_size+10), A_mat(stencil_size+10, stencil_size+10), B_mat(stencil_size+10, 4))
    ALLOCATE(knn_ids(stencil_size), knn_dists(stencil_size), w_l(stencil_size), w_x(stencil_size), w_y(stencil_size), w_h(stencil_size))
    
    g_wall_idx = 1; g_norm_x = 0.0_DP; g_norm_y = 0.0_DP
    n_points = 0; nx = NINT(Lx/h) + 1; ny = NINT(Ly/h) + 1
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
                node_type(idx) = 4; g_wall_idx(idx) = id_map(i, 0); g_norm_x(idx) = 0.0_DP; g_norm_y(idx) = -1.0_DP
            ELSE IF (j == 0 .AND. (x <= obs_x1 .OR. x >= obs_x2)) THEN
                node_type(idx) = 3 
            ELSE IF (x >= obs_x1-1e-4 .AND. x <= obs_x2+1e-4 .AND. y >= -h-1e-4 .AND. y <= obs_y+1e-4) THEN
                IF (ABS(y - obs_y) < 1e-4 .OR. ABS(x - obs_x1) < 1e-4 .OR. ABS(x - obs_x2) < 1e-4) THEN
                    IF (j >= 0) THEN; node_type(idx) = 3 
                    ELSE; node_type(idx) = 4; g_wall_idx(idx) = id_map(i, 0); g_norm_x(idx) = 0.0_DP; g_norm_y(idx) = -1.0_DP; END IF
                ELSE
                    node_type(idx) = 4
                    IF (ABS(y - (obs_y - h)) < 1e-4 .AND. id_map(i, j+1) > 0) THEN
                        g_wall_idx(idx) = id_map(i, j+1); g_norm_x(idx) = 0.0_DP; g_norm_y(idx) = -1.0_DP
                    ELSE IF (ABS(x - (obs_x1 + h)) < 1e-4 .AND. id_map(i-1, j) > 0) THEN
                        g_wall_idx(idx) = id_map(i-1, j); g_norm_x(idx) = 1.0_DP; g_norm_y(idx) = 0.0_DP
                    ELSE IF (ABS(x - (obs_x2 - h)) < 1e-4 .AND. id_map(i+1, j) > 0) THEN
                        g_wall_idx(idx) = id_map(i+1, j); g_norm_x(idx) = -1.0_DP; g_norm_y(idx) = 0.0_DP
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
    
    M_lap%row_ptr(1) = 1; M_dx%row_ptr(1) = 1; M_dy%row_ptr(1) = 1; M_hyp%row_ptr(1) = 1

    DO i = 1, n_points
        CALL find_knn(root, coords(:, i), stencil_size, knn_ids, knn_dists)
        DO j = 1, stencil_size; coords_build(:, j) = coords(:, knn_ids(j)); END DO
        CALL compute_weights(coords_build(:, 1:stencil_size), coords(:, i), w_l, w_x, w_y, w_h, A_mat, B_mat, ipiv)
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
    CALL solve_ns(Re, V0, Ly, dt_rk, n_points, coords, M_lap, M_dx, M_dy, M_hyp, node_type, out_mirror, g_wall_idx, g_norm_x, g_norm_y, psi, omega, &
                  omg_n, k1, k2, k3, k4, omg_temp, b, r, r0_hat, p, v, s, t, lap_x, dx_x, dy_x)

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
    PRINT *, 'Simulación RK4 Completada HPC.'
END PROGRAM main_navier_stokes
