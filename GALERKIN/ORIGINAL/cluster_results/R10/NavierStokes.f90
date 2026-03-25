! ==============================================================================
! SOLVER NAVIER-STOKES 2D - GALERKIN P1 (Ψ-ξ)
! ==============================================================================
module GlobalData
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    ! ==========================================================================
    ! >>> PANEL DE CONTROL PRINCIPAL <<<
    ! Modifica estos valores para ajustar la simulación a tu gusto
    ! ==========================================================================

    ! 1. FÍSICA DEL PROBLEMA
    real(dp) :: V0 = 1.0_dp          ! Velocidad del flujo a la entrada
    real(dp) :: nu = 0.05_dp          ! Viscosidad cinemática (Re = V0*D/nu. Ej: D=0.2, Re=2 -> nu=0.1)

    ! 2. ESTABILIDAD NUMÉRICA (Crucial al cambiar el tamaño de la malla)
    !    - Mallas gruesas  (ej. 50x20):   usa 0.1_dp a 0.2_dp
    !    - Mallas medias   (ej. 150x30):  usa 0.02_dp a 0.05_dp
    !    - Mallas ultra finas (ej. 400x80): usa 0.001_dp a 0.005_dp
    real(dp) :: omega_thom = 0.002_dp ! Factor de sub-relajación para la vorticidad en la pared

    ! 3. LÍMITES DE ITERACIÓN Y CONVERGENCIA
    !    - Mallas finas requieren MÁS iteraciones en ambos parámetros
    integer :: max_picard = 25000    ! Máximo de iteraciones globales (Bucle principal de la física)
    integer :: max_gs     = 15000    ! Máximo de iteraciones del solver interno (Ecuación de Poisson)
    real(dp) :: tol       = 1.0d-6   ! Criterio de parada (Tolerancia de error)

    ! ==========================================================================
    ! >>> FIN DEL PANEL DE CONTROL <<<
    ! ==========================================================================

    ! Variables del sistema (No tocar)
    integer :: n_nodes, n_elements
    real(dp), allocatable :: x(:), y(:)
    integer, allocatable :: conn(:,:)
    integer, allocatable :: node_type(:)

    integer, parameter :: max_adj = 50 
    real(dp), allocatable :: A_val(:,:)
    integer, allocatable :: A_col(:,:)
    integer, allocatable :: num_adj(:)
    real(dp), allocatable :: b_vec(:)

    real(dp), allocatable :: Psi(:), Xi(:), Psi_old(:)
end module GlobalData

module MeshOps
    ! ... (Mantén tu subrutina read_mesh() exactamente igual) ...
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
        allocate(conn(n_elements, 3))
        do i = 1, n_elements
            read(20, *) id_tmp, conn(i,1), conn(i,2), conn(i,3)
        end do
        close(20)
    end subroutine read_mesh
end module MeshOps

module SolverOps
    use GlobalData
    implicit none
contains

    ! ... (Mantén init_sparse, add_to_sparse, assemble_vorticity_transport y assemble_poisson iguales) ...
    subroutine init_sparse()
        A_val = 0.0_dp
        A_col = 0
        num_adj = 0
        b_vec = 0.0_dp
    end subroutine init_sparse

    subroutine add_to_sparse(row, col, val)
        integer :: row, col
        real(dp) :: val
        integer :: k
        do k = 1, num_adj(row)
            if (A_col(k, row) == col) then
                A_val(k, row) = A_val(k, row) + val
                return
            end if
        end do
        num_adj(row) = num_adj(row) + 1
        A_col(num_adj(row), row) = col
        A_val(num_adj(row), row) = val
    end subroutine add_to_sparse

    ! --- FÓRMULA DE THOM (Conectada al Panel de Control) ---
    subroutine update_thom_bc()
        integer :: i, j, closest_int
        real(dp) :: min_dist, d, h_dist, Xi_new

        do i = 1, n_nodes
            if (node_type(i) == 1 .or. node_type(i) >= 5) then
                min_dist = 1.0d30
                closest_int = -1
                do j = 1, n_nodes
                    if (node_type(j) == 0) then
                        d = (x(i)-x(j))**2 + (y(i)-y(j))**2
                        if (d < min_dist) then
                            min_dist = d
                            closest_int = j
                        end if
                    end if
                end do
                
                if (closest_int /= -1) then
                    h_dist = max(sqrt(min_dist), 1.0d-8) 
                    Xi_new = -2.0_dp * (Psi(closest_int) - Psi(i)) / (h_dist**2)
                    
                    ! AQUI USAMOS LA VARIABLE DEL PANEL DE CONTROL
                    Xi(i) = (1.0_dp - omega_thom) * Xi(i) + omega_thom * Xi_new
                end if
            else if (node_type(i) == 2 .or. node_type(i) == 3) then
                Xi(i) = 0.0_dp
            end if
        end do
    end subroutine update_thom_bc

    ! ... (Ensamblajes omitidos por brevedad, usa los que ya tienes que funcionan perfecto) ...
    subroutine assemble_vorticity_transport()
        integer :: e, tmp, i, j
        real(dp) :: x1, y1, x2, y2, x3, y3, detJ, Area
        real(dp) :: b(3), c(3), u, v, C_e, K_e, val
        integer :: n(3)

        call init_sparse()

        do e = 1, n_elements
            n(1)=conn(e,1); n(2)=conn(e,2); n(3)=conn(e,3)
            x1=x(n(1)); y1=y(n(1)); x2=x(n(2)); y2=y(n(2)); x3=x(n(3)); y3=y(n(3))

            detJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            if (detJ < 0.0_dp) then
                tmp=n(2); n(2)=n(3); n(3)=tmp
                x2=x(n(2)); y2=y(n(2)); x3=x(n(3)); y3=y(n(3))
                detJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            end if
            
            Area = 0.5_dp * detJ
            if (Area < 1.0d-12) Area = 1.0d-12 

            b(1)=y2-y3; b(2)=y3-y1; b(3)=y1-y2
            c(1)=x3-x2; c(2)=x1-x3; c(3)=x2-x1

            u =  (Psi(n(1))*c(1) + Psi(n(2))*c(2) + Psi(n(3))*c(3)) / (2.0_dp*Area)
            v = -(Psi(n(1))*b(1) + Psi(n(2))*b(2) + Psi(n(3))*b(3)) / (2.0_dp*Area)

            do i = 1, 3
                do j = 1, 3
                    C_e = (1.0_dp/6.0_dp) * (u*b(j) + v*c(j))
                    K_e = (nu / (4.0_dp*Area)) * (b(i)*b(j) + c(i)*c(j))
                    val = C_e + K_e
                    call add_to_sparse(n(i), n(j), val)
                end do
            end do
        end do
    end subroutine assemble_vorticity_transport

    subroutine assemble_poisson()
        integer :: e, tmp, i, j
        real(dp) :: x1, y1, x2, y2, x3, y3, detJ, Area
        real(dp) :: b(3), c(3), K_e, f_e(3)
        integer :: n(3)

        call init_sparse()

        do e = 1, n_elements
            n(1)=conn(e,1); n(2)=conn(e,2); n(3)=conn(e,3)
            x1=x(n(1)); y1=y(n(1)); x2=x(n(2)); y2=y(n(2)); x3=x(n(3)); y3=y(n(3))

            detJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            if (detJ < 0.0_dp) then
                tmp=n(2); n(2)=n(3); n(3)=tmp
                x2=x(n(2)); y2=y(n(2)); x3=x(n(3)); y3=y(n(3))
                detJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            end if
            
            Area = 0.5_dp * detJ
            if (Area < 1.0d-12) Area = 1.0d-12 

            b(1)=y2-y3; b(2)=y3-y1; b(3)=y1-y2
            c(1)=x3-x2; c(2)=x1-x3; c(3)=x2-x1

            f_e(1) = (Area/12.0_dp) * (2.0_dp*Xi(n(1)) + Xi(n(2)) + Xi(n(3)))
            f_e(2) = (Area/12.0_dp) * (Xi(n(1)) + 2.0_dp*Xi(n(2)) + Xi(n(3)))
            f_e(3) = (Area/12.0_dp) * (Xi(n(1)) + Xi(n(2)) + 2.0_dp*Xi(n(3)))

            do i = 1, 3
                b_vec(n(i)) = b_vec(n(i)) + f_e(i)
                do j = 1, 3
                    K_e = (1.0_dp / (4.0_dp*Area)) * (b(i)*b(j) + c(i)*c(j))
                    call add_to_sparse(n(i), n(j), K_e)
                end do
            end do
        end do
    end subroutine assemble_poisson


    ! --- SOLVER GAUSS-SEIDEL (Conectado al Panel de Control) ---
    subroutine solve_gs(var)
        real(dp), intent(inout) :: var(:)
        integer :: it, i, k, j
        real(dp) :: sum_j, diag_val
        
        ! AQUI USAMOS LA VARIABLE DEL PANEL DE CONTROL
        do it = 1, max_gs
            do i = 1, n_nodes
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
                    else
                        var(i) = var(i)
                    end if
                end if
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
    real(dp) :: error

    call read_mesh()
    
    allocate(Psi(n_nodes), Xi(n_nodes), Psi_old(n_nodes))
    allocate(A_val(max_adj, n_nodes), A_col(max_adj, n_nodes))
    allocate(num_adj(n_nodes), b_vec(n_nodes))
    
    Psi = 0.0_dp
    Xi = 0.0_dp

    ! Condiciones iniciales usando la velocidad del panel de control
    do i = 1, n_nodes
        if (node_type(i) == 3) Psi(i) = V0 * y(i)      
        if (node_type(i) == 2) Psi(i) = V0 * 1.0_dp    
        if (node_type(i) == 1 .or. node_type(i) >= 5) Psi(i) = 0.0_dp 
    end do

    print*, "--------------------------------------------------------"
    print*, "Iniciando Simulación CFD (Galerkin P1 - Psi/Xi)"
    print*, "Nodos:", n_nodes, "| Elementos:", n_elements
    print*, "--------------------------------------------------------"

    ! --- BUCLE DE PICARD (Conectado al Panel de Control) ---
    do iter = 1, max_picard
        Psi_old = Psi

        call update_thom_bc()
        call assemble_vorticity_transport()
        call solve_gs(Xi)
        call assemble_poisson()
        call solve_gs(Psi)

        error = sqrt(sum((Psi - Psi_old)**2)/n_nodes)
        if (mod(iter, 20) == 0) print*, "Picard Iter:", iter, "| Error global:", error
        
        if (error > 1.0d3 .or. error /= error) then
            print*, "¡ERROR FATAL! Divergencia detectada en iteración", iter
            exit
        end if
        
        if (error < tol) then
            print*, " Convergió exitosamente en iteración: ", iter
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
    
    print*, "Resultados guardados en FLUJO y VORT."

end program NavierStokes_Galerkin

