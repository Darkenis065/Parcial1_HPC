C=======================================================================
C GUIA DE COMPILACION Y ENTORNO HPC (OPENMP WAVEFRONT TASKING)
C=======================================================================
C Compilador GNU (GCC / gfortran):
C   gfortran -O3 -march=native -fopenmp -ffast-math -funroll-loops navier.f -o sim
C
C Compilador Intel (ifort / ifx):
C   ifort -O3 -xHost -qopenmp -ipo -fp-model fast=2 navier.f -o sim
C
C Variables de Entorno Recomendadas para Maxima Eficiencia:
C   export OMP_NUM_THREADS=$(nproc)    ! Utilizar todos los nucleos
C   export OMP_PROC_BIND=true          ! Fijar hilos a nucleos (afinidad)
C   export OMP_PLACES=cores            ! Evitar migracion de hilos
C=======================================================================

      PROGRAM NAST2D
      USE omp_lib
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NX=140, NY=140)
      DIMENSION P(NX,NY), XSI(NX,NY)
      INTEGER I, J, NITER, JJJ, LENGTH, J1, I1, I2

      JJJ = 10000

C     POSICION DEL OBSTACULO
      LENGTH = NY/4
      J1 = NY/3
      I1 = NX/3
      I2 = I1 + LENGTH - 1

C     PARAMETROS DE SIMULACION
      RE = 2.D0
      VW = 0.1D0
      SW = VW

      OPEN (33, FILE='FLUJO')
      OPEN (34, FILE='VORT')

C     INICIALIZACION EMBARRASSINGLY PARALLEL
!$omp parallel do collapse(2) private(I, J)
      DO J = 1, NY
         DO I = 1, NX
            XSI(I,J) = 0.D0
            P(I,J) = (FLOAT(J) - 1.D0)
         END DO
      END DO

C     RELAJACION ITERATIVA
      DO NITER = 1, JJJ
         WRITE (6,*) NITER
         CALL PRELAX(P, XSI, NX, NY, SW, J1, I1, I2)
         CALL VRELAX(P, XSI, NX, NY, VW, J1, I1, I2, RE)
      END DO

C     ESCRITURA DE RESULTADOS (Secuencial para evitar condiciones de carrera I/O)
      DO I = 1, NX
         DO J = 1, NY
            WRITE (33, '(2I10, F18.13)') I, J, P(I,J)
            WRITE (34, '(2I10, F18.13)') I, J, XSI(I,J)
         END DO
         WRITE (33, '()')
         WRITE (34, '()')
      END DO

      CLOSE (33)
      CLOSE (34)
      END

C=======================================================================
      SUBROUTINE PRELAX(P, XSI, NX, NY, SW, J1, I1, I2)
      USE omp_lib
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(NX,NY), XSI(NX,NY)
      
      INTEGER BS, JJ, II, J_MAX, I_MAX, I, J
      PARAMETER (BS = 32)

      CSW_A = SW * 0.25D0
      CSW_B = 1.D0 - SW

C     FRONTERA G (Superior) - PARALELIZACION DIRECTA
!$omp parallel do private(I, TEMP)
      DO I = 2, NX-1
         TEMP = 1.D0 + 2.D0*P(I,NY-1) + P(I+1,NY) + P(I-1,NY)
         P(I,NY) = CSW_A*TEMP + CSW_B*P(I,NY)
      END DO

C     FRONTERA F (Izquierda) - PARALELIZACION DIRECTA
!$omp parallel do private(J, TEMP)
      DO J = 2, NY-1
         TEMP = 2.D0*P(2,J) + P(1,J+1) + P(1,J-1) - XSI(1,J)
         P(1,J) = CSW_A*TEMP + CSW_B*P(1,J)
      END DO

C     FRONTERA H (Derecha) - PARALELIZACION DIRECTA
!$omp parallel do private(J, TEMP)
      DO J = 2, NY-1
         TEMP = 2.D0*P(NX-1,J) + P(NX,J-1) + P(NX,J+1) - XSI(NX,J)
         P(NX,J) = CSW_A*TEMP + CSW_B*P(NX,J)
      END DO

C     FRONTERA E, A (Inferior)
!$omp parallel do private(I)
      DO I = 1, NX
         P(I,1) = 0.D0
      END DO

C     PUNTOS INTERIORES: OPENMP WAVEFRONT TASKING
!$omp parallel
!$omp single
      DO JJ = 2, NY-1, BS
         DO II = 2, NX-1, BS
C           Se asegura que los indices negativos referencien a la frontera (memoria segura y satisfecha)
!$omp task depend(in: P(MAX(1,II-BS),JJ), P(II,MAX(1,JJ-BS)))
!$omp&     depend(out: P(II,JJ))
!$omp&     private(I, J, TEMP, I_MAX, J_MAX) 
!$omp&     firstprivate(II, JJ, CSW_A, CSW_B)
            
            J_MAX = MIN(JJ + BS - 1, NY-1)
            I_MAX = MIN(II + BS - 1, NX-1)
            
            DO J = JJ, J_MAX
               DO I = II, I_MAX
                  TEMP = P(I,J-1) + P(I,J+1) + P(I-1,J) + P(I+1,J) 
     &                 - XSI(I,J)
                  P(I,J) = CSW_A*TEMP + CSW_B*P(I,J)
               END DO
            END DO
!$omp end task
         END DO
      END DO
!$omp end single
!$omp end parallel

C     FRONTERA DEL OBSTACULO (D, C, B)
!$omp parallel do collapse(2) private(I, J)
      DO J = 1, J1
         DO I = I1, I2
            P(I,J) = 0.D0
         END DO
      END DO

C     ESQUINAS
      P(1,NY) = (P(2,NY) + P(1,NY-1)) * 0.5D0
      P(NX,NY) = (P(NX-1,NY) + P(NX,NY-1)) * 0.5D0

      RETURN
      END

C=======================================================================
      SUBROUTINE VRELAX(P, XSI, NX, NY, VW, J1, I1, I2, RE)
      USE omp_lib
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(NX,NY), XSI(NX,NY)
      
      INTEGER BS, JJ, II, J_MAX, I_MAX, I, J
      PARAMETER (BS = 32)

      CVW_A = VW * 0.25D0
      CVW_B = 1.D0 - VW
      CRE4  = RE * 0.25D0

C     FRONTERA H (Derecha)
!$omp parallel do private(J, TEMP)
      DO J = 2, NY-1
         TEMP = 2.D0*XSI(NX-1,J) + XSI(NX,J+1) + XSI(NX,J-1)
         XSI(NX,J) = CVW_A*TEMP + CVW_B*XSI(NX,J)
      END DO

C     PUNTOS INTERIORES: OPENMP WAVEFRONT TASKING
!$omp parallel
!$omp single
      DO JJ = 2, NY-1, BS
         DO II = 2, NX-1, BS
C           La dependencia recae estrictamente sobre la matriz que se actualiza iterativamente (XSI)
!$omp task depend(in: XSI(MAX(1,II-BS),JJ), XSI(II,MAX(1,JJ-BS)))
!$omp&     depend(out: XSI(II,JJ))
!$omp&     private(I, J, TEMP, TEMP2, TEMP3, I_MAX, J_MAX) 
!$omp&     firstprivate(II, JJ, CVW_A, CVW_B, CRE4)
            
            J_MAX = MIN(JJ + BS - 1, NY-1)
            I_MAX = MIN(II + BS - 1, NX-1)
            
            DO J = JJ, J_MAX
               DO I = II, I_MAX
                  TEMP = XSI(I,J+1) + XSI(I,J-1) 
     &                 + XSI(I+1,J) + XSI(I-1,J)
                  TEMP2 = CRE4 * (P(I,J+1) - P(I,J-1))
     &                  * (XSI(I+1,J) - XSI(I-1,J))
                  TEMP3 = CRE4 * (P(I+1,J) - P(I-1,J))
     &                  * (XSI(I,J+1) - XSI(I,J-1))
                  XSI(I,J) = CVW_A*(TEMP - TEMP2 + TEMP3) 
     &                     + CVW_B*XSI(I,J)
               END DO
            END DO
!$omp end task
         END DO
      END DO
!$omp end single
!$omp end parallel

C     FRONTERAS G, E, A
!$omp parallel do private(I)
      DO I = 1, NX
         XSI(I,1) = 0.D0
         XSI(I,NY) = 0.D0
      END DO

C     FRONTERA F (Izquierda)
!$omp parallel do private(J)
      DO J = 2, NY-1
         XSI(1,J) = 0.D0
      END DO

C     CONDICION DE FRONTERA EN LOS BORDES DEL OBSTACULO
!$omp parallel do private(J)
      DO J = 1, J1
         XSI(I1,J) = 2.D0*(P(I1-1,J) - P(I1,J))
         XSI(I2,J) = 2.D0*(P(I2+1,J) - P(I2,J))
      END DO

!$omp parallel do private(I)
      DO I = I1+1, I2-1
         XSI(I,J1) = 2.D0*(P(I,J1+1) - P(I,J1))
      END DO

C     INTERIOR DEL OBSTACULO
!$omp parallel do collapse(2) private(I, J)
      DO J = 1, J1-1
         DO I = I1+1, I2-1
            XSI(I,J) = 0.D0
         END DO
      END DO

      RETURN
      END
