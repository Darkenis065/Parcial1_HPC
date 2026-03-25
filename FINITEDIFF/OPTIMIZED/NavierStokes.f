      PROGRAM NAST2D
      IMPLICIT REAL*8 (A-H,O-Z)
C=======================================================================
C SOLUCION ESTACIONARIA DE LAS ECUACIONES DE NAVIER STOKES EN DOS 
C DIMENSIONES PARA UN FLUIDO VISCOSO INCOMPRESIBLE
C REFACTORIZACION HPC: COLUMN-MAJOR, TILING, STRENGTH REDUCTION Y CSE
C=======================================================================
      PARAMETER (NX=140, NY=140)
      DIMENSION P(NX,NY), XSI(NX,NY)
      INTEGER I, J, NITER, JJJ, LENGTH, J1, I1, I2

      WRITE (6,*) 'NUMERO DE ITERACIONES'
      READ (5,*) JJJ

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

C     INICIO: FLUJO (P) Y VORTICIDAD (XSI) MEDIANTE OPERACIONES NATIVAS
      XSI = 0.D0
      DO J = 1, NY
         DO I = 1, NX
            P(I,J) = (FLOAT(J) - 1.D0)
         END DO
      END DO

C     RELAJACION DE STREAM Y LUEGO VORTICIDAD
      DO NITER = 1, JJJ
         WRITE (6,*) NITER
         CALL PRELAX(P, XSI, NX, NY, SW, J1, I1, I2)
         CALL VRELAX(P, XSI, NX, NY, VW, J1, I1, I2, RE)
      END DO

C     ESCRITURA DE RESULTADOS
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
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(NX,NY), XSI(NX,NY)
      
C     DECLARACION EXPLICITA DE ENTEROS PARA EVITAR CONFLICTOS IMPLICITOS
      INTEGER BS, JJ, II, J_MAX, I_MAX, I, J
      PARAMETER (BS = 32)

C     ERRADICACION DE SUBEXPRESIONES COMUNES (CSE / LICM)
      CSW_A = SW * 0.25D0
      CSW_B = 1.D0 - SW

C     FRONTERA G (Superior)
      DO I = 2, NX-1
         TEMP = 1.D0 + 2.D0*P(I,NY-1) + P(I+1,NY) + P(I-1,NY)
         P(I,NY) = CSW_A*TEMP + CSW_B*P(I,NY)
      END DO

C     FRONTERA F (Izquierda)
      DO J = 2, NY-1
         TEMP = 2.D0*P(2,J) + P(1,J+1) + P(1,J-1) - XSI(1,J)
         P(1,J) = CSW_A*TEMP + CSW_B*P(1,J)
      END DO

C     FRONTERA H (Derecha)
      DO J = 2, NY-1
         TEMP = 2.D0*P(NX-1,J) + P(NX,J-1) + P(NX,J+1) - XSI(NX,J)
         P(NX,J) = CSW_A*TEMP + CSW_B*P(NX,J)
      END DO

C     FRONTERA E, A (Inferior)
      DO I = 1, NX
         P(I,1) = 0.D0
      END DO

C     PUNTOS INTERIORES: BLOCKING / TILING Y ORDEN COLUMN-MAJOR
      DO JJ = 2, NY-1, BS
         J_MAX = MIN(JJ + BS - 1, NY-1)
         DO II = 2, NX-1, BS
            I_MAX = MIN(II + BS - 1, NX-1)
            DO J = JJ, J_MAX
               DO I = II, I_MAX
                  TEMP = P(I,J-1) + P(I,J+1) + P(I-1,J) + P(I+1,J)
     &                 - XSI(I,J)
                  P(I,J) = CSW_A*TEMP + CSW_B*P(I,J)
               END DO
            END DO
         END DO
      END DO

C     FRONTERA DEL OBSTACULO (D, C, B)
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
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(NX,NY), XSI(NX,NY)
      
C     DECLARACION EXPLICITA DE ENTEROS PARA TILING
      INTEGER BS, JJ, II, J_MAX, I_MAX, I, J
      PARAMETER (BS = 32)

C     ERRADICACION DE SUBEXPRESIONES COMUNES (CSE / LICM)
      CVW_A = VW * 0.25D0
      CVW_B = 1.D0 - VW
      CRE4  = RE * 0.25D0

C     FRONTERA H (Derecha)
      DO J = 2, NY-1
         TEMP = 2.D0*XSI(NX-1,J) + XSI(NX,J+1) + XSI(NX,J-1)
         XSI(NX,J) = CVW_A*TEMP + CVW_B*XSI(NX,J)
      END DO

C     PUNTOS INTERIORES: BLOCKING / TILING Y ORDEN COLUMN-MAJOR
      DO JJ = 2, NY-1, BS
         J_MAX = MIN(JJ + BS - 1, NY-1)
         DO II = 2, NX-1, BS
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
         END DO
      END DO

C     FRONTERAS G, E, A
      DO I = 1, NX
         XSI(I,1) = 0.D0
         XSI(I,NY) = 0.D0
      END DO

C     FRONTERA F (Izquierda)
      DO J = 2, NY-1
         XSI(1,J) = 0.D0
      END DO

C     CONDICION DE FRONTERA EN LOS BORDES DEL OBSTACULO
      DO J = 1, J1
         XSI(I1,J) = 2.D0*(P(I1-1,J) - P(I1,J))
         XSI(I2,J) = 2.D0*(P(I2+1,J) - P(I2,J))
      END DO

      DO I = I1+1, I2-1
         XSI(I,J1) = 2.D0*(P(I,J1+1) - P(I,J1))
      END DO

C     INTERIOR DEL OBSTACULO
      DO J = 1, J1-1
         DO I = I1+1, I2-1
            XSI(I,J) = 0.D0
         END DO
      END DO

      RETURN
      END

