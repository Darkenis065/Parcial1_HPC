      PROGRAM NAST2D
      IMPLICIT REAL*8 (A-H,O-Z)
C=======================================================================
C SOLUCION ESTACIONARIA DE LAS ECUACIONES DE NAVIER STOKES EN DOS 
C DIMENSIONES PARA UN FLUIDO VISCOSO INCOMPRESIBLE
C=======================================================================
C NX      : NUMERO DE CELDAS EN LA DIRECCION X
C NY      : NUMERO DE CELDAS EN LA DIRECCION Y
C P       : ARREGLO DE DIMENSION NX x NY QUE REPRESENTA EL FLUJO (STREAM)
C XSI     : ARREGLO DE DIMENSION NX x NY QUE REPRESENTA LA VORTICIDAD
C JJJ     : NUMERO DE ITERACIONES
C VW      : PARAMETRO DE RELAJACION DE LA VORTICIDAD
C SW      : PARAMETRO DE RELAJACION DEL FLUJO
C LENGTH  : ANCHO DEL OBSTACULO
C J1      : LARGO DEL OBSTACULO (ALTURA)
C RE      : NUMERO DE REYNOLDS DE RED
C=======================================================================

      PARAMETER (NX=140, NY=140)
      DIMENSION P(NX,NY), XSI(NX,NY)

      WRITE (6,*) 'NUMERO DE ITERACIONES'
      READ (5,*) JJJ

C     POSICION DEL OBSTACULO
      LENGTH = NY/4  ! ANCHO
      J1 = NY/3      ! ALTURA DEL OBSTACULO
      I1 = NX/3
      I2 = I1 + LENGTH - 1

C     PARAMETROS DE SIMULACION
      RE = 2.D0
      VW = 0.1D0
      SW = VW

      OPEN (33, FILE='FLUJO')
      OPEN (34, FILE='VORT')

C     INICIO: FLUJO (P) Y VORTICIDAD (XSI) CON VALORES ARBITRARIOS
      DO I = 1, NX
         DO J = 1, NY
            XSI(I,J) = 0.D0
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
C     RELAJACION DE LA FUNCION STREAM (P)
      DIMENSION P(NX,NY), XSI(NX,NY)

C     FRONTERA G (Superior)
      DO I = 2, NX-1
         TEMP = 1.D0 + 2.D0*P(I,NY-1) + P(I+1,NY) + P(I-1,NY)
         P(I,NY) = SW*(TEMP/4.D0) + (1.D0-SW)*P(I,NY)
      END DO

C     FRONTERA F (Izquierda)
      DO J = 2, NY-1
         TEMP = 2.D0*P(2,J) + P(1,J+1) + P(1,J-1) - XSI(1,J)
         P(1,J) = SW*(TEMP/4.D0) + (1.D0-SW)*P(1,J)
      END DO

C     FRONTERA H (Derecha)
      DO J = 2, NY-1
         TEMP = 2.D0*P(NX-1,J) + P(NX,J-1) + P(NX,J+1) - XSI(NX,J)
         P(NX,J) = SW*(TEMP/4.D0) + (1.D0-SW)*P(NX,J)
      END DO

C     FRONTERA E, A (Inferior)
      DO I = 1, NX
         P(I,1) = 0.D0
      END DO

C     PUNTOS INTERIORES
      DO I = 2, NX-1
         DO J = 2, NY-1
            TEMP = P(I,J-1) + P(I,J+1) + P(I-1,J) + P(I+1,J) - XSI(I,J)
            P(I,J) = SW*(TEMP/4.D0) + (1.D0-SW)*P(I,J)
         END DO
      END DO

C     FRONTERA DEL OBSTACULO (D, C, B)
      DO I = I1, I2
         DO J = 1, J1
            P(I,J) = 0.D0
         END DO
      END DO

C     ESQUINAS
      P(1,NY) = (P(2,NY) + P(1,NY-1)) / 2.D0
      P(NX,NY) = (P(NX-1,NY) + P(NX,NY-1)) / 2.D0

      RETURN
      END

C=======================================================================
      SUBROUTINE VRELAX(P, XSI, NX, NY, VW, J1, I1, I2, RE)
      IMPLICIT REAL*8 (A-H,O-Z)
C     RELAJA LA VORTICIDAD (XSI)
      DIMENSION P(NX,NY), XSI(NX,NY)
      RE4 = RE/4.D0

C     FRONTERA H (Derecha)
      DO J = 2, NY-1
         TEMP = 2.D0*XSI(NX-1,J) + XSI(NX,J+1) + XSI(NX,J-1)
         XSI(NX,J) = VW*(TEMP/4.D0) + (1.D0-VW)*XSI(NX,J)
      END DO

C     PUNTOS INTERIORES
      DO I = 2, NX-1
      DO J = 2, NY-1
      TEMP = XSI(I,J+1) + XSI(I,J-1) + XSI(I+1,J) + XSI(I-1,J)
      TEMP2 = RE4*(P(I,J+1) - P(I,J-1))*(XSI(I+1,J) - XSI(I-1,J))
      TEMP3 = RE4*(P(I+1,J) - P(I-1,J))*(XSI(I,J+1) - XSI(I,J-1))
      XSI(I,J) = VW*((TEMP - TEMP2 + TEMP3)/4.D0) + (1.D0-VW)*XSI(I,J)
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
      DO I = I1+1, I2-1
         DO J = 1, J1-1
            XSI(I,J) = 0.D0
         END DO
      END DO

      RETURN
      END
