C     Subroutine to calculate the Kirchhoff stress given by the fibre reinforced 
C     modified anisotropic (MA) hyperelastic constitutive law.
C
C     15-Aug-2018
C     David Nolan      
C      
      SUBROUTINE kstress_calc(DGRAD, C10, D1, xkap, xk1, xk2, XKIRCH, 
     1 XJ, XAMAT, XI1, XI2, XI4, NPT, NANISO, iter)
C     
      INCLUDE 'ABA_PARAM.INC'
C      
      INTENT(IN) :: DGRAD, C10, D1, xkap, xk1, xk2, XAMAT, NANISO, iter
      INTENT(OUT):: XKIRCH, XJ, XI1, XI2, XI4      
C      
      DIMENSION DGRAD(3,3), BMAT(3,3), XKIRCH(3,3), XAMAT(3,2),
     1 B2MAT(3,3), XA(3), a_curnt(3), aoa(3,3), XANISOK(3,3),
     2 XANISOTOT(3,3), XI4(2)
C      
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)
C
C LOCAL ARRAYS
C ----------------------------------------------------------------
C DGRAD     - the deformation gradient
C BMAT      - left Cauchy-Green deformation tensor
C XKIRCH    - Kirchhoff stress
C XAMAT     - array containing fibre vectors
C B2MAT     - square of the B tensor
C XA        - fibre vector in the reference configuration
C a_curnt   - fibre vector in the current configuration
C aoa       - fibre structural tensor
C XANISOK   - anisotropic Kirchhoff stress for a given fibre family
C XANISOTOT - total anisotropic Kirchhoff stress
C XI4       - array storing the fibre invariant      
C ----------------------------------------------------------------      
C     
C      
C     Determinant of the deformation gradient
      XJ=DGRAD(1, 1)*DGRAD(2, 2)*DGRAD(3, 3)
     1   -DGRAD(1, 2)*DGRAD(2, 1)*DGRAD(3, 3)
     2   +DGRAD(1, 2)*DGRAD(2, 3)*DGRAD(3, 1)
     3   +DGRAD(1, 3)*DGRAD(3, 2)*DGRAD(2, 1)
     4   -DGRAD(1, 3)*DGRAD(3, 1)*DGRAD(2, 2)
     5   -DGRAD(2, 3)*DGRAD(3, 2)*DGRAD(1, 1)
C     
C        
C     Calculate left Cauchy-Green deformation tensor     
      BMAT = MATMUL(DGRAD,TRANSPOSE(DGRAD))
C
C       
C-----------------------------------------------------------------------
C     CALCULCATE THE INVARIANTS                
C-----------------------------------------------------------------------
C     
C     I1 invariant      
      XI1 = BMAT(1,1)+BMAT(2,2)+BMAT(3,3)
C        
C     Square of the B matrix
      B2MAT = MATMUL(BMAT,BMAT)
C
C     Trace of B squared   
      TRB2 = B2MAT(1,1)+B2MAT(2,2)+B2MAT(3,3)
C     
C     I2 invariant      
      XI2 = 0.5*((XI1**2.0) - TRB2)
C
C     I3 invariant      
      XI3=BMAT(1, 1)*BMAT(2, 2)*BMAT(3, 3)
     1   -BMAT(1, 2)*BMAT(2, 1)*BMAT(3, 3)
     2   +BMAT(1, 2)*BMAT(2, 3)*BMAT(3, 1)
     3   +BMAT(1, 3)*BMAT(3, 2)*BMAT(2, 1)
     4   -BMAT(1, 3)*BMAT(3, 1)*BMAT(2, 2)
     5   -BMAT(2, 3)*BMAT(3, 2)*BMAT(1, 1)
C                                   
C      
C--------------------------------------------------------------------
C     CALCULATE THE ISOTROPIC PORTION OF THE KIRCHHOFF STRESS
C--------------------------------------------------------------------
C
C     PART 1
      COEFF1 = TWO*C10/(XJ**(TWO/THREE))
C
C     KIRCHHOFF STRESS PART 2
C
      TRBMAT= COEFF1*(BMAT(1,1)+BMAT(2,2)+BMAT(3,3))/THREE
C
C     KIRCHHOFF STRESS PART 1
C      
      XKIRCH = COEFF1 * BMAT
C
C     SUBTRACT THE PART 2
C    
      DO I = 1,3
        XKIRCH(I,I) = XKIRCH(I,I) - TRBMAT
      END DO
C
C     FORM VOLUMETRIC PART, 3 
C      
      COEFF3 = 2*(XJ-ONE)*XJ/D1
C
C     ADD TO THE PREVIOUS PARTS
C      
      DO I = 1,3
        XKIRCH(I,I) = XKIRCH(I,I) + COEFF3
      END DO
C           
C      
C      
C--------------------------------------------------------------------
C     CALCULATE THE ANISOTROPIC PORTION OF THE KIRCH STRESS
C--------------------------------------------------------------------       
C  
C     Zero the matrix      
      DO I=1,3
        DO J=1,3
          XANISOTOT(I,J) = ZERO
        END DO
      END DO
C                  
C--------------------------------------------------------------------
C     Loop over the number of fibre families    
C--------------------------------------------------------------------        
      DO IANISO = 1,NANISO                                     
C         
C         Fibre vector in reference configuration  
          XA = XAMAT(1:3,IANISO);
C                                              
C         Calculate fibre vector in current configuration: a = F*A
          a_curnt = MATMUL(DGRAD,XA)
C
C         Calculate structural tensor: a \otimes a     
          DO i=1,3
              DO j=1,3
                aoa(i,j) = a_curnt(i) * a_curnt(j)
              ENDDO
          ENDDO            
C            
C         Calculate the fibre invariant
          XI4(IANISO) = aoa(1,1) + aoa(2,2) + aoa(3,3)
C            
C         The fibres only contribute to stress when in tension
          IF (XI4(IANISO)>ONE) THEN
C          
C             Calculate the various parts of the equation for kirch stress
              p1 = 2.0*xk1;                      
C
              p2 = XI4(IANISO) - 1.0
C
              p3 = xk2*(p2**TWO)                       
C
              p4 = exp(p3)
C                
              XANISOK = p1*p2*p4*aoa
C                
C             Sum the stress over the fibre families              
              XANISOTOT = XANISOTOT + XANISOK     
C            
          END IF
C        
      END DO
C     
C     Add the isotropic and anisotropic parts      
      XKIRCH = XKIRCH + XANISOTOT
C        
      RETURN
      END SUBROUTINE kstress_calc