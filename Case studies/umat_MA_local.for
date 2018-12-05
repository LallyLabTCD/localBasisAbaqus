C
C ----------------------------------------------------------------
C UMAT FOR COMPRESSIBLE HYPERELASTIC MODEL DEVELOPED BY NOLAN, 
C GOWER, DESTRADE, OGDEN, MCGARRY (2014)
C
C 3D CONTINUUM ELEMENTS
C ----------------------------------------------------------------
C      
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C ----------------------------------------------------------------
C
C LOCAL ARRAYS
C ----------------------------------------------------------------
C XKIRCH      - Kirchhoff stress
C DUMSTRSS    - Dummy cauchy stress tensor
C XAMAT       - Array of the unit vectors describing the fibre direction
C XIF         - array storing the fibre invariant  
C U, R, RT    - matrices from the polar decomposition of the def. grad.      
C ----------------------------------------------------------------
C
      DIMENSION XKIRCH(3,3), DUMSTRSS(3,3), XAMAT(3,2), XI4(2), R(3,3),
     1 U(3,3), RT(3,3), XADUM(3,1), XADUM1(3,1)
C           
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0,
     1 SIX=6.D0)
C
C    ** PROPERTIES **
C
C PROPS(1) - C10
C PROPS(2) - D1
C PROPS(3) - K1 
C PROPS(4) - K2
C PROPS(5) - KAPPA (FIBRE DISPERSION)
C PROPS(6) - THETA1 (ORIENTATION ANGLE OF THE 1ST FIBRE FAMILY)
C PROPS(7) - THETA2 (ORIENTATION ANGLE OF THE 2ND FIBRE FAMILY)
C
C    ** STATE DEPENDENT VARIABLES **   
C
C STATEV(1) - I1 
C STATEV(2) - I2
C STATEV(3) - I3
C STATEV(4) - I4
C STATEV(5) - I6
C STATEV(6) - J
C
C ----------------------------------------------------------------
C     NUMBER OF FIBRE FAMILIES
C ----------------------------------------------------------------
C
      NANISO = 2        
C
C ----------------------------------------------------------------
C     MATERIAL PROPERTIES
C ----------------------------------------------------------------
      C10 = PROPS(1)
      D1 = PROPS(2)
      xk1 = PROPS(3)
      xk2 = PROPS(4)
      xkap = PROPS(5)   ! not used in this UMAT
      THETAD1 = PROPS(6)
      THETAD2 = PROPS(7)
        
C       Convert degrees to radians        
      XPI = 3.14159265359
      THETAR1 = THETAD1*XPI/180.0
      THETAR2 = THETAD2*XPI/180.0
C        
C     Fibre vectors in the reference configuration
      XAMAT(1,1)=COS(thetar1); XAMAT(2,1)= SIN(thetar1); XAMAT(3,1)=0.
      XAMAT(1,2)=COS(thetar2); XAMAT(2,2)= SIN(thetar2); XAMAT(3,2)=0.
C
C     This is related to the calculation of the numberical Jacobian
C     0 indicates that the stress calculation is not part of the 
C     stress calculation.        
      iter=0       
C
C---------------------------------------------------------------------------
C         POLAR DECOMPOSITION
C---------------------------------------------------------------------------            
      call kpolarDecomp(DFGRD1, U, R)
C     Calculate the transpose of the rotation matrix      
      RT = transpose(R)
      !if(NOEL.EQ.2276)then
      !    write(6,*) 'F ', NPT
      !    call kprinter(DFGRD1,3,3)
      !    write(6,*) 'R'
      !    call kprinter(R,3,3)
      !    write(6,*) '*********'
      !endif
C
C-----------------------------------------------------------------
C        CALCULATE LOCAL FIBRE VECTOR IN THE REFERENCE CONFIG.
C-----------------------------------------------------------------                 
C    
      DO IAN = 1,NANISO
C         Pick out fibre vector          
          XADUM(:,1) = XAMAT(:,IAN)
C          
C         Update to local basis system       
          XADUM1 = MATMUL(RT,XADUM)
C
C         Store result          
          XAMAT(:,IAN) = XADUM1(:,1)         
      END DO        
C	  
C-----------------------------------------------------------------
C     ZERO THE TANGENT MATRIX
C-----------------------------------------------------------------        
C
      DO I=1,NTENS
      DO J=1,NTENS
        DDSDDE(I,J) = ZERO
      END DO
      END DO          
C      
C-----------------------------------------------------------------
C     CALCULATE THE KIRCHHOFF STRESS 
C-----------------------------------------------------------------
C     
      call kstress_calc(DFGRD1, C10, D1, xkap, xk1, xk2, XKIRCH, XJ, 
     1 XAMAT, XI1, XI2, XI4, NPT, NANISO, iter)
C     
        STATEV(1) = XI1
        STATEV(2) = XI2
        STATEV(3) = XJ
        STATEV(4) = XI4(1)
        STATEV(5) = XI4(2)    
C             
C-----------------------------------------------------------------
C CONVERT KIRCHHOFF STRESS TO CAUCHY STRESS               
C-----------------------------------------------------------------      
      DUMSTRSS = XKIRCH / XJ           
C     
C     Convert the stress tensor to a Voigt vector      
      call kmatrix2vector(DUMSTRSS, STRESS, nshr)         
      
C      
C-----------------------------------------------------------------
C CALCULATE THE TANGENT MATRIX
C-----------------------------------------------------------------
C      
      call kCTM(XKIRCH,DFGRD1,NTENS,PROPS,XAMAT,NPT,ITER,
     + NANISO,DDSDDE,NPROPS,XJ,NSHR)      
C
      RETURN
C
C      
      ENDSUBROUTINE UMAT
C
C-----------------------------------------------------------------------------
C           SUBROUTINES
C-----------------------------------------------------------------------------
C
C     Subroutine to calculate the Kirchhoff stress      
      INCLUDE 'kstress_calc.for'
C      
C     Subroutine to perform a polar decomposition of the deformation gradient      
      INCLUDE 'kpolarDecomp.for'
C      
C     Subroutine to calculate the tangent matrix 
      INCLUDE 'kCTM.for'
C      
C     Helper subroutines      
      INCLUDE 'utilities.for'
      