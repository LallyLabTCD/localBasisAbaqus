C     Subroutine to calculate the consistent tangent matrix using the 
C     perturbation method outlined by Miehe and Sun
C
C     15-Aug-2018
C     David Nolan      
C          
      subroutine kCTM(XKIRCH1,DFGRD1,NTENS,PROPS,XAMAT,NPT,ITER,
     + NANISO,DDSDDE,NPROPS,XJ1,NSHR)
C      
      INCLUDE 'ABA_PARAM.INC'  
C      
      DIMENSION XKIRCH1(3,3), DFP(3,3), DFGRD_PERT(3,3),XKIRCH_PERT(3,3)
     1 , CMJ(3,3), CMJVEC(NTENS), ilist(6), jlist(6), XAMAT(3,2), XI4(2)
     2 , DFGRD1(3,3), DDSDDE(NTENS,NTENS), PROPS(NPROPS)
C
C      
C LOCAL ARRAYS
C ----------------------------------------------------------------
C XKIRCH1         - True Kirchhoff stress (i.e. not based on perturbed def. grad.)
C DFGRD1          - True deformation gradient for this increment      
C DFP             - Increment of the perturbed deformation gradient
C DFGRD_PERT      - Perturbed def. grad.
C XKIRCH_PERT     - Perturbed Kirchhoff stress
C CMJ             - (:,:,I,J) Components of material jacobian
C CMJVEC          - Above in vector form
C ILIST, JLIST    - Set of index components to be perturbed upon
C DUMSTRSS        - DUMMY STRESS TENSOR
C XAMAT           - array containing reference configuration fibre vectors
C XI4             - array storing the fibre invariant (will not be used in this subroutine)
C DDSDDE          - Voigt notation matrix to store the tangent matrix 
C PROPS           - Properties read in from Abaqus .inp file
C ----------------------------------------------------------------
C          
C     Perturbation parameter
      eps = 1.0e-08     
C
C     Material properties      
      C10 = PROPS(1)
      D1 = PROPS(2)
      xk1 = PROPS(3)
      xk2 = PROPS(4)
      xkap = PROPS(5)          
C                
C     The Kirchhoff stress is perturbed six times over the
C     i, j'th components of the deformation gradient. This
C     array lists each of the componets to be perturbed over.
C         i        |      j      
      ilist(1) = 1; jlist(1) = 1;
      ilist(2) = 2; jlist(2) = 2;
      ilist(3) = 3; jlist(3) = 3;
      ilist(4) = 1; jlist(4) = 2;
      ilist(5) = 1; jlist(5) = 3;
      ilist(6) = 2; jlist(6) = 3;        
C      
C     Loop over each of the six components of the deformation gradient
C     which are to be perturbed, each time calculating a new column of
C     the tangent matrix.      
      Perturbation: do iter = 1,NTENS          
C      
C-----------------------------------------------------------------
C         CREATE THE PERTURBED DEFORMATION GRADIENT
C-----------------------------------------------------------------     
C         Pick out the (i,j) component of the deformation gradient
C         to be perturbed      
          ii = ilist(iter)
          jj = jlist(iter)     
C     
C         Calculate the increment of the perturbed deformation gradient          
          call kdelF(ii, jj, DFGRD1, eps, DFP)
C          
C         Add to the "true" deformation gradient to create the perturbed
C         deformation gradient          
          DFGRD_PERT = DFGRD1 + DFP
C
C-----------------------------------------------------------------
C         CALCULATE KIRCHHOFF STRESS BASED ON THE PERTURBED DEFORMATION GRADIENT
C-----------------------------------------------------------------
C                             N.B.
          call kstress_calc(DFGRD_PERT, C10, D1, xkap, xk1, xk2, 
     1    XKIRCH_PERT, XJP, XAMAT, XI1, XI2, XI4, NPT, NANISO,
     2    iter)                    
C
C-----------------------------------------------------------------
C          
C         Difference between the perturbed(i,j) and unpert. stress
          CMJ = XKIRCH_PERT - XKIRCH1          
C         
C         Normalize by perturbation parameter and return to Cauchy stress       
          CMJ = CMJ/XJ1/eps
C      
C      
C-----------------------------------------------------------------
C         FORM THE TANGENT MATRIX
C-----------------------------------------------------------------      
C
C         Convert tensor to Voigt vector          
          call kmatrix2vector(CMJ, CMJVEC, NSHR)
C
C         Insert this vector into the iter'th column of DDSDDE          
          do insert = 1,NTENS
C          
              DDSDDE(insert,iter) = CMJVEC(insert)
C            
          end do
C		  
C
       end do Perturbation
       
       end subroutine kCTM