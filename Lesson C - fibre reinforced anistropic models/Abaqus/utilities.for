C * KDELF        -   Calculate the increment of the deformation gradient for
C                   a given perturbation in (i,j), with epsilon
C
C * KPRINTER     -   Print out a matrix of any size
C
C * KMTMS     -   Multiply two matrices
C
C * KMATRIX2VECTOR     -   Convert a 3x3 matrix to a 6x1 vector
C                   
C * KDOTPROD           - Dot product of two vectors
C
C------------------------------------------------------------      
      subroutine kdotprod(A, B, dotp, n)
      
      INCLUDE 'ABA_PARAM.INC'
      
      intent(in) :: A, B, n
      intent(out):: dotp      
      
      dimension A(n), B(n)
      dotp = 0.0
      
      do i = 1,n
        dotp = dotp + A(i)*B(i)
      end do
      
      end subroutine kdotprod    
C      
C-------------------------------------------------------      
C
C-------------------------------------------------------
C      
      subroutine kmatrix2vector(XMAT, VEC, NSHR)
      
      INCLUDE 'ABA_PARAM.INC'
      
      intent(in) :: XMAT, NSHR
      intent(out):: VEC
      
      dimension xmat(3,3), vec(6)
  
        do i=1,3
            vec(i) = xmat(i,i);
        end do
               
        vec(4) = xmat(1,2);
        
        IF (NSHR==3) then
            vec(5) = xmat(1,3);
            vec(6) = xmat(2,3);
        END IF
      
      end subroutine kmatrix2vector      
C
C-------------------------------------------------------      
C
C-------------------------------------------------------
C      
      SUBROUTINE KMTMS (M, N, L, A, KA, B, KB, C, KC)
      
      INCLUDE 'ABA_PARAM.INC'
C      
      intent(in) :: M, N, L, A, KA, B, KB, KC
      intent(out):: C      
C      
C
C    PRODUCT OF REAL MATRICES
C
      DIMENSION A(KA,N), B(KB,L), C(KC,L)
      DOUBLE PRECISION W
C       
C
      DO 30 J = 1,L
         DO 20 I = 1,M
            W = 0.D0
            DO 10 K = 1,N
               W = W + A(I,K) * B(K,J)
   10       CONTINUE
            C(I,J) = W
   20    CONTINUE
   30 CONTINUE
      RETURN
      END SUBROUTINE KMTMS      
C
C-------------------------------------------------------      
C
C-------------------------------------------------------
C      
      subroutine kprinter(tens, m, n)
      
      INCLUDE 'ABA_PARAM.INC'
      
      intent(in):: tens, m, n      
      
      dimension tens(m,n)
        
        write(6,*)
        do i = 1,m
        do j = 1,n
            write(6,'(e19.9)',advance='no'),tens(i,j)
        end do
        write(6,*)
        end do
        write(6,*)
      return
      end subroutine kprinter      
C
C-------------------------------------------------------      
C
C-------------------------------------------------------
C  
      subroutine kdelF(m, n, DGRAD, eps, DF)
      
      INCLUDE 'ABA_PARAM.INC'
      
      intent (in) :: DGRAD, eps, m, n
      intent (out):: DF      
        
C Input: the index's i & j; The current deformation gradient (DGRAD). The perturbation
C        increment (eps)
C
C Output: The perturbed increment DF
C
      dimension dyad1(3,3), dyad2(3,3), DGRAD(3,3), DF(3,3), DFp1(3,3)
        
c Zero the dyad matrices
c
      do i = 1,3
        do j = 1,3
            dyad1(i,j) = zero
            dyad2(i,j) = zero
        end do
      end do
      
c Place the 1's in the correct location        
      dyad1(m,n) = 1.0;
c
      dyad2(n,m) = 1.0;
c      
c      KMTMS (M, N, L, A, KA, B, KB, C, KC)  
      call KMTMS(3, 3, 3, dyad1, 3, DGRAD, 3, DFp1, 3)
      DF = DFp1
      
      
      call KMTMS(3, 3, 3, dyad2, 3, DGRAD, 3, DFp1, 3)
      DF = DF + DFp1
           
      
      DF = 0.5*DF*eps            
      
      end subroutine kdelF      
C
C-------------------------------------------------------      
C
C-------------------------------------------------------
C       