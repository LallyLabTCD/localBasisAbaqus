C     Subroutine to perform a polar decomposition on the deformation gradient
C
C     Adapted from Sommer for an Abaqus subroutine. Downloaded from the
C     International Society for Biomechanics' website. This subroutine is
C     free for non-commercial use.      
C     https://isbweb.org/software/movanal/sommer.txt      
C
C     H.J. Sommer III, Professor of Mechanical Engineering, 207 Reber Building
C     The Pennsylvania State University, University Park, PA  16802
C     (814)865-2519, FAX (814)863-4848, Bitnet HJS@PSUECL, Internet HJS@ECL.PSU.EDU    
C      
c********1*********2*********3*********4*********5*********6*********7**
c
c	kpolarDecomp(G,R,S)
c
c	positive polar decomposition of 3x3 matrix  G = R * S
c	forcing positive orthonormal rotation matrix
c
c	INPUTS
c	G = 3x3 general matrix
c
c	OUTPUTS
c	R = 3x3 positive orthonormal matrix
c	S = 3x3 symmetric matrix
c
c	PRECISION:	single
c	COMMONS:	none
c	CALLS:		none
c	FUNCTIONS:	ABS, SQRT
c	REFERENCE:	Veldpaus, F.E., H.J. Woltring, and L.J.M.G. Dortmans,
c				A Least-Squares Algorithm for the Equiform
c				Transformation from Spatial Marker Coordinates, 
c				J. Biomechanics, 21(1):45-54 (1988).
c	DATE:		10/8/92 - HJSIII
c
c
      SUBROUTINE kpolarDecomp(G,S,R)
c
c	declarations
      INCLUDE 'ABA_PARAM.INC'
      dimension G(3,3),R(3,3),S(3,3),COG(3,3),XP(3,3),ADP(3,3),
     1 XPBI(3,3)
c
c	constants
      EPS=1.0E-5
c
c	cofactors and determinant of g
      COG(1,1)=G(2,2)*G(3,3)-G(2,3)*G(3,2)
      COG(2,1)=G(1,3)*G(3,2)-G(1,2)*G(3,3)
      COG(3,1)=G(1,2)*G(2,3)-G(1,3)*G(2,2)
      COG(1,2)=G(2,3)*G(3,1)-G(2,1)*G(3,3)
      COG(2,2)=G(1,1)*G(3,3)-G(1,3)*G(3,1)
      COG(3,2)=G(1,3)*G(2,1)-G(1,1)*G(2,3)
      COG(1,3)=G(2,1)*G(3,2)-G(2,2)*G(3,1)
      COG(2,3)=G(1,2)*G(3,1)-G(1,1)*G(3,2)
      COG(3,3)=G(1,1)*G(2,2)-G(1,2)*G(2,1)
      G3=G(1,1)*COG(1,1)+G(2,1)*COG(2,1)+G(3,1)*COG(3,1)
c
c	P = trans(G) * G = S * S
      DO 10000 I=1,3
      XP(I,1)=G(1,I)*G(1,1)+G(2,I)*G(2,1)+G(3,I)*G(3,1)
      XP(I,2)=G(1,I)*G(1,2)+G(2,I)*G(2,2)+G(3,I)*G(3,2)
      XP(I,3)=G(1,I)*G(1,3)+G(2,I)*G(2,3)+G(3,I)*G(3,3)
10000 CONTINUE
c
c	adjoint of P
      ADP(1,1)=XP(2,2)*XP(3,3)-XP(2,3)*XP(3,2)
      ADP(2,2)=XP(1,1)*XP(3,3)-XP(1,3)*XP(3,1)
      ADP(3,3)=XP(1,1)*XP(2,2)-XP(1,2)*XP(2,1)
c
c	G invariants
      G1SQ=XP(1,1)+XP(2,2)+XP(3,3)
      G1=SQRT(G1SQ)
      G2SQ=ADP(1,1)+ADP(2,2)+ADP(3,3)
      G2=SQRT(G2SQ)
c
c	initialize iteration
      H1=G2/G1SQ
      H2=G3*G1/G2SQ
      X=1.0
      Y=1.0
c
c	iteration loop
10001 CONTINUE
      DEN=2.0*(X*Y-H1*H2)
      RES1=1.0-X*X+2.0*H1*Y
      RES2=1.0-Y*Y+2.0*H2*X
      DX=(Y*RES1+H1*RES2)/DEN
      DY=(H2*RES1+X*RES2)/DEN
      X=X+DX
      Y=Y+DY
      IF(ABS(DX/X).GT.EPS.OR.ABS(DY/Y).GT.EPS)GO TO 10001
c
c	BETA invariants
      BETA1=X*G1
      BETA2=Y*G2
c
c	invert ( trans(G) * G + BETA2 * identity )
      XP(1,1)=XP(1,1)+BETA2
      XP(2,2)=XP(2,2)+BETA2
      XP(3,3)=XP(3,3)+BETA2
      XPBI(1,1)=XP(2,2)*XP(3,3)-XP(2,3)*XP(3,2)
      XPBI(1,2)=XP(1,3)*XP(3,2)-XP(1,2)*XP(3,3)
      XPBI(1,3)=XP(1,2)*XP(2,3)-XP(1,3)*XP(2,2)
      XPBI(2,1)=XP(2,3)*XP(3,1)-XP(2,1)*XP(3,3)
      XPBI(2,2)=XP(1,1)*XP(3,3)-XP(1,3)*XP(3,1)
      XPBI(2,3)=XP(1,3)*XP(2,1)-XP(1,1)*XP(2,3)
      XPBI(3,1)=XP(2,1)*XP(3,2)-XP(2,2)*XP(3,1)
      XPBI(3,2)=XP(1,2)*XP(3,1)-XP(1,1)*XP(3,2)
      XPBI(3,3)=XP(1,1)*XP(2,2)-XP(1,2)*XP(2,1)
      DETPBI=XP(1,1)*XPBI(1,1)+XP(2,1)*XPBI(1,2)+XP(3,1)*XPBI(1,3)
c
c	R = (cofac(G)+BETA1*G) * inv(trans(G)*G+BETA2*identity)
      DO 10002 I=1,3
      R(I,1)=((COG(I,1)+BETA1*G(I,1))*XPBI(1,1)
     1       +(COG(I,2)+BETA1*G(I,2))*XPBI(2,1)
     2       +(COG(I,3)+BETA1*G(I,3))*XPBI(3,1))/DETPBI
      R(I,2)=((COG(I,1)+BETA1*G(I,1))*XPBI(1,2)
     1       +(COG(I,2)+BETA1*G(I,2))*XPBI(2,2)
     2       +(COG(I,3)+BETA1*G(I,3))*XPBI(3,2))/DETPBI
      R(I,3)=((COG(I,1)+BETA1*G(I,1))*XPBI(1,3)
     1       +(COG(I,2)+BETA1*G(I,2))*XPBI(2,3)
     2       +(COG(I,3)+BETA1*G(I,3))*XPBI(3,3))/DETPBI
10002 CONTINUE
c
c	S = trans(R) * G
      DO 10003 I=1,3
      S(I,1)=R(1,I)*G(1,1)+R(2,I)*G(2,1)+R(3,I)*G(3,1)
      S(I,2)=R(1,I)*G(1,2)+R(2,I)*G(2,2)+R(3,I)*G(3,2)
      S(I,3)=R(1,I)*G(1,3)+R(2,I)*G(2,3)+R(3,I)*G(3,3)
10003 CONTINUE
c
c	done
      RETURN
      END