MODULE Utilities
USE Global
IMPLICIT NONE
!
CONTAINS
!
FUNCTION Norm(g) RESULT(normresult)
USE Global
IMPLICIT NONE
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN) :: g
REAL(KIND=r8) :: normresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: tmp
!
mx(1) = SIZE(g,1); mx(2) = SIZE(g,2)
tmp = REAL(mx(1)*mx(2)*numvars,KIND=r8)
!
normresult = SQRT(SUM(g*g)/tmp)
!
END FUNCTION Norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION UDiv2D(f1,f2) RESULT(udivresult)
USE Global
IMPLICIT NONE
!
! UNDIVIDED divergence of the 2D flux function 
!
!  f = (f1(0:mx(1),1:mx(2)),f2(1:mx(1),0:mx(2))).  
!
! The result is stored in udivresult(1:mx(1),1:mx(2)).
!
REAL(KIND=r8), DIMENSION(0:,1:), INTENT(IN):: f1
REAL(KIND=r8), DIMENSION(1:,0:), INTENT(IN):: f2
REAL(KIND=r8), DIMENSION(1:SIZE(f2,1), &
                         1:SIZE(f1,2)):: udivresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(f2,1)
mx(2) = SIZE(f1,2)
!
udivresult(1:mx(1),1:mx(2)) = f1(1:mx(1)  ,1:mx(2)  ) &
                            - f1(0:mx(1)-1,1:mx(2)  ) &
                            + f2(1:mx(1)  ,1:mx(2)  ) &
                            - f2(1:mx(1)  ,0:mx(2)-1)
!
END FUNCTION UDiv2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULap2D(a) RESULT(ulapresult)
USE Global
IMPLICIT NONE
!
! 2D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1). The result is stored in
! ulapresult(1:mx(1),1:mx(2)).
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2, &
                         1:SIZE(a,2)-2):: ulapresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(a,1)-2
mx(2) = SIZE(a,2)-2
!
ulapresult(1:mx(1),1:mx(2)) =        a(2:mx(1)+1,1:mx(2)  ) &
                            +        a(0:mx(1)-1,1:mx(2)  ) &
                            +        a(1:mx(1)  ,2:mx(2)+1) &
                            +        a(1:mx(1)  ,0:mx(2)-1) &
                            - 4.0_r8*a(1:mx(1)  ,1:mx(2)  )
!
END FUNCTION ULap2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Restriction2D(a) RESULT(restrictionresult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(:,:,:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)/2, &
                         1:SIZE(a,2)/2, &
                         1:SIZE(a,3)):: restrictionresult
!
INTEGER:: nvars
INTEGER, DIMENSION(1:2):: cmx, mx
!
mx(1)  = SIZE(a,1)
mx(2)  = SIZE(a,2)
nvars  = SIZE(a,3)
cmx(1:2) = mx(1:2)/2
!
restrictionresult(1:cmx(1),1:cmx(2),1:nvars) &
  = 0.25_r8*(a(2:mx(1)  :2,2:mx(2)  :2,1:nvars) &
  +          a(1:mx(1)-1:2,2:mx(2)  :2,1:nvars) &
  +          a(2:mx(1)  :2,1:mx(2)-1:2,1:nvars) &
  +          a(1:mx(1)-1:2,1:mx(2)-1:2,1:nvars))
!
END FUNCTION Restriction2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Prolongation2D(a) RESULT(presult)
USE Global
IMPLICIT NONE
!
! This prolongation algorithm assumes no ghost layers are present.
!
REAL(KIND=r8), DIMENSION(:,:,:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)*2, &
                         1:SIZE(a,2)*2, &
                         1:SIZE(a,3)):: presult
!
INTEGER:: nvars
INTEGER, DIMENSION(1:2):: cmx, mx
!
cmx(1) = SIZE(a,1)
cmx(2) = SIZE(a,2)
nvars  = SIZE(a,3)
mx(1:2) = cmx(1:2)*2
!
presult(2:mx(1)  :2,2:mx(2)  :2,1:nvars) = a(1:cmx(1),1:cmx(2),1:nvars)
presult(1:mx(1)-1:2,2:mx(2)  :2,1:nvars) = a(1:cmx(1),1:cmx(2),1:nvars)
presult(2:mx(1)  :2,1:mx(2)-1:2,1:nvars) = a(1:cmx(1),1:cmx(2),1:nvars)
presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nvars) = a(1:cmx(1),1:cmx(2),1:nvars)
!
END FUNCTION Prolongation2D
END MODULE Utilities