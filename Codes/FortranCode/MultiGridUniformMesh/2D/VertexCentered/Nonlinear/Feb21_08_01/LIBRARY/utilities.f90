MODULE Utilities
USE Global
IMPLICIT NONE
!
CONTAINS
!
FUNCTION Norm(g) RESULT(normresult)
USE Global
IMPLICIT NONE
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN) :: g
REAL(KIND=r8) :: normresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: tmp
REAL(KIND=r8), DIMENSION(0:SIZE(g,1)-1,0:SIZE(g,2)-1):: g2
!
mx(1) = SIZE(g,1)-1
mx(2) = SIZE(g,2)-1
tmp = REAL(mx(1)*mx(2),KIND=r8)
!
g2 = g*g
!
g2(    0,    0) = 0.25_r8*g2(    0,    0)
g2(mx(1),    0) = 0.25_r8*g2(mx(1),    0)
g2(    0,mx(2)) = 0.25_r8*g2(    0,mx(2))
g2(mx(1),mx(2)) = 0.25_r8*g2(mx(1),mx(2))
!
g2(        0,1:mx(2)-1) = 0.5_r8*g2(        0,1:mx(2)-1)
g2(  mx(1)  ,1:mx(2)-1) = 0.5_r8*g2(  mx(1)  ,1:mx(2)-1)
g2(1:mx(1)-1,        0) = 0.5_r8*g2(1:mx(1)-1,        0)
g2(1:mx(1)-1,  mx(2)  ) = 0.5_r8*g2(1:mx(1)-1,  mx(2)  )
!
normresult = SQRT(SUM(g2)/tmp)
!
END FUNCTION Norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Restriction2(a) RESULT(restrictionresult)
USE Global
IMPLICIT NONE
!
! Injection:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-1)/2, &
                         0:(SIZE(a,2)-1)/2):: restrictionresult
!
INTEGER, DIMENSION(1:2):: cmx, mx
!
mx(1)  = SIZE(a,1)-1
mx(2)  = SIZE(a,2)-1
cmx(1:2) = mx(1:2)/2
!
! Boundary points first:
restrictionresult(0:cmx(1),       0  ) = a(0:mx(1):2,      0    )
restrictionresult(0:cmx(1),  cmx(2)  ) = a(0:mx(1):2,  mx(2)    )
restrictionresult(       0,1:cmx(2)-1) = a(      0  ,2:mx(2)-2:2)
restrictionresult(  cmx(1),1:cmx(2)-1) = a(  mx(1)  ,2:mx(2)-2:2)
!
! Interior points second:
restrictionresult(1:cmx(1)-1,1:cmx(2)-1) = (4.0_r8*a(2:mx(1)-2:2,2:mx(2)-2:2) &
                                         +  2.0_r8*a(1:mx(1)-3:2,2:mx(2)-2:2) &
                                         +  2.0_r8*a(3:mx(1)-1:2,2:mx(2)-2:2) &
                                         +  2.0_r8*a(2:mx(1)-2:2,1:mx(2)-3:2) &
                                         +  2.0_r8*a(2:mx(1)-2:2,3:mx(2)-1:2) &
                                         +         a(1:mx(1)-3:2,1:mx(2)-3:2) &
                                         +         a(3:mx(1)-1:2,3:mx(2)-1:2) &
                                         +         a(3:mx(1)-1:2,1:mx(2)-3:2) &
                                         +         a(1:mx(1)-3:2,3:mx(2)-1:2))/16.0_r8
!
END FUNCTION Restriction2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Restriction1(a) RESULT(restrictionresult)
USE Global
IMPLICIT NONE
!
! Injection:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-1)/2, &
                         0:(SIZE(a,2)-1)/2):: restrictionresult
!
INTEGER, DIMENSION(1:2):: cmx, mx
!
mx(1)  = SIZE(a,1)-1
mx(2)  = SIZE(a,2)-1
cmx(1:2) = mx(1:2)/2
!
! Boundary points first:
restrictionresult(0:cmx(1),       0  ) = a(0:mx(1):2,      0    )
restrictionresult(0:cmx(1),  cmx(2)  ) = a(0:mx(1):2,  mx(2)    )
restrictionresult(       0,1:cmx(2)-1) = a(      0  ,2:mx(2)-2:2)
restrictionresult(  cmx(1),1:cmx(2)-1) = a(  mx(1)  ,2:mx(2)-2:2)
!
! Interior points second:
restrictionresult(1:cmx(1)-1,1:cmx(2)-1) = (4.0_r8*a(2:mx(1)-2:2,2:mx(2)-2:2) &
                                         +         a(1:mx(1)-3:2,2:mx(2)-2:2) &
                                         +         a(3:mx(1)-1:2,2:mx(2)-2:2) &
                                         +         a(2:mx(1)-2:2,1:mx(2)-3:2) &
                                         +         a(2:mx(1)-2:2,3:mx(2)-1:2))/8.0_r8
!
END FUNCTION Restriction1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Restriction(a) RESULT(restrictionresult)
USE Global
IMPLICIT NONE
!
! Injection:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-1)/2, &
                         0:(SIZE(a,2)-1)/2):: restrictionresult
!
INTEGER, DIMENSION(1:2):: cmx, mx
!
mx(1)  = SIZE(a,1)-1
mx(2)  = SIZE(a,2)-1
cmx(1:2) = mx(1:2)/2
!
restrictionresult(0:cmx(1),0:cmx(2)) = a(0:mx(1):2,0:mx(2):2)
!
END FUNCTION Restriction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Prolongation(a) RESULT(presult)
USE Global
IMPLICIT NONE
!
! Interpolation:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-1)*2, &
                         0:(SIZE(a,2)-1)*2):: presult
!
INTEGER, DIMENSION(1:2):: cmx, mx
!
cmx(1) = SIZE(a,1)-1
cmx(2) = SIZE(a,2)-1
mx(1:2) = cmx(1:2)*2
!
presult(0:mx(1):2,0:mx(2):2) = a(0:cmx(1),0:cmx(2))
!
presult(1:mx(1)-1:2,0:mx(2):2) &
  = 0.5_r8*(presult(0:mx(1)-2:2,0:mx(2):2)+presult(2:mx(1):2,0:mx(2):2))
!
presult(0:mx(1),1:mx(2)-1:2) &
  = 0.5_r8*(presult(0:mx(1),0:mx(2)-2:2)+presult(0:mx(1),2:mx(2):2))
!
END FUNCTION Prolongation
END MODULE Utilities