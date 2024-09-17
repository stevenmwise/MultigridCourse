MODULE Utilities
USE Global
IMPLICIT NONE
!
CONTAINS
!
FUNCTION Norm(g) RESULT(normresult)
USE Global
IMPLICIT NONE
REAL(KIND=r8), DIMENSION(:,:), INTENT(IN) :: g
REAL(KIND=r8) :: normresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: tmp
!
mx(1) = SIZE(g,1); mx(2) = SIZE(g,2)
tmp = REAL(mx(1),KIND=r8)*REAL(mx(2),KIND=r8)
!
normresult = SQRT(SUM(g*g)/tmp)
!
END FUNCTION Norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Restriction(a) RESULT(restrictionresult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)/2, &
                         1:SIZE(a,2)/2):: restrictionresult
!
INTEGER, DIMENSION(1:2):: cmx, mx
!
mx(1)  = SIZE(a,1)
mx(2)  = SIZE(a,2)
cmx(1:2) = mx(1:2)/2
!
restrictionresult(1:cmx(1),1:cmx(2)) &
  = 0.25_r8*(a(2:mx(1)  :2,2:mx(2)  :2) &
  +          a(1:mx(1)-1:2,2:mx(2)  :2) &
  +          a(2:mx(1)  :2,1:mx(2)-1:2) &
  +          a(1:mx(1)-1:2,1:mx(2)-1:2))
!
END FUNCTION Restriction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Prolongation(a) RESULT(presult)
USE Global
IMPLICIT NONE
!
! This prolongation algorithm assumes no ghost layers are present in a.
!
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)*2, &
                         1:SIZE(a,2)*2):: presult
!
INTEGER, DIMENSION(1:2):: cmx, mx
!
cmx(1) = SIZE(a,1)
cmx(2) = SIZE(a,2)
mx(1:2) = cmx(1:2)*2
!
presult(2:mx(1)  :2,2:mx(2)  :2) = a(1:cmx(1),1:cmx(2))
presult(1:mx(1)-1:2,2:mx(2)  :2) = a(1:cmx(1),1:cmx(2))
presult(2:mx(1)  :2,1:mx(2)-1:2) = a(1:cmx(1),1:cmx(2))
presult(1:mx(1)-1:2,1:mx(2)-1:2) = a(1:cmx(1),1:cmx(2))
!
END FUNCTION Prolongation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION BiLinProlongation(a) RESULT(presult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0:(SIZE(a,2)-2)*2+1):: presult
!
INTEGER, DIMENSION(1:2):: cmx, mx
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0: SIZE(a,2)-2   +1):: b
!
cmx(1) = SIZE(a,1)-2
cmx(2) = SIZE(a,2)-2
mx(1:2) = cmx(1:2)*2
!
! Linear interpolation in the x-direction first:
            b(0: mx(1)  :2,0:cmx(2)+1  ) &
  =  3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ) &
  +         a(1:cmx(1)+1  ,0:cmx(2)+1  )
            b(1: mx(1)+1:2,0:cmx(2)+1  ) &
  =         a(0:cmx(1)    ,0:cmx(2)+1  ) &
  +  3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  )
!
! Linear interpolation in the y-direction second:
      presult(0: mx(1)+1  ,0: mx(2)  :2) &
  = (3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ) &
  +         b(0: mx(1)+1  ,1:cmx(2)+1  ))/16.0_r8
      presult(0: mx(1)+1  ,1: mx(2)+1:2) &
  = (       b(0: mx(1)+1  ,0:cmx(2)    ) &
  +  3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ))/16.0_r8
!
END FUNCTION BiLinProlongation
END MODULE Utilities