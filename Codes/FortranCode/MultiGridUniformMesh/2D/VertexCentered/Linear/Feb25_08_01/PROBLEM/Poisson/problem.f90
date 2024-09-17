MODULE Problem
USE Global
!
CONTAINS
!
FUNCTION Energy(u,f,h) RESULT(energyresult)
USE Global
IMPLICIT NONE
!
! Dirichlet Energy:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8):: energyresult
!
energyresult = 0.0_r8
!
END FUNCTION Energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION EnergyGradient(u,f,h) RESULT(energygradientresult)
USE Global
IMPLICIT NONE
!
! Dirichlet Energy Gradient:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2, &
                         1:SIZE(u,2)-2):: energygradientresult
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
!
CALL BoundaryConditions(u)
!
mx(1) = SIZE(u,1)-1; mx(2) = SIZE(u,2)-1
!
energygradientresult(1:mx(1)-1,1:mx(2)-1) = Operator(u(0:mx(1)  ,0:mx(2)  ),h) &
                                          -          f(1:mx(1)-1,1:mx(2)-1)
!
END FUNCTION EnergyGradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Operator(u,h) RESULT(operatorresult)
USE Global
IMPLICIT NONE
!
! Laplacian Operator:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2, &
                         1:SIZE(u,2)-2):: operatorresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: h2
!
CALL BoundaryConditions(u)
!
mx(1) = SIZE(u,1)-1; mx(2) = SIZE(u,2)-1
!
h2 = h*h
!
operatorresult(1:mx(1)-1,1:mx(2)-1) = -(     u(2:mx(1)  ,1:mx(2)-1) &
                                    +        u(0:mx(1)-2,1:mx(2)-1) &
                                    +        u(1:mx(1)-1,2:mx(2)  ) &
                                    +        u(1:mx(1)-1,0:mx(2)-2) &
                                    - 4.0_r8*u(1:mx(1)-1,1:mx(2)-1))/h2
!
END FUNCTION Operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Relaxation(u,f,h,passes)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: h
INTEGER, INTENT(IN):: passes
!
INTEGER:: i, its, j
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: h2
!
mx(1) = SIZE(u,1)-1; mx(2) = SIZE(u,2)-1
h2 = h*h
!
CALL BoundaryConditions(u)
!
DO its = 1, passes
!
  DO j = 1, mx(2)-1
    DO i = 1, mx(1)-1
!
      u(i,j) = (h2*f(i,j)+u(i+1,j  )+u(i-1,j  ) &
             +            u(i  ,j+1)+u(i  ,j-1))/4.0_r8
!
    END DO
  END DO
!
END DO
!
END SUBROUTINE Relaxation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE BoundaryConditions(u)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(u,1)-1; mx(2) = SIZE(u,2)-1
!
u(0:mx(1),    0) = 0.0_r8
u(0:mx(1),mx(2)) = 0.0_r8
!
u(    0,1:mx(2)-1) = 0.0_r8
u(mx(1),1:mx(2)-1) = 0.0_r8
!
END SUBROUTINE BoundaryConditions
END MODULE Problem