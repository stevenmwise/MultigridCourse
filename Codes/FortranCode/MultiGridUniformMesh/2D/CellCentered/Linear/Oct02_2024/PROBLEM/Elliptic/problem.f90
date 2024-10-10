MODULE Problem
USE Global
!
CONTAINS
!
FUNCTION Residual(u,f,hf) RESULT(residualResult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: hf
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2, &
                         1:SIZE(u,2)-2):: residualResult
!
INTEGER, DIMENSION(1:2):: nf
!
CALL BoundaryConditions(u)
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
!
residualResult(1:nf(1),1:nf(2)) =            f(1:nf(1)  ,1:nf(2)  ) &
                                - FDOperator(u(0:nf(1)+1,0:nf(2)+1),hf)
!
END FUNCTION Residual
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION FDOperator(u,hf) RESULT(operatorResult)
USE Global
IMPLICIT NONE
!
! Helmholtz Operator:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), INTENT(IN):: hf
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2, &
                         1:SIZE(u,2)-2):: operatorResult
!
INTEGER, DIMENSION(1:2):: nf
REAL(KIND=r8):: hf2
!
CALL BoundaryConditions(u)
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
hf2 = hf*hf
!
operatorResult(1:nf(1),1:nf(2)) = -(       u(2:nf(1)+1,1:nf(2)  ) &
                                  +        u(0:nf(1)-1,1:nf(2)  ) &
                                  +        u(1:nf(1)  ,2:nf(2)+1) &
                                  +        u(1:nf(1)  ,0:nf(2)-1) &
                                  - 4.0_r8*u(1:nf(1)  ,1:nf(2)  ))/hf2 &
                                  +        u(1:nf(1)  ,1:nf(2)  )
!
END FUNCTION FDOperator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Smooth(u,f,hf,passes) RESULT(smoothResult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: hf
INTEGER, INTENT(IN):: passes
!
INTEGER:: i, its, j
INTEGER, DIMENSION(1:2):: nf
REAL(KIND=r8):: hf2, omega2, tmp
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-1,0:SIZE(u,2)-1):: smoothResult
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
hf2 = hf*hf
omega2 = 1.0_r8-omega
!
CALL BoundaryConditions(u)
!
DO its = 1, passes
!
DO j = 1, nf(2)
  DO i = 1, nf(1)
!
    tmp = (hf2*f(i,j)+u(i+1,j  )+u(i-1,j  ) &
        +             u(i  ,j+1)+u(i  ,j-1))/(4.0_r8+hf2)
    u(i,j) = omega*tmp+omega2*u(i,j)
!
  END DO
END DO
!
CALL BoundaryConditions(u)
!
END DO
!
smoothResult = u
!
END Function Smooth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE BoundaryConditions(u)
USE Global
IMPLICIT NONE
!
! Homogeneous Neumann boundary conditions:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
!
INTEGER, DIMENSION(1:2):: nf
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
!
u(1:nf(1),      0) = u(1:nf(1),    1)
u(1:nf(1),nf(2)+1) = u(1:nf(1),nf(2))
!
u(      0,0:nf(2)+1) = u(    1,0:nf(2)+1)
u(nf(1)+1,0:nf(2)+1) = u(nf(1),0:nf(2)+1)
!
END SUBROUTINE BoundaryConditions
END MODULE Problem
