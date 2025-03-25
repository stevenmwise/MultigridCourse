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
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2):: cocc
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,1:SIZE(u,2)-2):: f1
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,0:SIZE(u,2)-2):: f2
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2,1:2):: xcc
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,1:SIZE(u,2)-2,1:2):: xew
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,0:SIZE(u,2)-2,1:2):: xns
!
CALL BoundaryConditions(u)
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
hf2 = hf*hf
!
xcc(1:nf(1),1:nf(2),1:2) = GetCCCoords(hf,nf)
xew(0:nf(1),1:nf(2),1:2) = GetEWCoords(hf,nf)
xns(1:nf(1),0:nf(2),1:2) = GetNSCoords(hf,nf)
!
f1(0:nf(1),1:nf(2)) = Df(xew(0:nf(1)  ,1:nf(2),1:2)) &
                    * (    u(1:nf(1)+1,1:nf(2)    )-u(0:nf(1),1:nf(2)  ))
f2(1:nf(1),0:nf(2)) = Df(xns(1:nf(1)  ,0:nf(2),1:2)) &
                    * (    u(1:nf(1)  ,1:nf(2)+1  )-u(1:nf(1)  ,0:nf(2)))
!
operatorResult(1:nf(1),1:nf(2)) &
  = -(f1(1:nf(1),1:nf(2))-f1(0:nf(1)-1,1:nf(2)  ))/hf2 &
  -  (f2(1:nf(1),1:nf(2))-f2(1:nf(1)  ,0:nf(2)-1))/hf2 &
  + Co(xcc(1:nf(1),1:nf(2),1:2))*u(1:nf(1),1:nf(2))
!
END FUNCTION FDOperator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Df(x) RESULT(dfResult)
USE Global
IMPLICIT NONE
!
! Diffusivity:
!
REAL(KIND=r8), DIMENSION(:,:,1:), INTENT(IN):: x
!
INTEGER, DIMENSION(1:2):: nf
REAL(KIND=r8), DIMENSION(1:SIZE(x,1),1:SIZE(x,2)):: dfResult
!
nf(1) = SIZE(x,1)
nf(2) = SIZE(x,2)

dfResult(1:nf(1),1:nf(2)) = x(1:nf(1),1:nf(2),1)*x(1:nf(1),1:nf(2),1) &
                          + x(1:nf(1),1:nf(2),2)*x(1:nf(1),1:nf(2),2)+1.0_r8

! dfResult(1:nf(1),1:nf(2)) = 1.0_r8
!
END FUNCTION Df
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Co(x) RESULT(coResult)
USE Global
IMPLICIT NONE
!
! Diffusivity:
!
INTEGER, DIMENSION(1:2):: nf
REAL(KIND=r8), DIMENSION(:,:,1:), INTENT(IN):: x
REAL(KIND=r8), DIMENSION(1:SIZE(x,1),1:SIZE(x,2)):: coResult
!
nf(1) = SIZE(x,1)
nf(2) = SIZE(x,2)
!
coResult(1:nf(1),1:nf(2)) = x(1:nf(1),1:nf(2),1)*x(1:nf(1),1:nf(2),1) &
                          + x(1:nf(1),1:nf(2),2)*x(1:nf(1),1:nf(2),2)+1.0_r8
!
! coResult(1:nf(1),1:nf(2)) = 1.0_r8
!
END FUNCTION Co
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION GetCCCoords(hf,nf) RESULT(coordResult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: hf
INTEGER, DIMENSION(1:2), INTENT(IN):: nf
!
INTEGER:: i, j
REAL(KIND=r8):: x, y
REAL(KIND=r8), DIMENSION(1:nf(1),1:nf(2),1:2):: coordResult
!
DO j = 1, nf(2)
  y = (j-0.5_r8)*hf+xLower(2)
  DO i = 1, nf(1)
    x = (i-0.5_r8)*hf+xLower(1)
!
    coordResult(i,j,1) = x
    coordResult(i,j,2) = y 
!
  END DO
END DO
!
END Function GetCCCoords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION GetEWCoords(hf,nf) RESULT(coordResult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: hf
INTEGER, DIMENSION(1:2), INTENT(IN):: nf
!
INTEGER:: i, j
REAL(KIND=r8):: x, y
REAL(KIND=r8), DIMENSION(0:nf(1),1:nf(2),1:2):: coordResult
!
DO j = 1, nf(2)
  y = (j-0.5_r8)*hf+xLower(2)
  DO i = 0, nf(1)
    x = i*hf+xLower(1)
!
    coordResult(i,j,1) = x
    coordResult(i,j,2) = y 
!
  END DO
END DO
!
END Function GetEWCoords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION GetNSCoords(hf,nf) RESULT(coordResult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: hf
INTEGER, DIMENSION(1:2), INTENT(IN):: nf
!
INTEGER:: i, j
REAL(KIND=r8):: x, y
REAL(KIND=r8), DIMENSION(1:nf(1),0:nf(2),1:2):: coordResult
!
DO j = 0, nf(2)
  y = j*hf+xLower(2)
  DO i = 1, nf(1)
    x = (i-0.5_r8)*hf+xLower(1)
!
    coordResult(i,j,1) = x
    coordResult(i,j,2) = y 
!
  END DO
END DO
!
END Function GetNSCoords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Smooth(u,f,hf,passes,direction) RESULT(smoothResult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: hf
INTEGER, INTENT(IN):: passes
CHARACTER(LEN=3), INTENT(IN):: direction
!
INTEGER:: i, ifirst, ilast, istride, its, j, jfirst, jlast, jstride
INTEGER, DIMENSION(1:2):: nf
REAL(KIND=r8):: hf2, omegaPrime, tmp
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-1,0:SIZE(u,2)-1):: smoothResult
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2):: cocc
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,1:SIZE(u,2)-2):: f1
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,0:SIZE(u,2)-2):: f2
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2,1:2):: xcc
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,1:SIZE(u,2)-2,1:2):: xew
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,0:SIZE(u,2)-2,1:2):: xns
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
hf2 = hf*hf
omegaPrime = 1.0_r8-omega
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
hf2 = hf*hf
!
xcc(1:nf(1),1:nf(2),1:2) = GetCCCoords(hf,nf)
xew(0:nf(1),1:nf(2),1:2) = GetEWCoords(hf,nf)
xns(1:nf(1),0:nf(2),1:2) = GetNSCoords(hf,nf)
!
f1(0:nf(1),1:nf(2)) = Df(xew(0:nf(1)  ,1:nf(2),1:2))
f2(1:nf(1),0:nf(2)) = Df(xns(1:nf(1)  ,0:nf(2),1:2))
cocc(1:nf(1),1:nf(2)) = Co(xcc(1:nf(1),1:nf(2),1:2))
!
SELECT CASE (direction)
  CASE ('fwd')
    jFirst = 1
    jLast = nf(2)
    jStride = 1
    iFirst = 1
    iLast = nf(1) 
    iStride = 1
  CASE DEFAULT
    jFirst = nf(2)
    jLast = 1
    jStride = -1
    iFirst = nf(1)
    iLast = 1
    iStride = -1
END SELECT
!
CALL BoundaryConditions(u)
!
DO its = 1, passes
!
DO j = jFirst, jLast, jStride
  DO i = iFirst, iLast, iStride
!
    tmp = (hf2*f(i,j)+u(i+1,j  )*f1(i,j)+u(i-1,j  )*f1(i-1,j) &
        +             u(i  ,j+1)*f2(i,j)+u(i  ,j-1)*f2(i,j-1)) &
        / (f1(i,j)+f1(i-1,j)+f2(i,j)+f2(i,j-1)+hf2*cocc(i,j))
    u(i,j) = omega*tmp+omegaPrime*u(i,j)
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
