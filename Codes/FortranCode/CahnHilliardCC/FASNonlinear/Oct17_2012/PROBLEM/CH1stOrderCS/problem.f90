MODULE Problem
USE Global
!
! Cahn-Hilliard Equation: 1st Order Scheme
!
! Semi-implicit, convexity splitting formulation.
!
CONTAINS
!
SUBROUTINE SetProblem
USE Global
USE ProblemDef
IMPLICIT NONE
!
INTEGER:: ierror
NAMELIST/problemdata/ eps
!
OPEN(UNIT=75,FILE='problemdata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file problemdata.dat. Program stop.'
  STOP
END IF
READ(75,NML=problemdata)
CLOSE(75)
OPEN(UNIT=76,FILE='output.dat',STATUS='OLD',ACTION='WRITE',FORM='FORMATTED', &
     POSITION='APPEND')
WRITE(76,NML=problemdata)
CLOSE(76)
!
END SUBROUTINE SetProblem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Initialize(u,aux,h,xlower)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(OUT):: aux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
!
INTEGER:: i, j, myseed, thesize
INTEGER, DIMENSION(1:2):: mx
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
REAL(KIND=r8):: r
REAL(KIND=r8), DIMENSION(1:2):: c, p, x
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
myseed = 1783221290
!
CALL RANDOM_SEED(SIZE=thesize)
ALLOCATE(seed(thesize))
DO j = 1, thesize
  seed(j) = ABS(myseed)+(j-1)
END DO
CALL RANDOM_SEED(PUT=seed)
DEALLOCATE(seed)
!
p(1:2) = mx(1:2)*h
c(1:2) = (xlower(1:2)+(xlower(1:2)+p(1:2)))/2.0_r8
!
DO i = 1, mx(1)
  x(1) = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
  DO j = 1, mx(2)
    x(2) = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
!
    CALL RANDOM_NUMBER(r)
    r = 0.005_r8*(2.0_r8*r-1.0_r8)
    u(i,j,1) = -0.0_r8+r
!    
!    u(i,j,1) = 0.5_r8*SIN(2.0_r8*pi*x(1))
!
!    r = SQRT(0.5_r8*(x(1)-0.5_r8)**2+(x(2)-0.5_r8)**2)
!    
!    u(i,j,1) = -tanh((r-0.25_r8)/(SQRT(2.0_r8)*eps))
!
    u(i,j,2) =  0.0_r8
!
  END DO
END DO
!
aux = 0.0_r8
!
CALL BoundaryConditions(u)
!
END SUBROUTINE Initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetForcingTerm(f,h,xlower)
USE Global
USE ProblemDef
IMPLICIT NONE
!
! Possibly time dependent forcing term:
!
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(OUT):: f
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8), DIMENSION(1:2):: c, p, r, x
!
mx(1) = SIZE(f,1); mx(2) = SIZE(f,2)
!
p(1:2) = mx(1:2)*h
c(1:2) = (xlower(1:2)+(xlower(1:2)+p(1:2)))/2.0_r8
!
f = 0.0_r8
!
END SUBROUTINE SetForcingTerm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Energy(u,f,aux,h) RESULT(energyresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
! Possibly time-dependent energy.
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8):: energyresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: h2
!
h2 = h*h
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
CALL BoundaryConditions(u(0:mx(1)+1,0:mx(2)+1,1:numvars))
!
energyresult = 0.0_r8
!
END FUNCTION Energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Operator(u,uo,uoo,aux,h) RESULT(operatorresult)
USE Global
USE Utilities, ONLY: ULap2D
USE ProblemDef
IMPLICIT NONE
!
! Operator:
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2,1:numvars):: operatorresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: h2, tmp1, tmp2
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
h2 = h*h
tmp1 = dt/h2
tmp2 = eps/h2
!
operatorresult(1:mx(1),1:mx(2),1) &
  = u(1:mx(1),1:mx(2),1) &
  - tmp1*ULap2D(u(0:mx(1)+1,0:mx(2)+1,2))
!
operatorresult(1:mx(1),1:mx(2),2) &
  =    u(1:mx(1),1:mx(2),2) &
  - Cube(u(1:mx(1),1:mx(2),1))/eps &
  + tmp2*ULap2D(u(0:mx(1)+1,0:mx(2)+1,1))
!
END FUNCTION Operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetSource(u,uo,uoo,f,aux,s,h)
USE Global
USE Utilities, ONLY: ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: u
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN OUT):: s
REAL(KIND=r8), INTENT(IN):: h
!
INTEGER, DIMENSION(1:2):: mx
!
REAL(KIND=r8):: h2
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
s(1:mx(1),1:mx(2),1) =  uo(1:mx(1),1:mx(2),1)
!
s(1:mx(1),1:mx(2),2) = -uo(1:mx(1),1:mx(2),1)/eps
!
END SUBROUTINE SetSource
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Relaxation(u,uo,uoo,s,aux,h,passes,level)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: s
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), INTENT(IN):: h
INTEGER, INTENT(IN):: passes
INTEGER, INTENT(IN):: level
!
INTEGER:: i, ipass, j, redblack
INTEGER, DIMENSION(1:2):: mx
!
REAL(KIND=r8):: det, h2, tmp1, tmp2, tmp3, tmp4
REAL(KIND=r8), DIMENSION(1:2,1:2):: a
REAL(KIND=r8), DIMENSION(1:2):: b
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
h2 = h*h
tmp1 = dt/h2
tmp2 = eps/h2
tmp4 = 4.0_r8*tmp2
!
a(1,1) = 1.0_r8
a(1,2) = 4.0_r8*tmp1
a(2,2) = 1.0_r8
!
CALL BoundaryConditions(u(0:mx(1)+1,0:mx(2)+1,1:numvars))
!
DO ipass = 1, passes
!
DO redblack = 1, 2
!
DO j = 1, mx(2)
  DO i = 1+MODULO(j+redblack,2), mx(1), 2
!
    tmp3 = Sqr(u(i,j,1))/eps
!
    a(2,1) = -tmp3-tmp4
!
    b(1) = s(i,j,1)+tmp1*(u(i-1,j,2)+u(i+1,j,2) &
         +                u(i,j-1,2)+u(i,j+1,2))
!
    b(2) = s(i,j,2)-tmp2*(u(i+1,j,1)+u(i-1,j,1) &
         +                u(i,j+1,1)+u(i,j-1,1))
!
! Solve using Cramer's rule:
    det = 1.0_r8-a(1,2)*a(2,1)
!
    u(i,j,1) = (b(1)-a(1,2)*b(2))/det
    u(i,j,2) = (b(2)-a(2,1)*b(1))/det
!
  END DO
END DO
!
CALL BoundaryConditions(u(0:mx(1)+1,0:mx(2)+1,1:numvars))
!
END DO
!
END DO
!
END SUBROUTINE Relaxation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE BoundaryConditions(u)
USE Global
USE ProblemDef
IMPLICIT NONE
!
! Boundary conditions:
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: u
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
! Periodic boundary conditions:
!
u(      0,1:mx(2),1:numvars) = u(mx(1),1:mx(2),1:numvars)
u(mx(1)+1,1:mx(2),1:numvars) = u(    1,1:mx(2),1:numvars)
!
u(0:mx(1)+1,      0,1:numvars) = u(0:mx(1)+1,mx(2),1:numvars)
u(0:mx(1)+1,mx(2)+1,1:numvars) = u(0:mx(1)+1,    1,1:numvars)
!
END SUBROUTINE BoundaryConditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE BoundaryConditionsAux(aux)
USE Global
USE ProblemDef
IMPLICIT NONE
!
! Boundary conditions for the auxiliary variables:
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: aux
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(aux,1)-2; mx(2) = SIZE(aux,2)-2
!
! Homogeneous Neumann boundary conditions:
!
aux(      0,1:mx(2),1:auxvars) = aux(    1,1:mx(2),1:auxvars)
aux(mx(1)+1,1:mx(2),1:auxvars) = aux(mx(1),1:mx(2),1:auxvars)
!
aux(0:mx(1)+1,      0,1:auxvars) = aux(0:mx(1)+1,    1,1:auxvars)
aux(0:mx(1)+1,mx(2)+1,1:auxvars) = aux(0:mx(1)+1,mx(2),1:auxvars)
!
END SUBROUTINE BoundaryConditionsAux
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Optional user supplied routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULapMob2D(p,a) RESULT(ulapmobresult)
USE Global
USE ProblemDef
USE Utilities, ONLY: UDiv2D
IMPLICIT NONE
!
! Level independent, 2D UNDIVIDED laplacian operator for non-constant mobility.
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: p
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,1:SIZE(a,2)-2):: ulapmobresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(a,1)-2,1:SIZE(a,2)-2):: f1
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,0:SIZE(a,2)-2):: f2
!
mx(1) = SIZE(a,1)-2; mx(2) = SIZE(a,2)-2
!
! Calculate the UNDIVIDED 2D flux function:
f1(0:mx(1),1:mx(2)) &
  = Mob(0.5_r8*(p(1:mx(1)+1,1:mx(2))+p(0:mx(1),1:mx(2)))) &
  *            (a(1:mx(1)+1,1:mx(2))-a(0:mx(1),1:mx(2)))
!
f2(1:mx(1),0:mx(2)) &
  = Mob(0.5_r8*(p(1:mx(1),1:mx(2)+1)+p(1:mx(1),0:mx(2)))) &
  *            (a(1:mx(1),1:mx(2)+1)-a(1:mx(1),0:mx(2)))
!
! Calculate the UNDIVIDED divergence of the flux:
ulapmobresult(1:mx(1),1:mx(2)) = UDiv2D(f1(0:mx(1),1:mx(2)), &
                                        f2(1:mx(1),0:mx(2)))
!
END FUNCTION ULapMob2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Mob(p) RESULT(mobresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: mobresult
!
mobresult = 1.0_r8
!
END FUNCTION Mob
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Sqr(p) RESULT(sqrresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: sqrresult
!
sqrresult = p*p
!
END FUNCTION Sqr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Cube(p) RESULT(cuberesult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: cuberesult
!
cuberesult = p*p*p
!
END FUNCTION Cube
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Quad(p) RESULT(quadresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: quadresult
!
quadresult = p*p*p*p
!
END FUNCTION Quad
!
END MODULE Problem
