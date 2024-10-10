PROGRAM LMGDriver
USE Global
USE LinearMultigridRoutines, ONLY: LinearMultigrid
USE Problem, ONLY: FDOperator, BoundaryConditions
IMPLICIT NONE
!
INTEGER:: i, iError, iter, j, level
INTEGER, DIMENSION(1:2):: nf, nc
REAL(KIND=r8):: pi = 3.1415926535897932_r8
REAL(KIND=r8):: c1, c2, tol, h, h1, h2, p1, p2, r, x, y
REAL(KIND=r8), DIMENSION(1:2):: xLower, xUpper
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: f, res, u
!
NAMELIST/inputData/tol, nf, maxLevel, postSmooth, preSmooth, pCycle, maxItrs, &
                   omega, xLower, xUpper
!
OPEN(UNIT=75,FILE='input.dat',STATUS='OLD',ACTION='READ',IOSTAT=iError)
IF(iError/=0) THEN
  PRINT *,'Error opening input file input.dat. Stopping.'
  STOP
END IF
READ(75,NML=inputData)
CLOSE(75)
OPEN(UNIT=76,FILE='output.dat',STATUS='UNKNOWN',ACTION='WRITE', &
     FORM='FORMATTED', POSITION='APPEND')
WRITE(76,NML=inputData)
CLOSE(76)
!
p1 = xUpper(1)-xLower(1)
p2 = xUpper(2)-xLower(2)
h1 = p1/nf(1)
h2 = p2/nf(2)
!
IF(maxLevel <= 0) THEN
  PRINT *, 'maxLevel must be a positive integer. Stopping.'
  STOP
END IF
!
nc = nf
DO level = maxLevel, 0, -1
  IF(MODULO(nc(1),2) == 1 .OR. MODULO(nc(2),2) == 1) THEN
    IF(0 < level) THEN
      PRINT *, 'Refinement Error. Stopping.'
      STOP
    END IF
  END IF
  nc = nc/2
END DO
!
IF(ABS(h1-h2) > 1.0E-14_r8) THEN
  PRINT *,'Stepsize error.  h1/=h2. Stopping.'
  STOP
ELSE
  h = h1
END IF
!
ALLOCATE(u(0:nf(1)+1,0:nf(2)+1),f(1:nf(1),1:nf(2)))
!
c1 = (xLower(1)+xUpper(1))/2.0_r8
c2 = (xLower(2)+xUpper(2))/2.0_r8
!
DO i = 1, nf(1)
  x = (REAL(i,KIND=r8)-0.5_r8)*h+xLower(1)
  DO j = 1, nf(2)
    y = (REAL(j,KIND=r8)-0.5_r8)*h+xLower(2)
!   
    u(i,j) = EXP(COS(2.0_r8*pi*x/p1)*COS(2.0_r8*pi*y/p2))
!
  END DO
END DO
!
CALL BoundaryConditions(u)
f = FDOperator(u,h)
u = 0.0_r8
!
u = LinearMultigrid(u,f,h,tol,iter)
!
CALL PrintOUT
!
DEALLOCATE(u,f)
!
CONTAINS
!
SUBROUTINE PrintOut
IMPLICIT NONE
!
CHARACTER(LEN=18):: fileName
!
fileName = './OUT/solution.dat'
OPEN(UNIT=9, FILE=fileName, STATUS='REPLACE', ACTION='READWRITE')
WRITE(9,'(F25.12)') 0.0
WRITE(9,'(I8)') 2
WRITE(9,'(F25.12)') h
WRITE(9,'(F25.12)') h
WRITE(9,'(F25.12)') xLower(1)
WRITE(9,'(F25.12)') xLower(2)
WRITE(9,'(F25.12)') xUpper(1)
WRITE(9,'(F25.12)') xUpper(2)
WRITE(9,'(I8)') nf(1)
WRITE(9,'(I8)') nf(2)
DO j = 1, nf(2)
  y = (REAL(j,KIND=r8)-0.5_r8)*h+xLower(2)
  DO i = 1, nf(1)
    x = (REAL(i,KIND=r8)-0.5_r8)*h+xLower(1)
    WRITE(9,'(3(f25.12,1x),f25.12)') x, y, u(i,j), f(i,j)
  END DO
END DO
CLOSE(9)
!
END SUBROUTINE PrintOut
!
END PROGRAM LMGDriver
