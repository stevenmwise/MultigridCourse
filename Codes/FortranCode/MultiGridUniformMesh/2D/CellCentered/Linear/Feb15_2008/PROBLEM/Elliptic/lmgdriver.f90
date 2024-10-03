PROGRAM NCGDriver
USE Global
USE LinearMultigridRoutines, ONLY: LinearMultigrid
IMPLICIT NONE
!
INTEGER:: i, ierror, iter, j, level
INTEGER, DIMENSION(1:2):: mx, cmx
REAL(KIND=r8):: pi = 3.1415926535897932_r8
REAL(KIND=r8):: c1, c2, tol, energyminimum, h, h1, h2, p1, p2, r, x, y
REAL(KIND=r8), DIMENSION(1:2):: xlower, xupper
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: f, res, u
!
NAMELIST/inputdata/tol, mx, minlevel, postsmooth, presmooth, maxits, omega, &
                   xlower, xupper
!
OPEN(UNIT=75,FILE='input.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file input.dat. Program stop.'
  STOP
END IF
READ(75,NML=inputdata)
CLOSE(75)
OPEN(UNIT=76,FILE='output.dat',STATUS='UNKNOWN',ACTION='WRITE', &
     FORM='FORMATTED', POSITION='APPEND')
WRITE(76,NML=inputdata)
CLOSE(76)
!
p1 = xupper(1)-xlower(1)
p2 = xupper(2)-xlower(2)
h1 = p1/mx(1)
h2 = p2/mx(2)
!
cmx = mx
DO level = 0, minlevel, -1
  IF(MODULO(cmx(1),2) == 1 .OR. MODULO(cmx(2),2) == 1) THEN
    IF(minlevel < level) THEN
      PRINT *, 'Refinement Error.  Stopping.'
      STOP
    END IF
  END IF
  cmx = cmx/2
END DO
!
IF(ABS(h1-h2) > 1.0E-14_r8) THEN
  PRINT *,'Stepsize error.  h1/=h2. Program stop.'
  STOP
ELSE
  h = h1
END IF
!
ALLOCATE(u(0:mx(1)+1,0:mx(2)+1),f(1:mx(1),1:mx(2)),res(1:mx(1),1:mx(2)))
!
u = 0.0_r8
!
c1 = (xlower(1)+xupper(1))/2.0_r8
c2 = (xlower(2)+xupper(2))/2.0_r8
!
DO i = 1, mx(1)
  x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
  DO j = 1, mx(2)
    y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
!
!    r = SQRT((x-c1)**2+(y-c1)**2)
!
!    f(i,j) = EXP(-1.0E+01_r8*r*r)
!    f(i,j) = 2.0_r8*SIN((x-xlower(1))*pi/p1)*SIN((y-xlower(2))*pi/p2)*pi*pi
    
    f(i,j) = EXP(COS(2.0_r8*pi*x/p1)*COS(2.0_r8*pi*y/p2))
!
!    CALL RANDOM_NUMBER(r)
!    r = 0.05_r8*(2.0_r8*r-1.0_r8)
!    u(i,j) = 0.0_r8+r
!    u(i,j) = SIN((x-xlower(1))*pi/p1)*SIN((y-xlower(2))*pi/p2)
    u(i,j) = 0.0_r8
!
  END DO
END DO
!
CALL LinearMultigrid(u,f,res,h,tol,iter,energyminimum)
!
CALL PrintOUT
!
DEALLOCATE(u,f,res)
!
CONTAINS
!
SUBROUTINE PrintOut
IMPLICIT NONE
!
CHARACTER(LEN=18):: file_name
!
file_name = './OUT/solution.dat'
OPEN(UNIT=9, FILE=file_name, STATUS='REPLACE', ACTION='READWRITE')
WRITE(9,'(F25.12)') 0.0
WRITE(9,'(I8)') 3
WRITE(9,'(F25.12)') h
WRITE(9,'(F25.12)') h
WRITE(9,'(F25.12)') xlower(1)
WRITE(9,'(F25.12)') xlower(2)
WRITE(9,'(F25.12)') xupper(1)
WRITE(9,'(F25.12)') xupper(2)
WRITE(9,'(I8)') mx(1)
WRITE(9,'(I8)') mx(2)
DO j = 1, mx(2)
  y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
  DO i = 1, mx(1)
    x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
    WRITE(9,'(3(f25.12,1x),f25.12)') x, y, u(i,j), f(i,j), res(i,j)
  END DO
END DO
CLOSE(9)
!
END SUBROUTINE PrintOut
!
END PROGRAM NCGDriver
