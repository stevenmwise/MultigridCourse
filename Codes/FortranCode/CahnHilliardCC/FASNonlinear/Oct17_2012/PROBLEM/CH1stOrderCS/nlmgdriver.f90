PROGRAM NLMGDriver
USE Global
USE NonLinearMultigridRoutines, ONLY: NonLinearMultigrid
USE Problem, ONLY: SetProblem, Initialize, Energy, SetForcingTerm
IMPLICIT NONE
!
LOGICAL:: restart
INTEGER:: i, ierror, iter, j, k, level, maxtimeiterations, numprintouts, &
          printiterations, restartnum, startprintout
INTEGER, DIMENSION(1:2):: mx, cmx
REAL(KIND=r8):: tol, energyvalue, h
REAL(KIND=r8), DIMENSION(1:2):: hh, p, xlower, xupper
REAL(KIND=r8), DIMENSION(:,:,:), ALLOCATABLE:: aux, f, res, u, uo, uoo
!
NAMELIST/inputdata/tol, mx, numvars, numauxvars, maxtimeiterations, dt, &
                   numprintouts, minlevel, postsmooth, presmooth, getenergy, &
                   maxits, omega, xlower, xupper, restart, restartnum
!
OPEN(UNIT=75, FILE='input.dat', STATUS='OLD', ACTION='READ', IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file input.dat. Program stop.'
  STOP
END IF
READ(75,NML=inputdata)
CLOSE(75)
OPEN(UNIT=76, FILE='output.dat', STATUS='UNKNOWN', ACTION='WRITE', &
     FORM='FORMATTED', POSITION='APPEND')
WRITE(76,NML=inputdata)
CLOSE(76)
!
p(1:2) = xupper(1:2)-xlower(1:2)
hh(1:2) = p(1:2)/REAL(mx(1:2),KIND=r8)
!
IF(ABS(hh(1)-hh(2)) > 1.0E-10_r8) THEN
  PRINT *, 'Space Step Size Error.  Stopping.'
  STOP
ELSE
  h = hh(1)
END IF
!
cmx(1:2) = mx(1:2)
DO level = 0, minlevel, -1
  IF(MODULO(cmx(1),2) == 1 .OR. MODULO(cmx(2),2) == 1) THEN
    IF(minlevel < level) THEN
      PRINT *, 'Refinement Error.  Stopping.'
      STOP
    END IF
  END IF
  cmx(1:2) = cmx(1:2)/2
END DO
!
IF(numauxvars <= 0) THEN
  auxvars = 1
ELSE
  auxvars = numauxvars
END IF
!
ALLOCATE(  u(0:mx(1)+1,0:mx(2)+1,1:numvars), uo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
         uoo(0:mx(1)+1,0:mx(2)+1,1:numvars),  f(1:mx(1)  ,1:mx(2)  ,1:numvars), &
         res(1:mx(1)  ,1:mx(2)  ,1:numvars),aux(0:mx(1)+1,0:mx(2)+1,1:auxvars))
!
u = 0.0_r8
!
CALL SetProblem
!
IF(restart) THEN
  CALL ReadIn(restartnum)
  startprintout = restartnum+1
ELSE
  CALL Initialize(  u(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                  aux(0:mx(1)+1,0:mx(2)+1,1:auxvars),h,xlower(1:2))
  CALL PrintOut(0)
  time = 0.0_r8
  startprintout = 1
END IF
!
IF(getenergy) OPEN(UNIT=91,FILE='energy.dat', STATUS='REPLACE', ACTION='WRITE')
!
uo = u
!
printiterations = maxtimeiterations/numprintouts
!
DO i = startprintout, numprintouts
  DO k = 1, printiterations
!
    time = REAL(k+(i-1)*printiterations,KIND=r8)*dt
    PRINT *, 'Time = ', time
!
    uoo = uo
    uo = u
    u = 2.0_r8*uo-uoo
!
    CALL SetForcingTerm(f(1:mx(1),1:mx(2),1:numvars),h,xlower(1:2))
!
    CALL NonLinearMultigrid(u,uo,uoo,f,aux,res,h,tol,iter,energyvalue)
!
    IF(getenergy) WRITE(91,'(2F25.12)') time, energyvalue
!
  END DO
!
  CALL PrintOut(i)
!
END DO
!
IF(getenergy) CLOSE(UNIT=91)
!
DEALLOCATE(u,uo,uoo,f,res,aux)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadIn(num)
USE Problem, ONLY: BoundaryConditions
IMPLICIT NONE
!
INTEGER, INTENT(IN):: num
!
CHARACTER(LEN=13):: file_name
CHARACTER(LEN=4):: string
INTEGER:: i, j, threenumvars
INTEGER, DIMENSION(1:2):: mmx
REAL(KIND=r8):: dummy, ttime, x, y
REAL(KIND=r8), DIMENSION(1:2):: hh, xxlower, xxupper
!
WRITE(string,'(i4)') num
DO i = 1, 4
  IF(string(i:i)==' ') string(i:i) = '0'
END DO
file_name = 'OUT/m'//string//'.dat'
!
OPEN(UNIT=9, FILE=file_name, STATUS='UNKNOWN', ACTION='READ')
READ(9,'(F25.12)') time
READ(9,'(I8)') threenumvars
READ(9,'(F25.12)') hh(1)
READ(9,'(F25.12)') hh(2)
READ(9,'(F25.12)') xxlower(1)
READ(9,'(F25.12)') xxlower(2)
READ(9,'(F25.12)') xxupper(1)
READ(9,'(F25.12)') xxupper(2)
READ(9,'(I8)') mmx(1)
READ(9,'(I8)') mmx(2)
!
IF(threenumvars /= 3*numvars) THEN 
  PRINT *, 'threenumvars /= 3*numvars'
  STOP
END IF
IF(ABS(hh(1)-h) > 1.0E-10_r8 .OR. ABS(hh(2)-h) > 1.0E-10_r8) THEN 
  PRINT *, 'hh /= h'
  STOP
END IF
IF(ANY(xxlower /= xlower)) THEN 
  PRINT *, 'xxlower /= xlower'
  STOP
END IF
IF(ANY(xxupper /= xupper)) THEN 
  PRINT *, 'xxupper /= xupper'
  STOP
END IF
IF(ANY(mmx /= mx)) THEN 
  PRINT *, 'mmx /= mx'
  STOP
END IF
!
DO j = 1, mx(2)
  DO i = 1, mx(1)
    READ(9,'(14(f25.12,1x),f25.12)') dummy, dummy, (  u(i,j,k),k=1,numvars), &
                                                   (  f(i,j,k),k=1,numvars), &
                                                   (res(i,j,k),k=1,numvars)
  END DO
END DO
!
CLOSE(9)
!
CALL BoundaryConditions(u(0:mx(1)+1,0:mx(2)+1,1:numvars))
!
END SUBROUTINE ReadIn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE PrintOut(num)
IMPLICIT NONE
!
INTEGER, INTENT(IN):: num
!
CHARACTER(LEN=13):: file_name
CHARACTER(LEN=4):: string
INTEGER:: i, j
REAL(KIND=r8):: x, y
!
WRITE(string,'(i4)') num
DO i = 1, 4
  IF(string(i:i)==' ') string(i:i) = '0'
END DO
file_name = 'OUT/m'//string//'.dat'
!
OPEN(UNIT=9, FILE=file_name, STATUS='REPLACE', ACTION='WRITE')
WRITE(9,'(F25.12)') time
WRITE(9,'(I8)') 3*numvars
WRITE(9,'(F25.12)') h
WRITE(9,'(F25.12)') h
WRITE(9,'(F25.12)') xlower(1)
WRITE(9,'(F25.12)') xlower(2)
WRITE(9,'(F25.12)') xupper(1)
WRITE(9,'(F25.12)') xupper(2)
WRITE(9,'(I8)') mx(1)
WRITE(9,'(I8)') mx(2)
!
DO j = 1, mx(2)
  y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
  DO i = 1, mx(1)
    x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
    WRITE(9,'(14(f25.12,1x),f25.12)') x, y, (  u(i,j,k),k=1,numvars), &
                                            (  f(i,j,k),k=1,numvars), &
                                            (res(i,j,k),k=1,numvars)
  END DO
END DO
!
CLOSE(9)
!
END SUBROUTINE PrintOut
!
END PROGRAM NLMGDriver