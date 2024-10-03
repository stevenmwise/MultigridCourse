MODULE LinearMultigridRoutines
USE Global
IMPLICIT NONE
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION LinearMultigrid(u,f,hf,tol,iter) RESULT(uApproxResult)
USE Global
USE Utilities, ONLY: Norm
USE Problem, ONLY: Residual
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: hf
REAL(KIND=r8), INTENT(IN):: tol
INTEGER, INTENT(OUT):: iter
!
INTEGER:: its
INTEGER, DIMENSION(1:2):: nf
REAL(KIND=r8):: errEst, residualValue
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-1,0:SIZE(u,2)-1):: uApproxResult
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
!
DO its = 1, maxItrs
!
  iter = its
!
  u = MGOperator(u,f,hf,maxLevel,errEst)
!
  residualValue = Norm(Residual(u,f,hf))
!
  PRINT 1001, its, errEst, residualValue
  1001 FORMAT('MG iteration =', i3, 3X, 'errEst =', ES15.7, 3X, &
              'residualValue =' ES15.7)
!
  IF(errEst < tol) EXIT
!
END DO
!
uApproxResult = u
!
END FUNCTION LinearMultigrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
RECURSIVE FUNCTION MGOperator(u,f,hf,level,errEst) RESULT(uApproxResult)
USE Global
USE Utilities, ONLY: Restriction, BiLinProlongation, Norm, Prolongation
USE Problem, ONLY: Residual, Smooth, BoundaryConditions
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: hf
INTEGER, INTENT(IN):: level
REAL(KIND=r8), INTENT(OUT):: errEst
!
INTEGER:: sig
INTEGER, DIMENSION(1:2):: nc, nf
REAL(KIND=r8):: hc
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: gc, cGc
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-1,0:SIZE(u,2)-1):: uApproxResult
!
nf(1) = SIZE(u,1)-2; nf(2) = SIZE(u,2)-2
!
u = Smooth(u,f,hf,presmooth)
!
errEst = 1.0_r8
!
IF(level > 0) THEN
!
  hc = 2.0_r8*hf
  nc(1:2) = nf(1:2)/2
!
  ALLOCATE(gc(1:nc(1),1:nc(2)),cGc(0:nc(1)+1,0:nc(2)+1))
!
  cGc(0:nc(1)+1,0:nc(2)+1) = 0.0_r8
!
! Calculate the coarse loading function gc:
  gc(1:nc(1),1:nc(2)) &
    = Restriction(Residual(u(0:nf(1)+1,0:nf(2)+1), &
                           f(1:nf(1)  ,1:nf(2)  ),hf))
!
! pCycle = 1 is VCycle, pCycle = 2 is WCycle:
  DO sig = 1, pCycle
    cGc = MGOperator(cGc,gc,hc,level-1,errEst)
  END DO
!
! Add the coarse-grid correction using bilinear prolongation:
!  u(0:nf(1)+1,0:nf(2)+1) = u(0:nf(1)+1,0:nf(2)+1) &
!                         + BiLinProlongation(cgc(0:nc(1)+1,0:nc(2)+1))
!
! Add the coarse-grid correction using injection prolongation:
  u(1:nf(1),1:nf(2)) = u(1:nf(1),1:nf(2)) &
                     + Prolongation(cGc(1:nc(1),1:nc(2)))
!
  IF(level == maxLevel) errEst = Norm(cGc(1:nc(1),1:nc(2)))
!
  u = Smooth(u,f,hf,postsmooth)
!
  DEALLOCATE(gc,cGc)
!
END IF
!
uApproxResult = u
!
END FUNCTION MGOperator
!
END MODULE LinearMultigridRoutines