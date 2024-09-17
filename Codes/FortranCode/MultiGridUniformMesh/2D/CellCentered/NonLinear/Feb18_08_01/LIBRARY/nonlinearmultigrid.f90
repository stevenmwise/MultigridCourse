MODULE NonLinearMultigridRoutines
USE Global
IMPLICIT NONE
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE NonLinearMultigrid(u,f,residualvector,h,tol,iter,energyminimum)
USE Global
USE Utilities, ONLY: Norm
USE Problem, ONLY: Energy, EnergyGradient
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(OUT):: residualvector
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), INTENT(IN):: tol
INTEGER, INTENT(OUT):: iter
REAL(KIND=r8), INTENT(OUT):: energyminimum
!
INTEGER:: its
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: residual
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
DO its = 1, maxits
!
  PRINT *, 'Multigrid iteration', its
!
  iter = its
!
  CALL FASVCycle(u(0:mx(1)+1,0:mx(2)+1),f(1:mx(1),1:mx(2)),h,0)
!
  residualvector(1:mx(1),1:mx(2)) = EnergyGradient(u(0:mx(1)+1,0:mx(2)+1), &
                                                   f(1:mx(1)  ,1:mx(2)  ),h)
!
  energyminimum = Energy(u(0:mx(1)+1,0:mx(2)+1),f(1:mx(1),1:mx(2)),h)
!
  residual = Norm(residualvector)
!
  PRINT *, 'Residual l2   =', residual
  PRINT *, 'Residual linf =', MAXVAL(residualvector)
  PRINT *, 'Energy value  =', energyminimum
  PRINT *, ' '
  PRINT *, ' '
!
  IF(residual < tol) EXIT

END DO
!
END SUBROUTINE NonLinearMultigrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
RECURSIVE SUBROUTINE FASVcycle(u,f,h,level)
USE Global
USE Utilities, ONLY: Restriction, Prolongation, BiLinProlongation
USE Problem, ONLY: BoundaryConditions, EnergyGradient, Operator, Relaxation
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: h
INTEGER, INTENT(IN):: level
!
INTEGER, DIMENSION(1:2):: cmx, mx
REAL(KIND=r8):: ch
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: cf, cs, cu
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
CALL Relaxation(u,f,h,presmooth)
!
IF(level>minlevel) THEN
!
  ch = 2.0_r8*h
  cmx(1:2) = mx(1:2)/2
  ALLOCATE(cf(1:cmx(1)  ,1:cmx(2)  ),cs(1:cmx(1)  ,1:cmx(2)  ), &
           cu(0:cmx(1)+1,0:cmx(2)+1))
!
  cu(1:cmx(1),1:cmx(2)) = Restriction(u(1:mx(1),1:mx(2)))
  cs(1:cmx(1),1:cmx(2)) = cu(1:cmx(1),1:cmx(2))
!
  CALL BoundaryConditions(cu)
!
! Calculate the coarse loading function cf:
  cf(1:cmx(1),1:cmx(2)) &
    = Restriction(-EnergyGradient(u(0:mx(1)+1,0:mx(2)+1), &
                                  f(1:mx(1)  ,1:mx(2)  ),h)) &
    + Operator(cu(0:cmx(1)+1,0:cmx(2)+1),ch)
!
  CALL FASVcycle(cu,cf,ch,level-1)
!
  cu(1:cmx(1),1:cmx(2)) = cu(1:cmx(1),1:cmx(2))-cs(1:cmx(1),1:cmx(2))
!
  CALL BoundaryConditions(cu)
  u(0:mx(1)+1,0:mx(2)+1) = u(0:mx(1)+1,0:mx(2)+1) &
                         + BiLinProlongation(cu(0:cmx(1)+1,0:cmx(2)+1))
!
!  u(1:mx(1),1:mx(2)) = u(1:mx(1),1:mx(2)) &
!                     + Prolongation(cu(1:cmx(1),1:cmx(2)))
!
  CALL Relaxation(u,f,h,postsmooth)
!
  DEALLOCATE(cf,cu,cs)
!
END IF
!
END SUBROUTINE FASVcycle
!
END MODULE NonLinearMultigridRoutines