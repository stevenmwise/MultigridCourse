MODULE LinearMultigridRoutines
USE Global
IMPLICIT NONE
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE LinearMultigrid(u,f,m,residualvector,h,tol,iter,energyminimum)
USE Global
USE Utilities, ONLY: Norm
USE Problem, ONLY: Energy, EnergyGradient
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: m
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
!  PRINT *, 'Multigrid iteration', its
!
  iter = its
!
  CALL VCycle(u(0:mx(1)+1,0:mx(2)+1),f(1:mx(1),1:mx(2)),m(0:mx(1)+1,0:mx(2)+1),h,0)
!
  residualvector(1:mx(1),1:mx(2)) = EnergyGradient(u(0:mx(1)+1,0:mx(2)+1), &
                                                   f(1:mx(1)  ,1:mx(2)  ), &
                                                   m(0:mx(1)+1,0:mx(2)+1),h)
!
  energyminimum = Energy(u(0:mx(1)+1,0:mx(2)+1),f(1:mx(1),1:mx(2)), &
                         m(0:mx(1)+1,0:mx(2)+1),h)
!
  residual = Norm(residualvector)
!
  PRINT *, 'Multigrid iteration', its, '  Residual l2   =', residual
!  PRINT *, 'Residual linf =', MAXVAL(ABS(residualvector))
!  PRINT *, 'Energy value  =', energyminimum
!  PRINT *, ' '
!  PRINT *, ' '
!
  IF(residual < tol) EXIT

END DO
!
END SUBROUTINE LinearMultigrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
RECURSIVE SUBROUTINE VCycle(u,f,m,h,level)
USE Global
USE Utilities, ONLY: Restriction, BiLinProlongation, Prolongation
USE Problem, ONLY: EnergyGradient, Relaxation, BoundaryConditions
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: m
REAL(KIND=r8), INTENT(IN):: h
INTEGER, INTENT(IN):: level
!
INTEGER, DIMENSION(1:2):: cmx, mx
REAL(KIND=r8):: ch
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: cf, cu, cm
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
CALL Relaxation(u,f,m,h,presmooth)
!
IF(level>minlevel) THEN
!
  ch = 2.0_r8*h
  cmx(1:2) = mx(1:2)/2
!
  ALLOCATE(cf(1:cmx(1),1:cmx(2)),cu(0:cmx(1)+1,0:cmx(2)+1),cm(0:cmx(1)+1,0:cmx(2)+1))
!
  cu(0:cmx(1)+1,0:cmx(2)+1) = 0.0_r8
  cm(1:cmx(1),1:cmx(2)) = Restriction(m(1:mx(1),1:mx(2)))
!
! Calculate the coarse loading function cf:
  cf(1:cmx(1),1:cmx(2)) &
    = Restriction(-EnergyGradient(u(0:mx(1)+1,0:mx(2)+1), &
                                  f(1:mx(1)  ,1:mx(2)  ), &
                                  m(0:mx(1)+1,0:mx(2)+1),h))
!
  CALL VCycle(cu,cf,cm,ch,level-1)
!
!  u(0:mx(1)+1,0:mx(2)+1) = u(0:mx(1)+1,0:mx(2)+1) &
!                         + BiLinProlongation(cu(0:cmx(1)+1,0:cmx(2)+1))
!
  u(1:mx(1),1:mx(2)) = u(1:mx(1),1:mx(2)) &
                     + Prolongation(cu(1:cmx(1),1:cmx(2)))
!
  CALL Relaxation(u,f,m,h,postsmooth)
!
  DEALLOCATE(cf,cu,cm)
!
END IF
!
END SUBROUTINE VCycle
!
END MODULE LinearMultigridRoutines