MODULE NonLinearMultigridRoutines
USE Global
IMPLICIT NONE
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE NonLinearMultigrid(u,uo,uoo,f,aux,residualvector,h,tol,iter,energyvalue)
USE Global
USE Utilities, ONLY: Norm
USE Problem, ONLY: Energy, SetSource!, SetAux
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: aux
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(OUT):: residualvector
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), INTENT(IN):: tol
INTEGER, INTENT(OUT):: iter
REAL(KIND=r8), INTENT(OUT):: energyvalue
!
INTEGER:: its
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: res
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2,1:numvars):: s
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
DO its = 1, maxits
!
!  PRINT *, ' '
!  PRINT *, 'Multigrid iteration', its
!
!  CALL SetAux(  u(0:mx(1)+1,0:mx(2)+1,1:numvars), &
!               uo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
!              uoo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
!                f(1:mx(1)  ,1:mx(2)  ,1:numvars), &
!              aux(0:mx(1)+1,0:mx(2)+1,1:auxvars),h)
!
  CALL SetSource(  u(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                  uo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                 uoo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                   f(1:mx(1)  ,1:mx(2)  ,1:numvars), &
                 aux(0:mx(1)+1,0:mx(2)+1,1:auxvars), &
                   s(1:mx(1)  ,1:mx(2)  ,1:numvars),h)
!
  iter = its
!
  CALL FASVCycle(  u(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                  uo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                 uoo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                   s(1:mx(1)  ,1:mx(2)  ,1:numvars), &
                 aux(0:mx(1)+1,0:mx(2)+1,1:auxvars),h,0)
!
    residualvector(1:mx(1)  ,1:mx(2)  ,1:numvars) &
    = Residual(  u(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                uo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
               uoo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
                 s(1:mx(1)  ,1:mx(2)  ,1:numvars), &
               aux(0:mx(1)+1,0:mx(2)+1,1:auxvars),h)
!
  res = Norm(residualvector)
!
  PRINT *, 'MG iteration', its, 'Residual l2   =', res
!  PRINT *, 'Residual linf =', MAXVAL(residualvector)
!
  IF(res < tol) THEN
!    IF(getenergy) energyvalue = Energy(  u(0:mx(1)+1,0:mx(2)+1,1:numvars), &
!                                        uo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
!                                         f(1:mx(1)  ,1:mx(2)  ,1:numvars), &
!                                       aux(0:mx(1)+1,0:mx(2)+1,1:auxvars),h)
    EXIT
  END IF
!
END DO
!
END SUBROUTINE NonLinearMultigrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
RECURSIVE SUBROUTINE FASVCycle(u,uo,uoo,s,aux,h,level)
USE Global
USE Utilities, ONLY: Restriction2D, Prolongation2D
USE Problem, ONLY: BoundaryConditions, BoundaryConditionsAux, Relaxation, &
                   Operator
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: s
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), INTENT(IN):: h
INTEGER, INTENT(IN):: level
!
INTEGER, DIMENSION(1:2):: cmx, mx
REAL(KIND=r8):: ch
REAL(KIND=r8), DIMENSION(:,:,:), ALLOCATABLE:: caux, cs, cu, cus, cuo, cuoo
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
CALL Relaxation(u,uo,uoo,s,aux,h,presmooth,level)
!
IF(level>minlevel) THEN
!
  ch = 2.0_r8*h
  cmx(1:2) = mx(1:2)/2
!
  ALLOCATE(caux(0:cmx(1)+1,0:cmx(2)+1,1:auxvars), &
             cs(1:cmx(1)  ,1:cmx(2)  ,1:numvars), &
             cu(0:cmx(1)+1,0:cmx(2)+1,1:numvars), &
            cus(0:cmx(1)+1,0:cmx(2)+1,1:numvars), &
            cuo(0:cmx(1)+1,0:cmx(2)+1,1:numvars), &
           cuoo(0:cmx(1)+1,0:cmx(2)+1,1:numvars))
!
  caux(1:cmx(1),1:cmx(2),1:auxvars) = Restriction2D(aux(1:mx(1),1:mx(2),1:auxvars))
    cu(1:cmx(1),1:cmx(2),1:numvars) = Restriction2D(  u(1:mx(1),1:mx(2),1:numvars))
   cuo(1:cmx(1),1:cmx(2),1:numvars) = Restriction2D( uo(1:mx(1),1:mx(2),1:numvars))
  cuoo(1:cmx(1),1:cmx(2),1:numvars) = Restriction2D(uoo(1:mx(1),1:mx(2),1:numvars))
!
  CALL BoundaryConditionsAux(caux(0:cmx(1)+1,0:cmx(2)+1,1:auxvars))
  CALL BoundaryConditions(     cu(0:cmx(1)+1,0:cmx(2)+1,1:numvars))
  CALL BoundaryConditions(    cuo(0:cmx(1)+1,0:cmx(2)+1,1:numvars))
  CALL BoundaryConditions(   cuoo(0:cmx(1)+1,0:cmx(2)+1,1:numvars))
!
  cus(1:cmx(1),1:cmx(2),1:numvars) = cu(1:cmx(1),1:cmx(2),1:numvars)
!
! Calculate the coarse source function cs.
  cs(1:cmx(1),1:cmx(2),1:numvars) &
    = Restriction2D(Residual(   u(0: mx(1)+1,0: mx(2)+1,1:numvars), &
                               uo(0: mx(1)+1,0: mx(2)+1,1:numvars), &
                              uoo(0: mx(1)+1,0: mx(2)+1,1:numvars), &
                                s(1: mx(1)  ,1: mx(2)  ,1:numvars), &
                              aux(0: mx(1)+1,0: mx(2)+1,1:auxvars),h)) &
    +               Operator(  cu(0:cmx(1)+1,0:cmx(2)+1,1:numvars), &
                              cuo(0:cmx(1)+1,0:cmx(2)+1,1:numvars), &
                             cuoo(0:cmx(1)+1,0:cmx(2)+1,1:numvars), &
                             caux(0:cmx(1)+1,0:cmx(2)+1,1:auxvars),ch)
!
  CALL FASVCycle(cu,cuo,cuoo,cs,caux,ch,level-1)
!
  cu(1:cmx(1),1:cmx(2),1:numvars) =  cu(1:cmx(1),1:cmx(2),1:numvars) &
                                  - cus(1:cmx(1),1:cmx(2),1:numvars)
!
  u(1:mx(1),1:mx(2),1:numvars) =                 u(1: mx(1),1: mx(2),1:numvars) &
                               + Prolongation2D(cu(1:cmx(1),1:cmx(2),1:numvars))
!
  CALL Relaxation(u,uo,uoo,s,aux,h,postsmooth,level)
!
  DEALLOCATE(caux,cs,cu,cus,cuo,cuoo)
!
END IF
!
END SUBROUTINE FASVCycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Residual(u,uo,uoo,s,aux,h) RESULT(residualresult)
USE Global
USE Problem, ONLY: Operator
IMPLICIT NONE
!
! Operator residual vector:
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: s
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2,1:numvars):: residualresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
  residualresult(1:mx(1)  ,1:mx(2)  ,1:numvars) &
  =            s(1:mx(1)  ,1:mx(2)  ,1:numvars) &
  - Operator(  u(0:mx(1)+1,0:mx(2)+1,1:numvars), &
              uo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
             uoo(0:mx(1)+1,0:mx(2)+1,1:numvars), &
             aux(0:mx(1)+1,0:mx(2)+1,1:auxvars),h)
!
END FUNCTION Residual
!
END MODULE NonLinearMultigridRoutines