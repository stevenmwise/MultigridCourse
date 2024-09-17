MODULE Global
IMPLICIT NONE
!
INTEGER, PARAMETER:: r8 = SELECTED_REAL_KIND(15,307)
INTEGER:: minlevel, maxits, postsmooth, presmooth
REAL(KIND=r8):: omega
!
END MODULE Global