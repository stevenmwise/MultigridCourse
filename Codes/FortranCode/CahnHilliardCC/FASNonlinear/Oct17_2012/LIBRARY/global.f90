MODULE Global
IMPLICIT NONE
!
LOGICAL:: getenergy
INTEGER, PARAMETER:: r8 = SELECTED_REAL_KIND(15,307)
INTEGER:: auxvars, minlevel, maxits, numauxvars, numvars, postsmooth, &
          presmooth
REAL(KIND=r8), PARAMETER:: pi = 3.1415926535897932_r8
REAL(KIND=r8):: dt, omega, time
!
END MODULE Global