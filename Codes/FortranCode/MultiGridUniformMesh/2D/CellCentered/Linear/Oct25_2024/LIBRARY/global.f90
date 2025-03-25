MODULE Global
IMPLICIT NONE
!
INTEGER, PARAMETER:: r8 = SELECTED_REAL_KIND(15,307)
INTEGER:: maxLevel, maxItrs, postSmooth, preSmooth, pCycle
REAL(KIND=r8):: omega
REAL(KIND=r8), DIMENSION(1:2):: xLower, xUpper
!
END MODULE Global