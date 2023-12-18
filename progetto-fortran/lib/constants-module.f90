MODULE CONSTANTS
    
    IMPLICIT NONE

    ! Grid dimension
    INTEGER, PARAMETER :: IMAX = 100
    INTEGER, PARAMETER :: JMAX = 100
    INTEGER, PARAMETER :: NMAX = 1e6


    ! Initial temperatures
    REAL, PARAMETER    :: TL   = 100
    REAL, PARAMETER    :: TR   = 50

    ! Asse x
    REAL, PARAMETER    :: xL   =  0.0
    REAL, PARAMETER    :: xR   =  2.0
    REAL, PARAMETER    :: xD   =  1.0

    ! Asse y 
    REAL, PARAMETER    :: yB   =  0.0
    REAL, PARAMETER    :: yT   =  2.0

   
    REAL, PARAMETER    :: kappa = 1.0
    REAL, PARAMETER    :: d     = 0.25
    REAL, PARAMETER    :: CFL   = 0.9

    ! Output file
    CHARACTER(LEN=200), PARAMETER :: DirName = "heat-dat-files/"
    CHARACTER(LEN=200), pARAMETER :: TestName = "Heat2D"
    
END MODULE 
