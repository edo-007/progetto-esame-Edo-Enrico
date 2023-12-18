

PROGRAM hd1omp

    IMPLICIT NONE
    CALL Heat1D_OpenMP()

END PROGRAM hd1omp


SUBROUTINE Heat1D_OpenMP
    ! ---------------------------------------------------------------------- !
    USE OMP_LIB 

    IMPLICIT NONE
    INTEGER :: i, n, timestep
    INTEGER :: IMAX    ! total number of cells
    REAL    :: xL, xR  ! left and right coords of the domain
    REAL    :: dx, dx2 ! mesh space (and its squared value)
    REAL    :: CFL     ! CFL number (<=1 for stability of EXPLICIT SCHEMES)
    REAL    :: time    ! current time
    REAL    :: dt      ! time step
    REAL    :: tend    ! final time of the simulation
    REAL    :: TL, TR  ! boundary conditions
    REAL    :: kappa   ! heat conduction coefficient 
    !
    REAL, ALLOCATABLE  :: T(:), Tnew(:), x(:)
    INTEGER, PARAMETER :: NMAX=1e6  ! max number of time steps
    
  
    INTEGER, PARAMETER :: NCPU = 4        ! number of threads for parallel simulations
    INTEGER            :: NPRCS           ! number of available threads
    ! ---------------------------------------------------------------------- !
    
      !$ WRITE(*,*) ' Parallel simulation with OpenMP directives. '
      NPRCS = OMP_GET_NUM_PROCS()
      CALL OMP_SET_NUM_THREADS(NCPU)
      !$ WRITE(*,*) ' Total number of available CPUs: ', NPRCS
      !$ WRITE(*,*) ' Total number of used CPUs:      ', NCPU
      CFL = 0.9
      WRITE(*,*) ' EXPLICIT discretization of the 1D heat equation. '
      WRITE(*,*) ' CFL number = ', CFL
  
    ! Domain definition
    IMAX = 5000
    xL   = -1.0
    xR   =  1.0
    ! Boundary conditions
    TL   = 100.
    TR   = 50.
    
    time = 0.0
    tend = 0.02
    
    kappa= 1.0
    
    ALLOCATE( x(IMAX)    )
    ALLOCATE( T(IMAX)    )
    ALLOCATE( Tnew(IMAX) )
  
    
    WRITE(*,*) ' Building the computational domain... '
    
    dx = (xR-xL)/REAL(IMAX-1)
    dx2= dx**2
    
    x(1) = xL
    DO i = 1, IMAX-1
      x(i+1) = x(i) + dx
    ENDDO  
    
    WRITE(*,*) ' Assigning the initial condition... '
    
    DO i = 1, IMAX
      IF(x(i).LE.0.0) THEN
        T(i) = TL
      ELSE
        T(i) = TR
      ENDIF  
    ENDDO  
    
    WRITE(*,*) ' START OF THE COMPUTATION '
  
  !$OMP PARALLEL
    DO n = 1, NMAX   ! main loop in time
      
      IF(time.GE.tend) EXIT
    
  !$OMP SINGLE
      dt = 0.5*CFL*dx2/kappa
      IF(time+dt.GT.tend) THEN
        dt = tend-time  ! adjust the last time step in order to exactly match tend
      ENDIF 
  
  !$OMP END SINGLE
  !$OMP BARRIER
      
      ! EXPLICIT finite difference method
      !$OMP DO
 
      DO i = 2, IMAX-1
        Tnew(i) = T(i) + kappa*dt/dx2 * (T(i+1) - 2.*T(i) + T(i-1)) 
      ENDDO 
      !$OMP END DO
      Tnew(1)    = TL
      Tnew(IMAX) = TR
      
  
      ! Update time and current solution
      !$OMP SINGLE
      time = time + dt
      T(:) = Tnew(:)
      timestep = n
      !$OMP END SINGLE
  
    ENDDO !n  
    !$OMP END PARALLEL
    WRITE(*,*) ' END OF THE COMPUTATION '
      
    ! Empty memory
    DEALLOCATE( T, Tnew, x )
    
  END SUBROUTINE Heat1D_OpenMP
    


