PROGRAM main
    
    IMPLICIT NONE

    CALL Heat2D_OpenMP()


CONTAINS 

    SUBROUTINE Heat2D_OpenMP
        ! ---------------------------------------------------------------------- !
        USE OMP_LIB 
        USE UTILS

        IMPLICIT NONE

        ! Time
        REAL    :: time, tend, tio, dtio, dt
        REAL    :: dx, dx2, dy, dy2
        REAL    :: exact

        INTEGER :: i, j, n
        !
        REAL, ALLOCATABLE  :: Tn(:,:), Tn1(:,:), Te(:,:), x(:), y(:)
    
        INTEGER, PARAMETER :: NCPU = 4        ! number of threads for parallel simulations
        
        WRITE(*,'(a)') ' | ====================================================== | '
        WRITE(*,'(a)') ' |      Finite difference code for 2D heat equation       | '
        WRITE(*,'(a)') ' | ====================================================== | '
        WRITE(*,'(a)') ' | '

        ! INTEGER            :: NPRCS           ! number of available threads
        ! ---------------------------------------------------------------------- !

        !$ WRITE(*,*) ' Parallel simulation with OpenMP directives. '
        ! NPRCS = OMP_GET_NUM_PROCS()
        ! CALL OMP_SET_NUM_THREADS(NCPU)
        !$ WRITE(*,*) ' Total number of available CPUs: ', NPRCS
        !$ WRITE(*,*) ' Total number of used CPUs:      ', NCPU
        
        
        ! WRITE(*,*) ' EXPLICIT discretization of the 1D heat equation. '
        ! WRITE(*,*) ' CFL number = ', CFL
        

        ! Boundary conditions
        time = 0.0
        tend = 0.05
        tio  = 0.0
        dtio = 0.01
    

        ! Matrix allocation
        
        WRITE(*,'(a)') ' | Building the computational domain... '
        CALL SECURE_ALLOCATION(x,y, Tn, Tn1, Te)
        
        dx = (xR-xL)/REAL(IMAX-1)
        dy = (yT-yB)/REAL(JMAX-1)
        dx2= dx**2
        dy2= dy**2

        ! Domain definition on dimension x and y
        x(1) = xL
        DO i = 1, IMAX-1
            x(i+1) = x(i) + dx
        ENDDO  
        y(1) = yB
        DO i = 1, JMAX-1
            y(i+1) = y(i) + dy
        ENDDO  

        WRITE(*,'(a)') ' | Assigning initial condition... '
        CALL INIT_CONDITIONS(Tn,x,xD) 
        
        WRITE(*,'(a)') ' | '
        WRITE(*,'(a)') ' | Explicit finite difference solver. '
        WRITE(*,'(a)') ' | '
        WRITE(*,'(a)') ' | START of the COMPUTATION '


        DO n = 1, NMAX   ! main loop in time
        
            ! Exit condition
            IF(time.GE.tend) THEN
                EXIT
            ENDIF

            dt = CFL*d*MAX(dx2, dy2)/kappa
            IF( (time+dt).GT.tend ) THEN
                dt = tend-time  ! adjust the last time step in order to exactly match tend
                tio = tend
            ENDIF 

#ifdef _OUTPUT
            IF( ( time + dt ).GT.tio) THEN
                dt = tio - time;
            ENDIF
#endif
    
        ! EXPLICIT finite difference method
            
            DO i = 2, IMAX-1

                DO j = 2, JMAX-1

                    Tn1(i,j) = Tn(i,j) + kappa*dt/dx2*( Tn(i + 1, j) - 2.0 * Tn(i,j) + Tn(i - 1,j) ) + &
                                         kappa*dt/dy2*( Tn(i, j + 1) - 2.0 * Tn(i,j) + Tn(i,j - 1) )

                ENDDO
                exact            = ( (TR + TL) / 2.0)  +  ( ERF((x(i) - 1) / (2 * SQRT(kappa * time))) * (TR - TL) / 2.0 );    
                Tn1(i, 1)        = exact
                Tn1(i, JMAX)     = exact                
            ENDDO 

            ! Update time and current solution
            time = time + dt
            CALL OVERWRITE_SOLUTION(Tn, Tn1)

#ifdef _OUTPUT
            IF ( (ABS(time - tio) < 1e-12) ) THEN

                WRITE(*,'(a,f15.7)') ' |   plotting data output at time ', time
                CALL EXACT_SOLUTION(Te, x, kappa, time)
                CALL DATA_OUTPUT(n, time, x, y, Tn, Te)
                tio = tio + dtio
            ENDIF
#endif

            
        ENDDO !n

        WRITE(*,'(a)') ' | END OF THE COMPUTATION '
        
        ! Empty memory
        CALL SECURE_DEALLOC(x, y, Tn, Tn1, Te)

        WRITE(*,'(a)') ' | '
        WRITE(*,'(a)') ' |         Finalization was successful. Bye :-)           | '
        WRITE(*,'(a)') ' | ====================================================== | '
        
        
    END SUBROUTINE Heat2D_OpenMP


END PROGRAM main



