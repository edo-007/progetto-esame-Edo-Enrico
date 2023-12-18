MODULE UTILS

    USE CONSTANTS


    CONTAINS


        SUBROUTINE DATA_OUTPUT(timestep, time, x, y, Tn, Te)
            !------------------------------------------------------------!
            IMPLICIT NONE
            !------------------------------------------------------------!
            INTEGER,            INTENT(IN) :: timestep
            REAL,               INTENT(IN) :: time, x(:), y(:), Tn(:,:), Te(:,:)
            !
            INTEGER                        :: i,j, DataUnit
            CHARACTER(LEN=10)              :: citer
            CHARACTER(LEN=200)             :: IOFileName
            !------------------------------------------------------------!
            !
            WRITE(citer,'(I4.4)') timestep                        ! convert iteration number to string
            IOFileName = TRIM(DirName)//TRIM(TestName)//'-'//TRIM(citer)//'.dat' ! name of output file
            DataUnit   = 100                                      ! unit for output file
            !
            OPEN(UNIT=DataUnit, FILE=TRIM(IOFilename), STATUS='UNKNOWN', ACTION='WRITE')
            !
            
            WRITE(DataUnit,*) 'TITLE = "CURRENT TIME ', time, ' "'   
            WRITE(DataUnit,*) 'ZONE T="Only Zone", I=',IMAX,' J=',JMAX
            WRITE(DataUnit,*) 'VARIABLES = "x" "y" "Tn" "Te"'
            !
            DO i = 1, IMAX
                DO j = 1, JMAX
                    WRITE(DataUnit,*) x(i), y(j), Tn(i,j), Te(i,j)
                ENDDO
            ENDDO  
            !
            CLOSE(DataUnit)
            !
        END SUBROUTINE DATA_OUTPUT


    
        ! Overwrite mat1 with mat2
        SUBROUTINE OVERWRITE_SOLUTION(mat1, mat2)
        
            IMPLICIT NONE

            REAL, ALLOCATABLE ::  mat1(:,:), mat2(:,:)
            INTEGER :: i,j

            DO i=2, IMAX-1
                DO j=1, JMAX
                    mat1(i,j) = mat2(i,j)
                ENDDO
            ENDDO

        END SUBROUTINE OVERWRITE_SOLUTION

        
        SUBROUTINE PRINT_MATRIX(matrix)
            
            IMPLICIT NONE
            
            REAL, ALLOCATABLE  ::  matrix(:,:)
            INTEGER :: i,j

            DO i=1, IMAX
                DO j=1, JMAX
                    WRITE(*,"(F10.1,$)") matrix(i,j)
                END DO
                WRITE (*,*) '' ! Print newline
            END DO

        END SUBROUTINE PRINT_MATRIX

        SUBROUTINE INIT_CONDITIONS(Tn, x, discPoint)
            
            IMPLICIT NONE

            INTEGER :: i,j
            
            REAL    :: Tn(IMAX,JMAX), x(IMAX)
            REAL    :: discPoint
            
            DO i = 1, IMAX
                DO j = 1, JMAX
                    IF(x(i).LE.discPoint) THEN
                        Tn(i,j) = TL
                    ELSE
                        Tn(i,j) = TR
                    ENDIF
                ENDDO  
            ENDDO 
            
        END SUBROUTINE INIT_CONDITIONS


        SUBROUTINE EXACT_SOLUTION(Te, x, kappa, time)
            
            REAL, ALLOCATABLE ::  Te(:,:), x(:) 
            REAL :: kappa, time
            INTEGER :: i,j


            DO i = 1, IMAX 
                DO j =1, JMAX
                    Te(i,j) = ((TR + TL) / 2.0) + &
                            (  ERF( (x(i) - 1) / (2 * sqrt(kappa * time)) ) * (TR - TL) / 2.0  )
                ENDDO
            ENDDO

        END SUBROUTINE EXACT_SOLUTION




        SUBROUTINE SECURE_DEALLOC(x,y, Tn, Tn1, Te)

            INTEGER           :: stat_var
            REAL, ALLOCATABLE :: x(:), y(:), Tn(:,:), Tn1(:,:), Te(:,:)

            IF (ALLOCATED(x)) THEN
                DEALLOCATE(x, STAT=stat_var)
                IF (stat_var.NE.0) THEN
                    WRITE(*,'(a)') ' | Error in deallocation of x '
                ENDIF
            ENDIF
            IF (ALLOCATED(y)) THEN
                DEALLOCATE(y, STAT=stat_var)
                IF (stat_var.NE.0) THEN
                    WRITE(*,'(a)') ' | Error in deallocation of y  '
                ENDIF
            ENDIF
            IF (ALLOCATED(Tn)) THEN
                DEALLOCATE(Tn, STAT=stat_var)
                if (stat_var.NE.0) THEN
                    WRITE(*,'(a)') ' | Error in deallocation of Tn  '
                ENDIF
            ENDIF
            IF (ALLOCATED(Tn1)) THEN
                DEALLOCATE(Tn1,STAT=stat_var)
                if (stat_var.NE.0) THEN
                    WRITE(*,'(a)') ' | Error in deallocation of Tn1  '
                ENDIF
            ENDIF
            IF (ALLOCATED(Te)) THEN
                DEALLOCATE(Te, STAT=stat_var)
                if (stat_var.NE.0) THEN
                    WRITE(*,'(a)') ' | Error in deallocation of Te  '
                ENDIF
            ENDIF

        END SUBROUTINE SECURE_DEALLOC


        SUBROUTINE SECURE_ALLOCATION(x,y, Tn, Tn1, Te)

            REAL, ALLOCATABLE  :: Tn(:,:), Tn1(:,:), Te(:,:), x(:), y(:)

            ALLOCATE( x(IMAX) )
            ALLOCATE( y(JMAX) )
            ALLOCATE( Tn( IMAX, JMAX))
            ALLOCATE( Tn1(IMAX, JMAX))
            ALLOCATE( Te(IMAX, JMAX))


        END SUBROUTINE SECURE_ALLOCATION

END MODULE UTILS