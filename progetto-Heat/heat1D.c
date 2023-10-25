#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
     
#define NMAX 1e6  // maximum number of time steps 
#define IMAX 500  // total number of cells
#define DAT_PATH "heat-dat-files"

char* trim(char *s) {

    char *ptr;
    if ( !s )
        return NULL;   // handle NULL string
    if ( !(*s) )
        return s;      // handle empty string

    // //  If an argument (character) passed to the isspace() function is a white-space character, 
    // //  it returns non-zero integer. If not, it returns 0.
    for (ptr = s + strlen(s)-1; (ptr >= s) && isspace(*ptr); --ptr);

    // Sets the first num bytes of the block of memory pointed by ptr to the specified valu).
    ptr[1] = '\0';
    return s;
}

void DataOutput( int timestep, char testname[200], int imax, float time, float *x, float *T ) {

  /*__________________________________*/
    int i, DataUnit;
    char citer[10];
    char IOFilename[200];
    FILE *fp;
  /*__________________________________*/

    trim( testname );
    sprintf(IOFilename, "%s/%s-%04d.dat", DAT_PATH, testname, timestep);

    fp = fopen( IOFilename, "w+t");
    if ( fp == NULL ) {

        /* Errore nell' apertura del file */
        char error_mex[250];
        sprintf( error_mex, "errore nell' apertura del file %s", IOFilename ); 
        perror( error_mex );

        exit( EXIT_FAILURE );
    }
    /* Scrittura su file */
    fprintf(fp, "CURRENT TIME : %f\n", time );
    fprintf(fp, "VARIABLES : 'x' 'T'\n" );
    fprintf(fp, "ZONE T='Only Zone', I= %d  F=POINT\n", IMAX  );

    for ( i = 0; i < IMAX; i++ ){
        fprintf(fp, "%f   %f\n",x[i], T[i]);
    }

    fclose(fp);

    

 }


int main(int argc, char **argv){

    float  xL, xR;          // left and right domain definition
    float  xD;              // location of discontinuity
    float  dx, dx2;         // mesh size and square of mesh size
    float *x;               // vertex coords
    
    float CFL;              // Courant-Friedrichs-Lewy number for stability condition (CFL<1) 
    float time;             // current time
    float dt;               // time step
    float dt_fix;           // user defined time step
    float tend;             // final time
    float *T;               // temperature 
    float *T1;          
    
    char TestName[200];     // name of test problem
    float tio;              // output time and time step 
    float dtio;             // output time step
    float d;                // stability parameter  
    float kappa;            // heat conduction coefficient
    float TL, TR;           // boundary conditions

// 
//  ======================================== 
//  &     SETTINGS OF THE COMPUTATION      &
//  ========================================
// 
    strcpy(TestName,"Heat1D_explicit"); 
    d = 0.45;  // for stability condition d<0.5
    
    
    xL     = -1.0;
    xR     = 1.0;
    xD     = 0.0 ; // location of the discontinuity
    
    time   = 0.0;
    tend   = 0.01;
    tio    = 0.0;
    dtio   = 1e-2;
    CFL    = 0.9;
    dt_fix = 1e-2;
    
    kappa  = 1;     // heat conduction coefficient 
    
    // Boundary condition
    TL = 100;
    TR = 50;
    
    // Allocate variables x, T, T1 of dim. IMAX         ( controlli da includere )
    x  = (float*) malloc( IMAX * sizeof ( float ) );
    T  = (float*) malloc( IMAX * sizeof ( float ) );    
    T1 = (float*) malloc( IMAX * sizeof ( float ) );

    //  1) Building computational domain  
    dx = ( xR - xL ) / (float)( IMAX -1 );
    dx2 = dx * dx;
    x[0] = xL;

    int i;
    for (i = 0; i < (IMAX-1); i++ ) {
        x[i+1] = x[i] + dx;
    }

    //  2) Initial condition  
    for (i = 0; i <= IMAX; i++ ){
        if ( x[i] < xD )
            T[i] = TL;
        else
            T[i] = TR;
    }

    // ___ 3) Computation: main loop in time _____________ 

        puts( " | ");
        puts( " | Explicit finite difference solver. ");
        puts( " | ");
        puts( " | START of the COMPUTATION ");

    int n;
    for ( n = 0; n <= NMAX; n++ ) {

        // 3.1) Compute the time step 
        if ( time >= tend ) 
            exit(EXIT_SUCCESS);

        dt = d * dx2 / kappa;
        if ( ( time + dt ) > tend ) {

            dt = tend - time;
            tio = tend;
        
        }
        if ( ( time + dt ) > tio ) {
            dt = tio - time;
        } 

        /* 3.2) Numerical scheme: FTCS */
    
        // EXPLICIT SOLVER
        for ( i = 1; i < (IMAX-1); i++ ) {
            T1[i] = T[i] + kappa * (dt/dx2) *(T[i-1] - 2*T[i] + T[i+1]);
        }
        T1[0] = TL;
        T1[IMAX] = TR;

        time = time + dt;   // Update time
        T = T1;             // Overwrite current Solution
         
        if ( fabs( time - tio ) < 1e-12 ){
            printf(" | \tPlotting data output at time %.7f\n", time);
            DataOutput(n, TestName, IMAX,time,x,T);
            tio += dtio;
        } 
    }
    return 0;
}