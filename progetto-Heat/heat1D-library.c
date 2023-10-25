void posix_memalign_all( void *G, void *C, void *D, void *P){


    if ( posix_memalign ( ( void *)G , 4096 , N*N*sizeof(int ) ) != 0 ) {
        perror("ERROR: allocation of G FAILED:" ) ;
        exit( -1 );
    }