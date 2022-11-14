#include <stdio.h>

#include <mpi.h>

int main(int argc, char *argv[]) {

    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);

    int left_rank = (rank==0)?(size-1):(rank-1);
    int right_rank = (rank==(size-1))?0:(rank+1)%size;

    int stop_cond_rank;
    MPI_Request stop_cond_request;
    int stop_cond= 0;

    int i, s=0;
    int smax = 1 << 20 + 10 * rank;
    printf( " -- %d -- smax = %d \n", rank, smax);

    while( s <= (1 << 22) )
    {

        /* Compute Here and set stop condition accordingly */
        for ( i = 0; i < (1 << 10); ++i )
            s += 1;

        if (s >= (1 << smax ))
            stop_cond = 1;

        if( stop_cond )
        {
            /* Cancel the left recv */
            MPI_Cancel( &stop_cond_request );
            if( rank != right_rank )
                MPI_Send( &rank, 1, MPI_INT, right_rank, 123, MPI_COMM_WORLD );

            break;
        }

        int did_recv = 0;
        MPI_Test( &stop_cond_request, &did_recv, MPI_STATUS_IGNORE );
        if( did_recv )
        {

            printf(" ---- DID RECV ---- \n");

            MPI_Irecv( &stop_cond_rank, 1, MPI_INT, left_rank, 123, MPI_COMM_WORLD, &stop_cond_request);

            stop_cond = 1;
            MPI_Wait( &stop_cond_request, MPI_STATUS_IGNORE );
            if( right_rank != stop_cond_rank )
                MPI_Send( &stop_cond_rank, 1, MPI_INT, right_rank, 123, MPI_COMM_WORLD );

            break;
        }
    }

    if( stop_cond )
    {
        /* Handle the stop condition */
        printf(" ---- Done ---- \n");
    }
    else
    {
        /* Cleanup */
        MPI_Cancel( &stop_cond_request );
    }

    MPI_Finalize();
    return(0);
}
