#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <sys/time.h>

#include "calculate_func.h"

int main(int argc, char *argv[])
{

int rank, size;
MPI_Status status;

/* Init */
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

if (rank != 0) { // Slaves
    int buf;
    srand(time(NULL) + rank);

int sum = 0;
    int flag = -1, res;
    int work = 0, nSteps=1000;
    double meanPI, varPI;
    MPI_Request request;
    MPI_Status status;
    while (1) {
        if(flag != 0)
        {
            MPI_Irecv(&res, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            flag = 0;
        }

        MPI_Test(&request, &flag, &status);

        if (flag != 0) {
            if (status.MPI_SOURCE != -1)
                sum += res;
            printf("recv : %d, master : %d sum=%d\n", res, status.MPI_SOURCE,sum);
            flag = -1;
        }

        mcstepPI(nSteps, &meanPI, &varPI);

        if (sum == 1)
            break;
    }
    sleep(2.0);

    MPI_Send(&meanPI, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

    printf("done rank = %d sum : %d work = %d meanPI=%10.5f varPI=%10.5f\n", rank, sum, work, meanPI, varPI);
}

else { // Master
    int sum = 0;
    int flag = -1, res;
    MPI_Request request;
    MPI_Status status;
    while (1) {
        if(flag != 0)
        {
            MPI_Irecv(&res, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            flag = 0;
        }

        MPI_Test(&request, &flag, &status);

        if (flag != 0) {
            if (status.MPI_SOURCE != -1)
                sum += res;
            flag = -1;
            printf("recv : %d, slave : %d sum = %d\n", res, status.MPI_SOURCE,sum);
        }


        if (sum == 6)
            break;
    }
    //sleep(4.0);

int buf;
    buf=1;
    MPI_Send(&buf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    buf=1;
    MPI_Send(&buf, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
    
    MPI_Recv(&buf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
    printf("res from 1 : %d\n", buf);
    MPI_Recv(&buf, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
    printf("res from 2 : %d\n", buf);

    printf("sum : %d\n", sum);
}

MPI_Finalize();
return 0;

}
