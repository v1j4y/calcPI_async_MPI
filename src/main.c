#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <sys/time.h>

#include "calculate_func.h"

# define M_PI           3.14159265358979323846  /* pi */

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
    srand( rank);

int sum = 0;
    int flag = -1, res;
    int stop=-1;
    int work = 0, nSteps=50000;
    double sumPI, sumPI2;
    double data[4];
    MPI_Request request;
    MPI_Request requeststop;
    MPI_Status status;
    MPI_Status statusstop;

    while (1) {

        if(flag != 0)
        {
            MPI_Irecv(&stop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &requeststop);
            flag = 0;
        }

        MPI_Test(&requeststop, &flag, &statusstop);

        if (flag != 0) {
            if (statusstop.MPI_SOURCE != -1)
                sum += 1;
            flag = 0;
        }

        if (stop == 1)
          break;

        sumPI  = 0.0;
        sumPI2 = 0.0;
        mcstepPI(nSteps, &sumPI, &sumPI2);

        data[0] = rank;
        data[1] = sumPI;
        data[2] = sumPI2;
        data[3] = nSteps;
        MPI_Isend(data, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
        sum++;

    }
}

else { // Master
    double sum = 0;
    int flag = -1, res;
    MPI_Request request;
    MPI_Status status;
    double data[4];
    int buf, i, nmax;
    double sumPItot  = 0.0;
    double sumPI2tot = 0.0;
    double nStepstot = 0.0;
    double errorPI = 10.0;
    double varPI   = 10.0;
    nmax = 12;
    while (1) {
        if(flag != 0)
        {
            MPI_Irecv(data, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &request);
            flag = 0;
        }

        MPI_Test(&request, &flag, &status);

        if (flag != 0) {
            if (status.MPI_SOURCE != -1)
                sum += data[0];
            flag = -1;
            sumPItot  += data[1];
            sumPI2tot += data[2];
            nStepstot += data[3];
            errorPI = M_PI/4 - sumPItot/nStepstot;
            varPI = (sumPI2tot - (sumPItot * sumPItot/nStepstot))/(nStepstot-1);
        }

        if (fabs(errorPI) < 0.00000001){
          for( i=1; i < size; ++i) {
            buf = 1;
            MPI_Send(&buf, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
          }
          break;
        }
    }

    printf("%10.9f (Error=%10.9f) Var=%10.9f\n", sumPItot/nStepstot, errorPI, varPI);
}

MPI_Finalize();
return 0;

}
