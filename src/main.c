#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <sys/time.h>

#include "calculate_func.h"

# define M_PI           3.14159265358979323846  /* pi */
# define M_WINDOW       10                      /* data window */

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
    int work = 0;
    int nSteps=50000;
    unsigned long long nStepstot=0;
    double sumPI = 0, var = 0, mean = 0, yim1 = 0;
    double data[4];
    MPI_Request request;
    MPI_Request requeststop;
    MPI_Status status;
    MPI_Status statusstop;

    MPI_Recv(&nSteps, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);

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

        mcstepPI(nSteps, &yim1, &mean, &var, &nStepstot);
        nStepstot += nSteps;
        //mean = sumPI/nStepstot;

        data[0] = rank;
        data[1] = mean;
        data[2] = var/nStepstot;
        data[3] = nStepstot;
        //printf(" rank=%d error=%10.5f\n",rank,M_PI/4 - mean);
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
    int buf, srs=0, i, nmax, idx=0, sumidx=0, j, k;
    double sumPItot  = 0.0;
    double sumPI2tot = 0.0;
    double nStepstot = 0.0;
    double errorPI = 10.0;
    double varPI   = 10.0;
    double meana  =  0.0;
    double vara   =  0.0;
    double nStepsa=  0.0;
    double meanb  =  0.0;
    double varb   =  0.0;
    double nStepsb=  0.0;
    double meanab =  0.0;
    double varab  =  0.0;
    double nStepsab=  0.0;
    //double datamean[M_WINDOW * (size-1)];
    //double datavar[M_WINDOW * (size-1)];
    //double datansteps[M_WINDOW * (size-1)];
    //int startidx[size-1];
    double *datamean;
    double *datavar;
    double *datansteps;
    int *startidx;
    int nItermax=  0;
    nItermax = 10;
    nItermax = nItermax - (nItermax % (size * (size-1)/2));
    //printf(" dim = %ld",(size-1)*sizeof(int));
    //printf(" dim = %ld",M_WINDOW * (size-1)*sizeof(double));
    //printf(" dim = %ld",M_WINDOW * (size-1)*sizeof(double));
    //printf(" dim = %ld\n",M_WINDOW * (size-1)*sizeof(double));
    printf("Nitermax = %d size=%d\n",nItermax,M_WINDOW * (size-1));
    startidx   = (int *)(malloc((size-1)*sizeof(int)));
    printf(" startidx=%d\n",startidx[0]);
    //datamean   = (double *)(malloc(M_WINDOW * (size-1)*sizeof(double)));
    //datavar    = (double *)(malloc(M_WINDOW * (size-1)*sizeof(double)));
    //datansteps = (double *)(malloc(M_WINDOW * (size-1)*sizeof(double)));
    
    int nSteps = 50000;
    for( i=1; i < size; ++i) {
      MPI_Send(&nSteps, 1, MPI_INT, i, 10, MPI_COMM_WORLD);
      //startidx[i-1] = -1;
    }

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
            //srs = (int)floor(data[0])-1;
            srs = floor(data[0])-1;
            //startidx[srs] += 1;
            //j = startidx[srs] % M_WINDOW;
            //datamean[srs * M_WINDOW + j] = data[1];
            //datavar[srs * M_WINDOW + j] = data[2];
            //datansteps[srs * M_WINDOW + j] = data[3];
            meanb   = data[1];
            varb    = data[2];
            nStepsb = data[3];
            //meanb  = data[1]/nStepsb;
            nStepsab = nStepsa + nStepsb;
            nStepstot += nSteps;
            meanab = (nStepsb * meanb + meana * nStepsa)/nStepsab;
            varab  = vara + varb + (nStepsa * nStepsb) * (meana - meanb) * (meana - meanb) / (nStepsab);
            vara    = varab;
            varab  = sqrt(varab / (nStepsab-1));
            //printf(" sum=%10.5f vara = %10.5f varb=%10.5f meana=%10.5f meanb=%10.5f meanab=%10.5f varab=%1.15f\n",sum,vara,varb,meana,meanb,meanab,varab);
            errorPI = M_PI/4 - meanab;
            varPI = varab;
            meana   = meanab;
            nStepsa = nStepsab;
        }

        //sumidx = 0;
        //for(i=0;i<size-1;++i) {
        //  if(startidx[i] > idx) sumidx += 1;
        //}
        //printf(" flag=%d idx=%d sumidx=%d sum=%10.5f srs=%d startidx=%d\n",flag, idx,sumidx,sum,srs,startidx[srs]);
        printf(" startidx=%d srs=%d\n", startidx[srs], srs);
        //startidx[srs] = 0;

        if(fabs(sum - 0.5 * size * (size-1)) < 1e-10 || (sum - 0.5 * size * (size-1)) > size*size) {
          meana = 0.0;
          vara  = -1.0;
          nStepsa = 0;
          //printf(" Setting to 0 ratio=%10.5f\n",nStepstot/nSteps);
          sum = 0.0;
          //printf("%10.9f (Error=%10.9f) Var=%10.9f\n", meanab, errorPI, varab);
        }

        //if (fabs(errorPI) < 0.00000001){
        if (nStepstot/nSteps >= nItermax){
          for( i=1; i < size; ++i) {
            buf = 1;
            MPI_Send(&buf, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
          }
          break;
        }
    }

    printf("%10.9f (Error=%10.9f) Var=%1.9f (%1.9f %1.9f) sum=%10.5f ratio=%10.5f\n", meanab, errorPI, varab, vara, varb, sum, nStepstot/nSteps);
}

MPI_Finalize();
return 0;

}
