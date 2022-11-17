#include "calculate_func.h"

void mcstepPI(int nSteps, double *yim1input, double *mean, double *var, unsigned long long *nStepstot) {
  int i;
  unsigned long long n0 = nStepstot[0];
  unsigned long long n =0;
  double x = 0;
  double yi = 0, yim1 = yim1input[0];
  double meani = 0, meanim1 = mean[0];
  double vari = 0, varim1 = var[0];
  for (i=1; i<=nSteps; ++i) {
    n = n0 + i;
    x = 1.0 * rand() / RAND_MAX;
    yi = sqrt(1 - x * x);
    meani = meanim1 + (yi - meanim1)/n;
    vari  = varim1  + (yi - meanim1) * (yi - meani);
    yim1 = yi;
    meanim1 = meani;
    varim1  = vari; 
    //printf(" vari=%10.5f ( yi-meani)=%10.5f meani=%10.5f\n", vari/(i-1),yi-meani, meani);
    //printf(" meani=%10.5f mean=%10.5f\n", meani, sumy/i);
  }
  var[0] = vari;
  mean[0] = meani;
  yim1input[0] = yim1;
}
