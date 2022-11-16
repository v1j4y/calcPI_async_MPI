#include "calculate_func.h"

void mcstepPI(int nSteps, double *sumPI, double *sumPI2) {
  int i;
  double y = 0, x = 0;
  for (i=0; i<nSteps; ++i) {
    x = 1.0 * rand() / RAND_MAX;
    y = sqrt(1 - x * x);
    sumPI[0]  += y;
    sumPI2[0] += y * y;
  }
}
