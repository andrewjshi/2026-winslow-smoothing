#include <math.h>

/******************************************************************************/

void flux2d (
  double U[2],
  double Q[4],
  double pars[3],
  double FF[4])
{
  double t11;
  double t2;
  double t3;
  double t5;
  double t6;
  t2 = Q[0];
  t3 = Q[3];
  t5 = pars[0] * (t2 + t3);
  t6 = pars[1];
  t11 = 0.1e1 / pars[2];
  FF[0] = -(0.2e1 * t2 * t6 + t5) * t11;
  FF[1] = -0.2e1 * t6 * (Q[1] / 0.2e1 + Q[2] / 0.2e1) * t11;
  FF[2] = FF[1];
  FF[3] = -(0.2e1 * t3 * t6 + t5) * t11;
  return;
}

/******************************************************************************/

void flux2d_deriv (
  double U[2],
  double Q[4],
  double pars[3],
  double FF[4],
  double dFFdU[8],
  double dFFdQ[16])
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t5;
  double t6;
  t1 = pars[0];
  t2 = Q[0];
  t3 = Q[3];
  t5 = t1 * (t2 + t3);
  t6 = pars[1];
  t11 = 0.1e1 / pars[2];
  FF[0] = -(0.2e1 * t2 * t6 + t5) * t11;
  FF[1] = -0.2e1 * t6 * (Q[1] / 0.2e1 + Q[2] / 0.2e1) * t11;
  FF[2] = FF[1];
  FF[3] = -(0.2e1 * t3 * t6 + t5) * t11;
  dFFdU[0] = 0.0e0;
  dFFdU[1] = 0.0e0;
  dFFdU[2] = 0.0e0;
  dFFdU[3] = 0.0e0;
  dFFdU[4] = 0.0e0;
  dFFdU[5] = 0.0e0;
  dFFdU[6] = 0.0e0;
  dFFdU[7] = 0.0e0;
  dFFdQ[0] = -(t1 + 0.2e1 * t6) * t11;
  dFFdQ[1] = 0.0e0;
  dFFdQ[2] = 0.0e0;
  dFFdQ[3] = -t6 * t11;
  dFFdQ[4] = 0.0e0;
  dFFdQ[5] = -t1 * t11;
  dFFdQ[6] = dFFdQ[3];
  dFFdQ[7] = 0.0e0;
  dFFdQ[8] = 0.0e0;
  dFFdQ[9] = dFFdQ[6];
  dFFdQ[10] = dFFdQ[5];
  dFFdQ[11] = 0.0e0;
  dFFdQ[12] = dFFdQ[9];
  dFFdQ[13] = 0.0e0;
  dFFdQ[14] = 0.0e0;
  dFFdQ[15] = dFFdQ[0];
  return;
}
