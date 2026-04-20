#include "dginterface.h"
#include "mathutil.h"
#include "linelast_autogen.cpp"

DGAPPLICATION(dglinelast)

/**************************************************************************/
/* Application properties */

void getappinfo(int D, phys& p, appinfo& app)
{
    app.N = 2;
    app.needq = true;
}

/* No sensors */
double sensvar(int D, darray& u, darray& dsdu, int sensor) {}

/**************************************************************************/
/* Fluxes */

void volflux(int D, darray& F, darray& dFdu, darray& dFdq,
             darray& u, darray& q, phys& p, int sensor)
{
    if (dFdu != 0)
        flux2d_deriv(u,q,p.pars,F,dFdu,dFdq);
    else
        flux2d(u,q,p.pars,F);
}

/**************************************************************************/
/* Sources */

void volsource(int D, darray& S, darray& dSdu, darray& dSdq,
               darray& u, darray& q, phys& p)
{
}

/**************************************************************************/
/* Numerical fluxes at interior boundaries */

void ibndflux(int D, darray& F, darray& dFdu1, darray& dFdu2,
              darray& dFdq, darray& u1, darray& u2, darray& q,
              darray& nn,phys& p, int sensor)
{
}

/**************************************************************************/
/* Boundary conditions for LDG variables */

void qbndcnds(int D, darray& uhat, darray& duhatdu, darray& u,
              darray& nn, phys& p, int bndnbr)
{
    uhat = u;
    if (duhatdu != 0) {
        duhatdu = 0.0;
        duhatdu(0,0) = 1;
        duhatdu(1,1) = 1;
    }
}

/**************************************************************************/
/* Boundary conditions */

void ubndcnds(int D, darray& F, darray& dFdu, darray& dFdq,
              darray& u, darray& q, darray& nn, phys& p,
              int bndnbr, int sensor)
{
    F = 0.0;
    dFdu = 0.0;
    dFdq = 0.0;
}

/**************************************************************************/
/* Boundary Post processing */

void bndposteval(int D, darray& val, darray& dvaldu, darray& dvaldq,
                 darray& u, darray& q, darray& nn, phys& p,
                 int bndnbr, int postqty, int sensor)
{
}
