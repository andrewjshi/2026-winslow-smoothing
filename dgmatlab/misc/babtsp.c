#include <stdlib.h>

int NN;

#define w(i,j) w[(i)+NN*(j)]

/* Local variables for babtsp: */
struct LOC_babtsp {
  int n, inf;
  int *w;
  int *tweight;
  int *backptr, *best, *fwdptr;
} ;

/* Local variables for explore: */
struct LOC_explore {
  struct LOC_babtsp *LINK;
  int *row, *col;
  int size;
} ;

int min(i, j, LINK)
int i, j;
struct LOC_explore *LINK;
{
  if (i <= j)
    return i;
  else
    return j;
}

int reduce(row, col, rowred, colred, LINK)
int *row, *col, *rowred, *colred;
struct LOC_explore *LINK;
{
  int i, j;
  int rvalue = 0;
  int temp, FORLIM, FORLIM1;

  FORLIM = LINK->size;
  for (i = 0; i < FORLIM; i++) {   /* reduce rows */
    temp = LINK->LINK->inf;
    FORLIM1 = LINK->size;
    for (j = 0; j < FORLIM1; j++)
      temp = min(temp, LINK->LINK->w(row[i] - 1,col[j] - 1), LINK);
    if (temp > 0) {
      FORLIM1 = LINK->size;
      for (j = 0; j < FORLIM1; j++) {
	if (LINK->LINK->w(row[i] - 1,col[j] - 1) < LINK->LINK->inf)
	  LINK->LINK->w(row[i] - 1,col[j] - 1) -= temp;
      }
      rvalue += temp;
    }
    rowred[i] = temp;
  }  /* reduce columns  */
  FORLIM = LINK->size;
  for (j = 0; j < FORLIM; j++) {
    temp = LINK->LINK->inf;
    FORLIM1 = LINK->size;
    for (i = 0; i < FORLIM1; i++)
      temp = min(temp, LINK->LINK->w(row[i] - 1,col[j] - 1), LINK);
    if (temp > 0) {
      FORLIM1 = LINK->size;
      for (i = 0; i < FORLIM1; i++) {
	if (LINK->LINK->w(row[i] - 1,col[j] - 1) < LINK->LINK->inf)
	  LINK->LINK->w(row[i] - 1,col[j] - 1) -= temp;
      }
      rvalue += temp;
    }
    colred[j] = temp;
  }
  return rvalue;
}  /* reduce */

void bestedge(r, c, most, LINK)
int *r, *c, *most;
struct LOC_explore *LINK;
{
  int i, j, k, mincolelt, minrowelt, zeroes, FORLIM, FORLIM1, FORLIM2;

  *most = -LINK->LINK->inf;
  FORLIM = LINK->size;
  for (i = 1; i <= FORLIM; i++) {
    FORLIM1 = LINK->size;
    for (j = 1; j <= FORLIM1; j++) {
      if (LINK->LINK->w(LINK->row[i-1] - 1,LINK->col[j-1] - 1) == 0) {
	minrowelt = LINK->LINK->inf;
	zeroes = 0;
	FORLIM2 = LINK->size;
	for (k = 0; k < FORLIM2; k++) {
	  if (LINK->LINK->w(LINK->row[i-1] - 1,LINK->col[k] - 1) == 0)
	    zeroes++;
	  else
	    minrowelt = min(minrowelt, LINK->LINK->w(LINK->row[i-1] - 1,LINK->col[k] - 1), LINK);
	}
	if (zeroes > 1)
	  minrowelt = 0;
	mincolelt = LINK->LINK->inf;
	zeroes = 0;
	FORLIM2 = LINK->size;
	for (k = 0; k < FORLIM2; k++) {
	  if (LINK->LINK->w(LINK->row[k] - 1,LINK->col[j-1] - 1) == 0)
	    zeroes++;
	  else
	    mincolelt = min(mincolelt, LINK->LINK->w(LINK->row[k] - 1,LINK->col[j-1] - 1), LINK);
	}
	if (zeroes > 1)
	  mincolelt = 0;
	if (minrowelt + mincolelt > *most) {
	      /* a better edge has been found */
		*most = minrowelt + mincolelt;
	  *r = i;   /* row index of better edge */
	  *c = j;   /* column index of better edge */
	}
      }
    }
  }
}  /* bestedge */

void explore(edges, cost, row_, col_, LINK)
int edges, cost;
int *row_, *col_;
struct LOC_babtsp *LINK;
{
  struct LOC_explore V;
  int avoid, c, colrowval, first, i, j, last, lowerbound, most, r;
  int *colred, *newcol, *newrow, *rowred;
  int FORLIM, FORLIM1;

  colred=(int*)calloc(NN,sizeof(int));
  newcol=(int*)calloc(NN,sizeof(int));
  newrow=(int*)calloc(NN,sizeof(int));
  rowred=(int*)calloc(NN,sizeof(int));

  V.LINK = LINK;
  V.row = row_;
  V.col = col_;
  V.size = LINK->n - edges;   /* number of rows, cols left in matrix */
  cost += reduce(V.row, V.col, rowred, colred, &V);
  if (cost < *LINK->tweight) {
    if (edges == LINK->n - 2) {  /* last two edges are forced */
      FORLIM = LINK->n;
      for (i = 0; i < FORLIM; i++)
	LINK->best[i] = LINK->fwdptr[i];
      if (LINK->w(V.row[0] - 1,V.col[0] - 1) == LINK->inf)
	avoid = 1;
      else
	avoid = 2;
      LINK->best[V.row[0] - 1] = V.col[2 - avoid];
      LINK->best[V.row[1] - 1] = V.col[avoid-1];
      *LINK->tweight = cost;
    } else {   /* record chosen edge */
      bestedge(&r, &c, &most, &V);
      lowerbound = cost + most;
      LINK->fwdptr[V.row[r-1] - 1] = V.col[c-1];
      LINK->backptr[V.col[c-1] - 1] = V.row[r-1];
      last = V.col[c-1];   /* prevent cycles */
      while (LINK->fwdptr[last-1] != 0)
	last = LINK->fwdptr[last-1];
      first = V.row[r-1];
      while (LINK->backptr[first-1] != 0)
	first = LINK->backptr[first-1];
      colrowval = LINK->w(last-1,first-1);
      LINK->w(last-1,first-1) = LINK->inf;
      for (i = 0; i <= r - 2; i++)   /* remove row */
	newrow[i] = V.row[i];
      FORLIM = V.size;
      for (i = r; i < FORLIM; i++)
	newrow[i-1] = V.row[i];
      for (i = 0; i <= c - 2; i++)   /* remove col */
	newcol[i] = V.col[i];
      FORLIM = V.size;
      for (i = c; i < FORLIM; i++)
	newcol[i-1] = V.col[i];
      explore(edges + 1, cost, newrow, newcol, LINK);
      LINK->w(last-1,first-1) = colrowval;   /* restore previous values */
      LINK->backptr[V.col[c-1] - 1] = 0;
      LINK->fwdptr[V.row[r-1] - 1] = 0;
      if (lowerbound < *LINK->tweight) {
	LINK->w(V.row[r-1] - 1,V.col[c-1] - 1) = LINK->inf;
	    /* exclude edge alredy chosen */
	explore(edges, cost, V.row, V.col, LINK);
	LINK->w(V.row[r-1] - 1,V.col[c-1] - 1) = 0;
	    /* restore excluded edge */
      }
    }
  }
  FORLIM = V.size;
  for (i = 0; i < FORLIM; i++) {   /* unreduce matrix */
    FORLIM1 = V.size;
    for (j = 0; j < FORLIM1; j++)
      LINK->w(V.row[i] - 1,V.col[j] - 1) += rowred[i] + colred[j];
  }

  free(colred); free(newcol); free(newrow); free(rowred);
}  /* explore */


static void babtsp(n_, inf_, w_, route, tweight_)
int n_, inf_;
int *w_;
int *route;
int *tweight_;
{
  struct LOC_babtsp V;
  int *col, *row;
  int i;
  int index = 1;
  int FORLIM;
  col=(int*)calloc(NN,sizeof(int));
  row=(int*)calloc(NN,sizeof(int));

  V.n = n_;
  V.inf = inf_;
  V.w = w_;
  V.tweight = tweight_;
  FORLIM = V.n;
  for (i = 1; i <= FORLIM; i++) {
    row[i-1] = i;
    col[i-1] = i;
    V.fwdptr[i-1] = 0;
    V.backptr[i-1] = 0;
  }
  *V.tweight = V.inf;
  explore(0L, 0L, row, col, &V);
  FORLIM = V.n;
  for (i = 0; i < FORLIM; i++) {
    route[i] = index;
    index = V.best[index-1];
  }

  free(col); free(row);
}  /* babtsp */

#include "mex.h"

/**************************************************************************/
/* Main MEX-function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* MATLAB syntax: p=twoopt(W) */
{
  int n=mxGetM(prhs[0]);
  int *w=(int*)mxGetPr(prhs[0]);
  NN=n;
  
  plhs[0]=mxCreateNumericMatrix(1,n,mxINT32_CLASS,mxREAL);
  int *p=(int*)mxGetPr(plhs[0]);
  int weight=0;

  babtsp(n,1000000,w,p,&weight);
  printf("%d\n",weight);
}
