#define w(i,j) w[(i)+n*(j)]

void twoopt(int n, int *w, int *route, int *tweight, int *tmp)
{
  int ahead, i, i1, i2;
  int index = 1;
  int j, j1, j2, last, limit, max, max1, next, s1, s2, t1, t2;
  int *ptr=tmp;

  for (i = 1; i < n; i++)
    ptr[route[i-1] - 1] = route[i];
  ptr[route[n-1] - 1] = route[0];
  do {
    max = 0;
    i1 = 1;
    for (i = 1; i <= n - 2; i++) {
      if (i == 1)
	limit = n - 1;
      else
	limit = n;
      i2 = ptr[i1-1];
      j1 = ptr[i2-1];
      for (j = i + 2; j <= limit; j++) {
	j2 = ptr[j1-1];
	max1 = w(i1-1,i2-1) + w(j1-1,j2-1) - w(i1-1,j1-1) - w(i2-1,j2-1);
	if (max1 > max) {  /* better pair of edges has been found */
	  s1 = i1;
	  s2 = i2;
	  t1 = j1;
	  t2 = j2;
	  max = max1;
	}
	j1 = j2;
      }
      i1 = i2;
    }
    if (max > 0) {  /* swap pair of edges */
      ptr[s1-1] = t1;
      next = s2;
      last = t2;
      do {   /* reverse appriopriate links */
	ahead = ptr[next-1];
	ptr[next-1] = last;
	last = next;
	next = ahead;
      } while (next != t2);
      *tweight -= max;   /* route is now shorter */
    }
  } while (max != 0);
  for (i = 0; i < n; i++) {
    route[i] = index;
    index = ptr[index-1];
  }
}  /* twoopt  */

#include "mex.h"

/**************************************************************************/
/* Main MEX-function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* MATLAB syntax: p=twoopt(W,p0) */
{
  int i;

  int n=mxGetM(prhs[0]);
  int *w=(int*)mxGetPr(prhs[0]);
  
  plhs[0]=mxDuplicateArray(prhs[1]);
  int *p=(int*)mxGetPr(plhs[0]);

  int weight;
  int *tmp=mxCalloc(n,sizeof(int));

  weight=w(n-1,0);
  for (i=1; i<n; i++)
    weight+=w(i-1,i);

  twoopt(n,w,p,&weight,tmp);

  //printf("%d\n",weight);

  mxFree(tmp);
}
