#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
//#include <stats.h>

extern void hcod_coef(double *xM, double *yM, int *nrow, int *ncol, int *hx, int *hy, double *rhoh) 
{
  double sx = 0.0, sy = 0.0, cross = 0.0;
  int i,j;
  double x,y,xh,yh;
  for (i = 0; i < *nrow - *hy; i++ )
  { 
    if (*hx >= 0)
    {
      for (j = 0; j < *ncol - *hx; j++)
      {
        x = xM[i + *nrow*j];
        y = yM[i + *nrow*j];
        xh = xM[i + *hy + (j + *hx)*(*nrow)];
        yh = yM[i + *hy + (j + *hx)*(*nrow)];
        cross += (x-xh)*(y-yh); 
        sx += (x-xh)*(x-xh); 
        sy += (y-yh)*(y-yh);
      }
    }
    else
    { 
      for (j = -*hx; j < *ncol; j++)
      {
        x = xM[i + *nrow*j];
        y = yM[i + *nrow*j];
        xh = xM[i + *hy + (j + *hx)*(*nrow)];
        yh = yM[i + *hy + (j + *hx)*(*nrow)];
        cross += (x-xh)*(y-yh);
        sx += (x-xh)*(x-xh);
        sy += (y-yh)*(y-yh);
      }
    } 
  }
  *rhoh = cross/sqrt(sx*sy);
}
