/* Computes the standard deviation of traveled distance and variation of momentum */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"


/* Function prototype */

__global__ void standard_aux( long n, double *var, double *var0, double *kinp);


/* Function */

void standard_dev( long n, long nblock, long nthread, double *r_gpu, double *p_gpu, double *r0_gpu, double *p0_gpu, double *sr, double *sp, double *kinp, double *kinp_gpu)
{
   long i;

   standard_aux<<<nblock,nthread>>>( n, r_gpu, r0_gpu, kinp_gpu);

   cudaMemcpy(kinp, kinp_gpu, nblock*sizeof(double), cudaMemcpyDeviceToHost);

   *sr=0.0;
   for (i=0;i<nblock;i++)
   {
       *sr+=kinp[i];
   };
   *sr/=(double) n;
   *sr=sqrt(*sr);

   standard_aux<<<nblock,nthread>>>( n, p_gpu, p0_gpu, kinp_gpu);

   cudaMemcpy(kinp, kinp_gpu, nblock*sizeof(double), cudaMemcpyDeviceToHost);

   *sp=0.0;
   for (i=0;i<nblock;i++)
   {
       *sp+=kinp[i];
   };
   *sp/=(double) n;
   *sp=sqrt(*sp);

   return;
}
