/* Computes kinetic energy */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"

/* FUnction prototype */

__global__ void moms_aux( long, long, double *, double *);


/* Function */

void moments(double *mom, long n, long order, long nblock, long nthread, double *p, double *p_gpu, double *kinp_gpu, double *kinp)
{
   long i;

   moms_aux<<<nblock,nthread>>>( n, order, p_gpu, kinp_gpu);

   cudaMemcpy(kinp, kinp_gpu, nblock*sizeof(double), cudaMemcpyDeviceToHost);

   *mom=0.0;
   for (i=0;i<nblock;i++)
   {
       *mom+=kinp[i];
   };
   *mom=*mom/(2.0*((double) n));

   return;
}
