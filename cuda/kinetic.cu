/* Computes kinetic energy */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"

/* FUnction prototype */

__global__ void kin_aux( long, double *, double *, double *, double *);


/* Function */

void kinetic(double *kin, double *m4, double *m6, long n, long nblock, long nthread, double *p_gpu,
             double *kinp_gpu, double *mom4_gpu, double *mom6_gpu, double *mom4, double *mom6, double *kinp)
{
   long i;
   double kk;

   kin_aux<<<nblock,nthread>>>( n, p_gpu, kinp_gpu, mom4_gpu, mom6_gpu);

   cudaMemcpy(kinp, kinp_gpu, nblock*sizeof(double), cudaMemcpyDeviceToHost);
   cudaMemcpy(mom4, mom4_gpu, nblock*sizeof(double), cudaMemcpyDeviceToHost);
   cudaMemcpy(mom6, mom6_gpu, nblock*sizeof(double), cudaMemcpyDeviceToHost);

   *kin=0.0;
   *m4=0.0;
   *m6=0.0;
   for (i=0;i<nblock;i++)
   {
       *kin+=kinp[i];
       *m4+=mom4[i];
       *m6+=mom6[i];
   };
   *kin=*kin/(2.0*((double) n));
   *m4=*m4/((double) n);
   *m6=*m6/((double) n);

   kk=2.0*(*kin);

   *m4/=kk*kk;
   *m6/=kk*kk*kk;

   return;
}
