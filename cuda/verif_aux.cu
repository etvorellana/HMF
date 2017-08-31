/* Auxiliary routine to compute kinetic energy */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"

__global__ void verif_aux( long n, double *var, double *mom4)
{

    int cind,tid,i;

    __shared__ double a1[MAXBLOCK];

    cind=threadIdx.x;
    tid=threadIdx.x+blockIdx.x*blockDim.x;

    a1[cind]=var[tid];

    __syncthreads();


    i=blockDim.x/2;
    while(i!=0)
    {
         if(cind<i)
         {
            a1[cind]+=a1[cind+i];
         };

         __syncthreads();

         i/=2;

   };

   if(cind==0)
   {
      mom4[blockIdx.x]=a1[0];
   };

    return;

}
