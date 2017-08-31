/* Auxiliary routine to compute kinetic energy */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"

__global__ void moms_aux( long n, long order, double *p, double *kinp)
{

    int cind,tid,i;

    __shared__ double kk[MAXBLOCK];

    cind=threadIdx.x;
    tid=threadIdx.x+blockIdx.x*blockDim.x;

    kk[cind]=0.0;

    while(tid<n)
    {
          kk[cind]=p[tid];
          for(i=0;i<order-1;i++)
          {kk[cind]=kk[cind]*p[tid];};
          tid+=blockDim.x*gridDim.x;
    };

    __syncthreads();


    i=blockDim.x/2;
    while(i!=0)
    {
         if(cind<i)
         {
            kk[cind]+=kk[cind+i];
         };

         __syncthreads();

         i/=2;

   };

   if(cind==0)
   {
      kinp[blockIdx.x]=kk[0];
   };

    return;

}
