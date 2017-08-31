
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"

__global__ void standard_aux( long n, double *var, double *var0, double *kinp)
{
    int cind,tid,i;
    double c;

    __shared__ double kk[MAXBLOCK];

    cind=threadIdx.x;
    tid=threadIdx.x+blockIdx.x*blockDim.x;

    kk[cind]=0.0;

    while(tid<n)
    {
          c=var[tid]-var0[tid];
          kk[cind]+=c*c;
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

