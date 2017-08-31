/* Auxiliary routine to compute kinetic energy */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"

__global__ void kin_aux( long n, double *p, double *kinp, double *mom4, double *mom6)
{

    int cind,tid,i;

    __shared__ double kk[MAXBLOCK],mm4[MAXBLOCK],mm6[MAXBLOCK],a1[MAXBLOCK];

    cind=threadIdx.x;
    tid=threadIdx.x+blockIdx.x*blockDim.x;

    kk[cind]=0.0;
    mm4[cind]=0.0;
    mm6[cind]=0.0;

    while(tid<n)
    {
          a1[cind]=p[tid]*p[tid];
          kk[cind]+=a1[cind];
          mm4[cind]+=a1[cind]*a1[cind];
          mm6[cind]+=a1[cind]*a1[cind]*a1[cind];
          tid+=blockDim.x*gridDim.x;
    };

    __syncthreads();


    i=blockDim.x/2;
    while(i!=0)
    {
         if(cind<i)
         {
            kk[cind]+=kk[cind+i];
            mm4[cind]+=mm4[cind+i];
            mm6[cind]+=mm6[cind+i];
         };

         __syncthreads();

         i/=2;

   };

   if(cind==0)
   {
      kinp[blockIdx.x]=kk[0];
      mom4[blockIdx.x]=mm4[0];
      mom6[blockIdx.x]=mm6[0];
   };

    return;

}
