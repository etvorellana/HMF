/* computes the force array */

#include <math.h>
#include <stdio.h>


__global__ void force_aux(long, double *, float *, float *, double *,double *);
__global__ void force_aux2(long, double, double, double *, double *, float *, float *);



void force(long n, long nblock, long nthread, double *mx, double *my, double *magx, double *magy, double *r_gpu,
           double *f_gpu, float *sinr_gpu, float *cosr_gpu, double *magx_gpu, double *magy_gpu)
{
   long i;

   force_aux<<<nblock,nthread>>>(n,r_gpu,sinr_gpu,cosr_gpu,magx_gpu,magy_gpu);

   cudaMemcpy(magx,magx_gpu,nblock*sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(magy,magy_gpu,nblock*sizeof(double),cudaMemcpyDeviceToHost);

   *mx=0.0;
   *my=0.0;
   for (i=0;i<nblock;i++)
   {
        *mx+=magx[i];
        *my+=magy[i];
   };

   *mx=*mx/((double) n);
   *my=*my/((double) n);


   force_aux2<<<nblock,nthread>>>(n,*mx,*my,r_gpu,f_gpu,sinr_gpu,cosr_gpu);


   return;
}
