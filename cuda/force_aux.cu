#include "precision.inc"


__global__ void force_aux(long n, double *r_gpu, float *sinr_gpu, float *cosr_gpu, double *magx_gpu, double *magy_gpu)
{
    long i,tid,cind;

    __shared__ float mmx[NDIM0],mmy[NDIM0];

    cind=threadIdx.x;
    tid=threadIdx.x+blockIdx.x*blockDim.x;

    mmx[cind]=0.0;
    mmy[cind]=0.0;

    while (tid<n)
    {

        cosr_gpu[tid]=cos(r_gpu[tid]);
        sinr_gpu[tid]=sin(r_gpu[tid]);

        mmx[cind]+=cosr_gpu[tid];
        mmy[cind]+=sinr_gpu[tid];

        tid+=blockDim.x*gridDim.x;
    };

    __syncthreads();

    i=blockDim.x/2;
    while (i!=0)
    {
          if (cind<i)
          {
             mmx[cind]+=mmx[cind+i];
             mmy[cind]+=mmy[cind+i];
          };

          __syncthreads();

          i/=2;

    };

    if (cind==0)
    {
       magx_gpu[blockIdx.x]=mmx[0];
       magy_gpu[blockIdx.x]=mmy[0];
    };

    return;

}
