
__global__ void force_aux0(long k, double *magx_gpu, double *magy_gpu)
{
     long i,cind;
     __shared__ double mx,my,mgx[512],mgy[512];

     cind=threadIdx.x;

     mgx[cind]=magx_gpu[cind];
     mgy[cind]=magy_gpu[cind];

     __syncthreads();

     mx=0.0;
     my=0.0;
     for (i=0;i<k;i++)
     {
          mx+=mgx[i];
          my+=mgy[i];
     };
     magx_gpu[0]=mx;
     magy_gpu[0]=my;

     return;
}
