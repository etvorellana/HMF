#include <stdio.h>

#define  B0   0.675603595979828813
#define  B1  -0.175603595979828813
#define  D0   1.35120719195965763
#define  D1  -1.70241438391931525


/* function prototypes */
void force(long, long, long, double *, double *, double *, double *,double *, double *, float *, float *, double *, double *);
__global__ void step_type1(long , double, double, double *, double *, double *);
__global__ void step_type2(long , double, double *, double *, double *);


void timestep(long n, long nblock, long nthread, double dt, double *mx, double *my, double *magx, double *magy,
              double *magx_gpu, double *magy_gpu, double *r_gpu, double *p_gpu, double *f_gpu, float *sinr_gpu, float *cosr_gpu)
{
   double bb0,bb1,dd0,dd1;

   bb0=B0*dt;
   bb1=B1*dt;
   dd0=D0*dt;
   dd1=D1*dt;

   step_type1<<<nblock,nthread>>>(n,bb0,dd0,r_gpu,p_gpu,f_gpu);

   force(n,nblock,nthread,mx,my,magx,magy,r_gpu,f_gpu,sinr_gpu,cosr_gpu,magx_gpu,magy_gpu);

   step_type1<<<nblock,nthread>>>(n,bb1,dd1,r_gpu,p_gpu,f_gpu);

   force(n,nblock,nthread,mx,my,magx,magy,r_gpu,f_gpu,sinr_gpu,cosr_gpu,magx_gpu,magy_gpu);

   step_type1<<<nblock,nthread>>>(n,bb1,dd0,r_gpu,p_gpu,f_gpu);

   force(n,nblock,nthread,mx,my,magx,magy,r_gpu,f_gpu,sinr_gpu,cosr_gpu,magx_gpu,magy_gpu);

   step_type2<<<nblock,nthread>>>(n,bb0,r_gpu,p_gpu,f_gpu);

   return;
}
