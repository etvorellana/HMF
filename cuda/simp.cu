#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"


/* Function prototypes */

void condini(long, long *, double, double, double [], double[]);
float ran2(long *);
void kinetic(double *, double *, double *, long, long, long, double *, double *, double *, double *, double *, double *, double *);
void timestep(long, long, long, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, float *, float *);
void force(long, long, long, double *, double *, double *, double *, double *, double *, float *, float *, double *, double *);
void verif(double *, double *, long, long, long, double *, double *, double *, double *);
void save_seq(double, long, double *, double *, double *, FILE *,FILE *, FILE *);
void standard_dev(long, long, long, double *, double *, double *, double *, double *, double *, double *, double *);
void correlations( long, long, long, double *, double *, double *, double *, double *, double *, double *, double *);
__global__ void compute_sojourn( long, double, double, double, double *, double *, double *, long *, long *, long *, double *, long *, long *, long *);

int main()
{

   long n,i,k,seed,idum,nblock,nthread,nblock2,nthread2,ifsave,ifverif,ifsojournp,ifsojournr,ifseq,ifstnd,ifsavesoj,nseq,
        device_num,nsojourn,numbin,soj_count;
   double dt,tf,ts,tseq,energ,energ0,p0,theta0,kin,potential,time,tcount,err,mx,my,mtot,m4,m6,
          tpi,r2,kurt_p,kurt_r,tseq_count,cf,cp,sr,sp,dp,dtheta,binini,binfin,deltabin,sumhist,soj_avg,soj_sig,
          sojmin,sojmax,binini_log,binfin_log,deltabin_log;
   FILE *arq1,*enrg,*magnet,*maginst,*fmoms,*phase0,*phase,*fverif,*fseq,*fseq2,*fseq3,*ferr,*difus,*correl,*fsoujp,*fsoujr,*fhistp,*fhistp_log,*fssig;

   double *r_gpu,*p_gpu,*f_gpu,*magx_gpu,*magy_gpu,*kinp_gpu,*mom4_gpu,*mom6_gpu,*kurt_gpu,*r0_gpu,*p0_gpu,*f0_gpu,
          *sojournp_gpu,*sojournr_gpu;
   float *sinr_gpu,*cosr_gpu;
   long *kp_gpu,*kpold_gpu,*ntp_gpu,*kr_gpu,*krold_gpu,*ntr_gpu;

/* Reads input parameters */

   arq1 = fopen("data.in","r");
   fscanf(arq1,"%lf",&dt);
   fscanf(arq1,"%lf",&tf);
   fscanf(arq1,"%lf",&ts);
   fscanf(arq1,"%lf",&tseq);
   fscanf(arq1,"%lf",&p0);
   fscanf(arq1,"%lf",&theta0);
   fscanf(arq1,"%li",&n);
   fscanf(arq1,"%li",&seed);
   fscanf(arq1,"%li",&ifsave);
   fscanf(arq1,"%li",&ifverif);
   fscanf(arq1,"%li",&ifsojournp);
   fscanf(arq1,"%li",&ifsojournr);
   fscanf(arq1,"%li",&ifsavesoj);
   fscanf(arq1,"%li",&ifseq);
   fscanf(arq1,"%li",&ifstnd);
   fscanf(arq1,"%li",&nsojourn);
   fscanf(arq1,"%li",&nseq);
   fscanf(arq1,"%lf",&dp);
   fscanf(arq1,"%lf",&dtheta);
   fscanf(arq1,"%li",&numbin);
   fscanf(arq1,"%lf",&binini);
   fscanf(arq1,"%lf",&binfin);
   fscanf(arq1,"%lf",&sojmin);
   fscanf(arq1,"%lf",&sojmax);
   fscanf(arq1,"%li",&nblock);
   fscanf(arq1,"%li",&nthread);
   fscanf(arq1,"%li",&nblock2);
   fscanf(arq1,"%li",&nthread2);
   fscanf(arq1,"%li",&device_num);
   fclose(arq1);

   cudaSetDevice(device_num);

   printf("=============================================\n");
   printf("  dt,tf,ts= %lf  %lf  %lf\n",dt,tf,ts);
   printf("  p0,theta0=  %lf %lf\n",p0,theta0);
   printf("  n, ifsave, ifverif, ifseq= %li  %li %li %li\n",n,ifsave,ifverif,ifseq);
   printf("  nblock, nthread, nblock2, nthread2 = %li  %li %li %li\n",nblock,nthread,nblock2,nthread2);
   printf("  device_num: %li\n",device_num);
   printf("=============================================\n\n");




/* Allocating memory on the CPU */

   double *r=(double*)malloc(n*sizeof(double));
   double *p=(double*)malloc(n*sizeof(double));
   double *f=(double*)malloc(n*sizeof(double));
   double *magx=(double*)malloc(nblock*sizeof(double));
   double *magy=(double*)malloc(nblock*sizeof(double));
   double *kinp=(double*)malloc(nblock2*sizeof(double));
   double *mom4=(double*)malloc(nblock2*sizeof(double));
   double *mom6=(double*)malloc(nblock2*sizeof(double));
   double *kurt=(double*)malloc(nblock2*sizeof(double));
   long *kold=(long*)malloc(nsojourn*sizeof(long));


/* Alocating memory on the GPU */

   cudaMalloc( (void**) &r_gpu, n*sizeof(double) );
   cudaMalloc( (void**) &p_gpu, n*sizeof(double) );
   cudaMalloc( (void**) &f_gpu, n*sizeof(double) );
   cudaMalloc( (void**) &sinr_gpu, n*sizeof(float) );
   cudaMalloc( (void**) &cosr_gpu, n*sizeof(float) );
   cudaMalloc( (void**) &magx_gpu, nblock*sizeof(double) );
   cudaMalloc( (void**) &magy_gpu, nblock*sizeof(double) );
   cudaMalloc( (void**) &kinp_gpu, nblock2*sizeof(double) );
   cudaMalloc( (void**) &mom4_gpu, nblock2*sizeof(double) );
   cudaMalloc( (void**) &mom6_gpu, nblock2*sizeof(double) );
   cudaMalloc( (void**) &kurt_gpu, nblock2*sizeof(double) );
   cudaMalloc( (void**) &kp_gpu, nsojourn*sizeof(long) );
   cudaMalloc( (void**) &kpold_gpu, nsojourn*sizeof(long) );
   cudaMalloc( (void**) &ntp_gpu, nsojourn*sizeof(long) );
   cudaMalloc( (void**) &kr_gpu, nsojourn*sizeof(long) );
   cudaMalloc( (void**) &krold_gpu, nsojourn*sizeof(long) );
   cudaMalloc( (void**) &ntr_gpu, nsojourn*sizeof(long) );

   if (ifstnd==1)
   {
       cudaMalloc( (void**) &r0_gpu, n*sizeof(double) );
       cudaMalloc( (void**) &p0_gpu, n*sizeof(double) );
       cudaMalloc( (void**) &f0_gpu, n*sizeof(double) );
   };


/* Initializing the random number generator */

   idum=-seed;
   kin=(float) ran2(&idum);


/* Prepares initial conditions and copy the arrays to the GPU */

   condini(n,&idum,p0,theta0,r,p);

   cudaMemcpy(r_gpu, r, n*sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(p_gpu, p, n*sizeof(double), cudaMemcpyHostToDevice);




/* Saves the initial distribution if ifsave=1 */

   if (ifsave==1)
   {
       phase0=fopen("initial.dat","w");
       for (i=0;i<n;i++)
       {
            fprintf(phase0,"%le %le\n",r[i],p[i]);
       };
       fclose(phase0);
   };


/* Compute initial kinetic and potential energies */

   force(n,nblock,nthread,&mx,&my,magx,magy,r_gpu,f_gpu,sinr_gpu,cosr_gpu,magx_gpu,magy_gpu);

   if (ifstnd==1)
   {
       cudaMemcpy(r0_gpu, r, n*sizeof(double), cudaMemcpyHostToDevice);
       cudaMemcpy(p0_gpu, p, n*sizeof(double), cudaMemcpyHostToDevice);
       cudaMemcpy(f0_gpu, f_gpu, n*sizeof(double), cudaMemcpyDeviceToDevice);
   }

   kinetic( &kin, &m4, &m6, n, nblock2, nthread2, p_gpu, kinp_gpu, mom4_gpu, mom6_gpu, mom4, mom6, kinp);
   potential=(1.0-mx*mx-my*my)/2.0;
   mtot=sqrt(mx*mx+my*my);

   printf("Initial kinetic and potential energies per particle:\n\n");
   printf("K= %lf  V= %lf\n\n",kin,potential);
   printf("Initial magnetization components:  %lf  %lf\n",mx,my);

   energ0=kin+potential;


   enrg=fopen("energy.dat","w");
   magnet=fopen("magnet.dat","w");
   maginst=fopen("maginst.dat","w");
   fmoms=fopen("moments.dat","w");
   ferr=fopen("error.dat","w");


   time=0.0;
   fprintf(enrg,"%lf  %lf  %lf\n",time,kin,potential);
   fprintf(magnet,"%lf  %lf\n",time,mtot);
   fprintf(maginst,"%lf  %lf  %lf\n",time,mx,my);
   fprintf(fmoms,"%lf  %lf  %lf\n",time,m4,m6);

   tcount=0.0;
   tseq_count=0.0;

/* time loop */

   while(time<=tf)
   {
      timestep(n,nblock,nthread,dt,&mx,&my,magx,magy,magx_gpu,magy_gpu,r_gpu,p_gpu,f_gpu,sinr_gpu,cosr_gpu);

      time+=dt;
      tcount+=dt;
      tseq_count+=dt;

      if (tcount>=ts)
      {
         kinetic( &kin, &m4, &m6, n, nblock2, nthread2, p_gpu, kinp_gpu, mom4_gpu, mom6_gpu, mom4, mom6, kinp);
         potential=(1.0-mx*mx-my*my)/2.0;
         energ=kin+potential;
         mtot=sqrt(mx*mx+my*my);
         err=(energ-energ0)/energ0;
         err=fabs(err);
         printf("  %lf  %e  %f\n",time,err,m4);
         tcount=0.0;
         fprintf(enrg,"%lf  %lf  %lf\n",time,kin,potential);
         fprintf(magnet,"%lf  %lf\n",time,mtot);
         fprintf(maginst,"%lf  %lf  %lf\n",time,mx,my);
         fprintf(fmoms,"%lf  %lf  %lf\n",time,m4,m6);
         fprintf(ferr,"%lf  %le\n",time,err);

      };
 
   };

   fclose(enrg);
   fclose(magnet);
   fclose(maginst);
   fclose(fmoms);
   fclose(fverif);

   if (ifsave==1)
   {
       tpi=6.2831853071795864770;
       cudaMemcpy(r, r_gpu, n*sizeof(double), cudaMemcpyDeviceToHost);
       cudaMemcpy(p, p_gpu, n*sizeof(double), cudaMemcpyDeviceToHost);
       phase=fopen("phase.dat","w");
       for (i=0;i<n;i++)
       {
            r2=r[i];
            while (r2>tpi) {r2-=tpi;};
            while (r2<0) {r2+=tpi;};
            fprintf(phase,"%le %le\n",r2,p[i]);
       };
       fclose(phase);
   };

   free(magx);
   free(magy);
   free(r);
   free(p);
   free(f);
   free(mom4);
   free(mom6);
   free(hist_sojp);

   cudaFree(magx_gpu);
   cudaFree(magy_gpu);
   cudaFree(r_gpu);
   cudaFree(p_gpu);
   cudaFree(f_gpu);
   cudaFree(mom4_gpu);
   cudaFree(mom6_gpu);

   return 0;

}
