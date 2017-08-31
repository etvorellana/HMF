// Save the position and momentum of particle 0

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void save_seq( double time, long nseq, double *r_gpu, double *p_gpu, double *f_gpu, FILE *fseq, FILE *fseq2, FILE *fseq3)
{
   long i;
   double pp[nseq],rr[nseq],ff[nseq],tpi,r2;

   tpi=6.2831853071795864770;

   cudaMemcpy(pp, p_gpu, nseq*sizeof(double), cudaMemcpyDeviceToHost);
   cudaMemcpy(rr, r_gpu, nseq*sizeof(double), cudaMemcpyDeviceToHost);
   cudaMemcpy(ff, f_gpu, nseq*sizeof(double), cudaMemcpyDeviceToHost);

   for (i=0;i<nseq;i++)
   {
        r2=rr[i];
        while (r2>tpi) {r2-=tpi;};
        while (r2<0) {r2+=tpi;};
        fprintf(fseq,"%lf %lf %lf",time,r2,pp[i]);
        fprintf(fseq2,"%lf %lf ",time,r2);
        fprintf(fseq3,"%lf %lf ",time,ff[i]);
   };

   fprintf(fseq,"\n");
   fprintf(fseq2,"\n");
   fprintf(fseq3,"\n");

   return;
}

