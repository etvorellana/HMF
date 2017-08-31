/* Computes kinetic energy */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "precision.inc"

/* FUnction prototype */

__global__ void verif_aux( long, double *, double *);


/* Function */

void verif(double *kurt_p, double *kurt_r, long n, long nblock, long nthread, double *p_gpu, double *r_gpu,
           double *kurt_gpu, double *kurt)
{
   long i;
   double m2,m4,avg,a;


   verif_aux<<<nblock,nthread>>>( n, p_gpu, kurt_gpu);

   cudaMemcpy(kurt, kurt_gpu, nblock*sizeof(double), cudaMemcpyDeviceToHost);

   m2=0.0;
   m4=0.0;
   avg=0.0;

   for (i=0;i<nblock;i++)
   {
        avg+=kurt[i];
   };
   avg/=(double) nblock;

   for (i=0;i<nblock;i++)
   {
       a=kurt[i]-avg;
       a*=a;
       m2+=a;
       m4+=a*a;
   };
   m2/=(double) nblock;
   m4/=(double) nblock;

   *kurt_p=m4/(m2*m2);



   verif_aux<<<nblock,nthread>>>( n, r_gpu, kurt_gpu);

   cudaMemcpy(kurt, kurt_gpu, nblock*sizeof(double), cudaMemcpyDeviceToHost);

   m2=0.0;
   m4=0.0;
   avg=0.0;

   for (i=0;i<nblock;i++)
   {
        avg+=kurt[i];
   };
   avg/=(double) nblock;

   for (i=0;i<nblock;i++)
   {
       a=kurt[i]-avg;
       a*=a;
       m2+=a;
       m4+=a*a;
   };
   m2/=(double) nblock;
   m4/=(double) nblock;

   *kurt_r=m4/(m2*m2);


   return;
}
