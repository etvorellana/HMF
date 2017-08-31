#include <stdlib.h>
#include <math.h>
#include <stdio.h>


/* Function prototypes */

float ran2(long *);


void condini(long n, long *idum, double p0, double theta0, double r[], double p[])
{

   long i;
   double lt,ptot;
 
   lt=0.0; 
   for (i=0;i<n;i++)
       {
       r[i]=((double) ran2(idum))*theta0;
       p[i]=(((double) ran2(idum))-0.5)*p0;
       lt+=p[i];
       };
   lt=lt/((double) n);

   ptot=0.0;
   for (i=0;i<n;i++)
       {
       p[i]=p[i]-lt;
       ptot+=p[i];     
       };

   printf("Total initial momentum: %lf\n",ptot);

   return;
}
