// Initializes the random number generators in the GPU

#include <stdlib.h>
#include <math.h>
#include <stdio.h>


__global__ void initialize_ran( long *primes)
{
     long idum0;
     int i0,j0;
     float rn0;

     i0=threadIdx.x;
     j0=blockIdx.x;

     idum0=-primes[blockDim.x*j0+i0];

#include "ran2.inc"

     primes[blockDim.x*j0+i0]=idum0;

     rn0=rn0;

     return;
}
