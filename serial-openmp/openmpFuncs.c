#include <stdio.h>
#include <stdlib.h>	
#include <math.h>
#include <omp.h>

#include "openmpLibs.h"


//Coloquei todos os cabeçalhos e definições para a openLibs.h

void WaterBag(long n, long *idum, double p0, double r0, double *r, double *p)
{
	long i;
	double aux = .0;
	
	for (i = 0; i < n; i++)
	{
		r[i] = ((double)ran2(idum))*r0;
		p[i] = ((double)ran2(idum) - .5)*2.*p0;
		aux += p[i];
	}
	aux = aux / ((double)n);
		
	//#pragma omp parallel
	#pragma omp parallel for 
	//{
		//#pragma omp for private(i) 
		for (i = 0; i < n; i++)
		{
			p[i] -= aux;
		}
	//}
	return;
}

void KineticEnergy(long n, double *energKin, double *p)
{

	long i;
	*energKin = .0;
	double aux = .0;
	
	
	//Usar construtores compostos é mais eficiente 
	//#pragma omp parallel
	#pragma omp parallel for reduction(+:aux)
	//{
		//#pragma omp for private(i) reduction(+:aux)
		for (i = 0; i < n; i++)
		{
			aux += p[i] * p[i];
		}
	//}
	*energKin = aux / ((double)2 * n);
	return;
}

void PotentialEnergy(long n, double *energPot, double *r, double magX, double magY)
{
	//Vc não tem uma região paralela para usar aqui um single
	//#pragma omp single
	{
		*energPot = .0;
		*energPot = .5*(1. - magX*magX - magY*magY);
	}
	return;
}

void Force(long n, double *force, double *r, double *magX, double *magY)
{

	long i;
	double *as = (double *)malloc((double)n * sizeof(double));
	double *ac = (double *)malloc((double)n * sizeof(double));
	double aux1, aux2;
	double magX_, magY_;
	aux1 = .0;
	aux2 = .0;
	*magX = .0;
	*magY = .0;
	magX_ = *magX;
	magY_ = *magY;
// vou muadar algumas coisas aqui para ficar mais imples
//#pragma omp parallel
	#pragma omp parallel private(aux1, aux2)
	{	
		//#pragma omp for private(i, aux1, aux2) reduction(+:magX_, magY_)
		#pragma omp for reduction(+:magX_, magY_)
		for (i = 0; i < n; i++)
		{
			aux1 = sin(r[i]);
			aux2 = cos(r[i]);
			magX_ += aux2;
			magY_ += aux1;
			as[i] = aux1;
			ac[i] = aux2;
		}
	//} Não pode fechar aqui

		//*magX = magX_;
		//*magY = magY_;

		#pragma omp single
		{
			//*magX = *magX / ((double)n);
			*magX = magX_ / ((double)n);
			//*magY = *magY / ((double)n);
			*magY = magY_ / ((double)n);
		}

//#pragma omp parallel 
	//{
		//#pragma omp for private(i)
		#pragma omp for
		for (i = 0; i < n; i++)
		{
			force[i] = ac[i] * (*magY) - as[i] * (*magX);
		}
	}
	free(as);
	free(ac);

	return;
}

void Integration(long n, double dt, double *magX, double *magY, double *r, double *p, double *f)
{
	long i;
	double mx, my;

	mx = .0;
	my = .0;
	
//#pragma omp parallel
	#pragma omp parallel for
	//{
		//#pragma omp for private(i)
		for (i = 0; i<n; i++)
		{
			p[i] += B0*dt*f[i];
			r[i] += D0*dt*p[i];
		}
	//}
	
	Force(n, f, r, &mx, &my);
	
//#pragma omp parallel
	#pragma omp parallel for
	//{
		//#pragma omp for private(i)
		for (i = 0; i<n; i++)
		{
			p[i] += B1*dt*f[i];
			r[i] += D1*dt*p[i];
		}
	//}
	
	Force(n, f, r, &mx, &my);
	
//#pragma omp parallel
	#pragma omp parallel for
	//{
		//#pragma omp for private(i)
		for (i = 0; i<n; i++)
		{
			p[i] += B1*dt*f[i];
			r[i] += D0*dt*p[i];
		}
	//}		
	Force(n, f, r, &mx, &my);
	
//#pragma omp parallel
	#pragma omp parallel for
	//{
		//#pragma omp for private(i)
		for (i = 0; i<n; i++)
		{
			p[i] += B0*dt*f[i];
		}
	//}
	*magX = mx;
	*magY = my;


	return;
}


float ran2(long *idum)
{
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		idum2 = (*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ1;
			*idum = IA1*(*idum - k*IQ1) - k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ1;
	*idum = IA1*(*idum - k*IQ1) - k*IR1;
	if (*idum < 0) *idum += IM1;
	k = idum2 / IQ2;
	idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j = ((int)iy / NDIV);
	iy = iv[j] - idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp = AM*iy) > RNMX) return RNMX;
	else return temp;
}

