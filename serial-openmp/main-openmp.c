/*
 * 
 * main.c
 * 
 * Copyright 2017 Esbel Tomas Valero Orellana <evalero@margot>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <stdio.h>
#include <stdlib.h>	
#include <math.h>
#include <omp.h>

#define	pi	3.14159265359
#define dpi	6.28318530718
#define B0   0.675603595979828813
#define B1  -0.175603595979828813
#define D0   1.35120719195965763
#define D1  -1.70241438391931525

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

void Integration(long n, double dt, double *magX, double *magY, double *r, double *p, double *f);
void WaterBag(long n, long *idum, double p0, double r0, double *r, double *p);
void KineticEnergy(long n, double *energKin, double *p);
void PotentialEnergy(long n, double *energPot, double *r, double magX, double magY);
void Force(long n, double *force, double *r, double *magX, double *magY);
float ran2(long *idum);

int main(int argc, char **argv)
{
	long n, seed, idum;
	double p0, r0;
	double energKin, energPot, magX, magY, energ, energ0, error;
	double time, finalTime, timeStep, timeOutput, timeCount;
	

	FILE *init = fopen("./initialPhase.dat", "w");
	FILE *enrg = fopen("./energy.dat", "w");
	FILE *fmag = fopen("./magnet.dat", "w");
	FILE *finalSpace = fopen("./finalPhase.dat", "w");


	FILE *in = fopen("./input.in", "r");
	fscanf(in, "%ld", &n);
	fscanf(in, "%lf", &finalTime);
	fscanf(in, "%lf", &timeStep);
	fscanf(in, "%lf", &timeOutput);
	fscanf(in, "%lf", &p0);
	fscanf(in, "%lf", &r0);
	fscanf(in, "%ld", &seed);
	fclose(in);

	idum = -seed;

	/* Inicialização de variáveis*/
	energKin = ran2(&idum);
	energPot = ran2(&idum);
	magX = ran2(&idum);
	magY = ran2(&idum);

	double *r = (double *)malloc((double)n * sizeof(double));
	double *p = (double *)malloc((double)n * sizeof(double));
	double *force = (double *)malloc((double)n * sizeof(double));
	
	WaterBag(n, &idum, p0, r0, r, p);

	#pragma omp parallel for
		for (long i = 0; i < n; i++)
		{
			fprintf(init, "%lf\t%lf\n", r[i], p[i]);
		}

	KineticEnergy(n, &energKin, p);
	Force(n, force, r, &magX, &magY);
	PotentialEnergy(n, &energPot, r, magX, magY);
	energ0 = energKin + energPot;

	//cout << "Energia Cinetica Inicial: " << energKin << endl;
	//cout << "Energia Potencial Inicial: " << energPot << endl;
	//cout << "Energia Total Inicial: " << energ0 << endl;
	//cout << "Magnetizacoes iniciais:   MagX: " << magX << "  MagY: " << magY << endl;

	error = .0;
	time = .0;
	timeCount = .0;

	while (time < finalTime)
	{
			
		Integration(n, timeStep, &magX, &magY, r, p, force);

		time += timeStep;
		timeCount += timeStep;

		if (timeCount >= timeOutput)
		{
			KineticEnergy(n, &energKin, p);
			PotentialEnergy(n, &energPot, r, magX, magY);
			energ = energKin + energPot;
			error = (energ - energ0) / energ0;
			error = fabs(error);
			// Colocar aqui um if para parar a simulação quandoo errofor grande 
			// Definir erro limite aceitavel
			//printf("%lf\t%1.2le\t%lf\t%lf\t%lf\n", time, error, npMean, var, statMoment4);
			timeCount = 0.0;			
			fprintf(enrg, "%lf\t%lf\t%lf\n", time, energKin, energPot);
			fprintf(fmag, "%lf %lf %lf %lf\n", time, magX, magY, sqrt(magX*magX + magY*magY));
		}
	}

	printf("Salvando os espacos de fase finais\n");

	double rr = .0;
    for (long i = 0; i < n; i++)
		{
			rr = r[i];
			while (rr > dpi / 2.)
			{
				rr -= dpi;
			}
			while (rr < -dpi / 2.)
			{
			rr += dpi;
			}
		fprintf(finalSpace, "%lf\t%lf\n", r[i], p[i]);
		}
 
	free(r);
	free(p);
	free(force);

	fclose(enrg);
	fclose(finalSpace);
	fclose(fmag);
	fclose(init);

	return 0;
}

void WaterBag(long n, long *idum, double p0, double r0, double *r, double *p)
{

	double aux = .0;
	#pragma omp parallel for reduction(+:aux)
		for (long i = 0; i < n; i++)
		{
			r[i] = ((double)ran2(idum))*r0;
			p[i] = ((double)ran2(idum) - .5)*2.*p0;
			aux += p[i];
		aux = aux / ((double)n);
	}
	#pragma omp parallel for
		for (long i = 0; i < n; i++)
		{
			p[i] -= aux;
		}
	return;
}

void KineticEnergy(long n, double *energKin, double *p)
{
	*energKin = .0;
	double aux = .0;
	#pragma omp parallel for reduction(+:aux)
		for (long i = 0; i < n; i++)
		{
			aux += p[i] * p[i];
		}
	*energKin = aux / ((double)2 * n);
	return;
}

void PotentialEnergy(long n, double *energPot, double *r, double magX, double magY)
{
	#pragma omp single
	{
		*energPot = .0;
		*energPot = .5*(1. - magX*magX - magY*magY);
	}
	return;
}

void Force(long n, double *force, double *r, double *magX, double *magY)
{
	double *as = (double *)malloc((double)n * sizeof(double));
	double *ac = (double *)malloc((double)n * sizeof(double));
	double aux1, aux2;

	aux1 = .0;
	aux2 = .0;
	*magX = .0;
	*magY = .0;

	#pragma omp parallel for reduction(+:magX,magY)
		for (long i = 0; i < n; i++)
		{
			aux1 = sin(r[i]);
			aux2 = cos(r[i]);
			*magX += aux2;
			*magY += aux1;
			as[i] = aux1;
			ac[i] = aux2;
		}
	#pragma omp single
	{
		*magX = *magX / ((double)n);
		*magY = *magY / ((double)n);
	}

	#pragma omp parallel for
		for (long i = 0; i < n; i++)
		{
			force[i] = ac[i] * (*magY) - as[i] * (*magX);
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
	
	#pragma omp parallel for
		for (i = 0; i<n; i++)
		{
			p[i] += B0*dt*f[i];
			r[i] += D0*dt*p[i];
		}
	
	Force(n, f, r, &mx, &my);
	#pragma omp parallel for
		for (i = 0; i<n; i++)
		{
			p[i] += B1*dt*f[i];
			r[i] += D1*dt*p[i];
		}
	Force(n, f, r, &mx, &my);
	#pragma omp parallel for
		for (i = 0; i<n; i++)
		{
			p[i] += B1*dt*f[i];
			r[i] += D0*dt*p[i];
		}
	Force(n, f, r, &mx, &my);
	#pragma omp parallel for
		for (i = 0; i<n; i++)
		{
			p[i] += B0*dt*f[i];
		}
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

