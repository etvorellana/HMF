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

#include "openmpLibs.h"

int main(int argc, char **argv)
{
	double t_0, t_F, dT;
	t_0 = omp_get_wtime();
	
	long n, seed, idum;
	double p0, r0;
	double energKin, energPot, magX, magY, energ, energ0, error;
	double time, finalTime, timeStep, timeOutput, timeCount;
	

	FILE *init = fopen("./initialPhaseOMP.dat", "w");
	FILE *enrg = fopen("./energyOMP.dat", "w");
	FILE *fmag = fopen("./magnetOMP.dat", "w");
	FILE *finalSpace = fopen("./finalPhaseOMP.dat", "w");

	FILE *in = fopen("./input.in", "r");
	fscanf(in, "%ld", &n);
	fscanf(in, "%lf", &finalTime);
	fscanf(in, "%lf", &timeStep);
	fscanf(in, "%lf", &timeOutput);
	fscanf(in, "%lf", &p0);
	fscanf(in, "%lf", &r0);
	fscanf(in, "%ld", &seed);
	fclose(in);

	printf("3\n");	
	
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
	//A escrita em disco não adianta ser paralelizada, pode dar m...
	//#pragma omp parallel for
	for (long i = 0; i < n; i++)
	{
		fprintf(init, "%lf\t%lf\n", r[i], p[i]);
	}
	
	KineticEnergy(n, &energKin, p);
	Force(n, force, r, &magX, &magY);
	PotentialEnergy(n, &energPot, r, magX, magY);
	energ0 = energKin + energPot;
	
	
	printf("Energia Cinetica Inicial: %lf \n", energKin);
	printf("Energia Potencial Inicial: %lf \n", energPot);
	printf("Energia Total Inicial: %lf \n", energ0);
	printf("Magnetizacoes iniciais:   MagX: %lf  MagY: %lf\n", magX, magY);
	

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
			if (error > RMAX){
				printf("%lf\t%1.2le\n", time, error);
			}
			
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
	t_F = omp_get_wtime();
	printf("d_T = %lf\n", t_F - t_0);
	
	return 0;
}
