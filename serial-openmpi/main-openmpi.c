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
#include <mpi.h>
#include <omp.h>

#define pi	3.14159265359
#define dpi	6.28318530718
#define  B0   0.675603595979828813
#define  B1  -0.175603595979828813
#define  D0   1.35120719195965763
#define  D1  -1.70241438391931525

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

#define RMAX 1.0E-5

//void Integration(long n, double dt, double *magX, double *magY, double *r, double *p, double *f);
void Integration(long n, long nLoc, double dt, double *magX, double *magY, double *r, double *p, double *f);
void WaterBag(long n, long *idum, double p0, double r0, double *r, double *p);
//void KineticEnergy(long n, double *energKin, double *p);
void KineticEnergy(long n, long nLoc, double *energKin, double *p);
//void PotentialEnergy(long n, double *energPot, double *r, double magX, double magY);
//void PotentialEnergy(double *energPot, double *r, double magX, double magY);
void PotentialEnergy(double *energPot, double magX, double magY);
//void Force(long n, double *force, double *r, double *magX, double *magY);
void Force(long n, long nLoc, double *force, double *r, double *magX, double *magY);
float ran2(long *idum);

int main(int argc, char **argv)
{
	
	double t_0, t_F, dT;
	t_0 = omp_get_wtime();
	
	long n, nLoc, seed, idum;
	double p0, r0;
	double energKin, energPot, magX, magY, energ, energ0, error;
	double time, finalTime, timeStep, timeOutput, timeCount;
	double *r, *p, *r_, *p_, *force;
	
	MPI_Init(NULL, NULL);
	
	FILE *init, *enrg, *fmag, *finalSpace;
	
	int world_rank; //Qual o rank do proceso no COMM_WORLD
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size; //Quantos processos no COMM_WORLD
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	//Abrindo os arquivos de saída no rank 0
	if(world_rank == 0){
		init = fopen("./initialPhase.dat", "w");
		enrg = fopen("./energy.dat", "w");
		fmag = fopen("./magnet.dat", "w");
		finalSpace = fopen("./finalPhase.dat", "w");
		
		FILE *in = fopen("input.in", "r");
	
		fscanf(in, "%ld", &n);			//Total de partículas
		fscanf(in, "%lf", &finalTime);	//Tempo total de simulação
		fscanf(in, "%lf", &timeStep);	//Passo de tempo
		fscanf(in, "%lf", &timeOutput);	//Passo para arquivo de saída
		fscanf(in, "%lf", &p0);			//?
		fscanf(in, "%lf", &r0);			//?
		fscanf(in, "%ld", &seed);		//?
		
		fclose(in);
		idum = -seed;
		
		energKin = ran2(&idum);
		energPot = ran2(&idum);
		magX = ran2(&idum);
		magY = ran2(&idum);
	}
	//Aguardamos o rank 0 concluir
	MPI_Barrier(MPI_COMM_WORLD);
	//Fazendo um broadcast para os outros rank
	MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&finalTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&timeStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&timeOutput, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&p0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&r0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&seed, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&energKin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&energPot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&magX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&magY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	nLoc = n/world_size;//Quantidade de partículas em cada processo
	
	//Alocando os vetores r, p e force
	if(world_rank == 0){ //no rank 0 eles tem n elementos
		r = (double *)malloc((size_t) (n * sizeof(double)));
		p = (double *)malloc((size_t) (n * sizeof(double)));
	}
	r_ = (double *)malloc((size_t) (nLoc * sizeof(double)));
	p_ = (double *)malloc((size_t) (nLoc * sizeof(double)));
	force = (double *)malloc((size_t) (nLoc * sizeof(double)));
	
	//Inicializando r e p em rank 0
	if(world_rank == 0){
		WaterBag(n, &idum, p0, r0, r, p);
		for (long i = 0; i < n; i++)
		{
			fprintf(init, "%lf\t%lf\n", r[i], p[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//distribuindo r e p para os outros ranks 
	MPI_Scatter(r, nLoc, MPI_DOUBLE, r_, nLoc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(p, nLoc, MPI_DOUBLE, p_, nLoc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	KineticEnergy(n, nLoc, &energKin, p_);
	Force(n, nLoc,force, r_, &magX, &magY);
	PotentialEnergy(&energPot, magX, magY);
	energ0 = energKin + energPot;
	if(world_rank == 0){
		//cout << "Energia Cinetica Inicial: " << energKin << endl;
		printf("Energia Cinetica Inicial: %lf \n", energKin);
		//cout << "Energia Potencial Inicial: " << energPot << endl;
		printf("Energia Potencial Inicial: %lf \n", energPot);
		//cout << "Energia Total Inicial: " << energ0 << endl;
		printf("Energia Total Inicial: %lf \n", energ0);
		//cout << "Magnetizacoes iniciais:   MagX: " << magX << "  MagY: " << magY << endl;
		printf("Magnetizacoes iniciais:   MagX: %lf  MagY: %lf\n", magX, magY);
	} 
	error = .0;
	time = .0;
	timeCount = .0;

	while (time < finalTime)
	{
			
		Integration(n, nLoc, timeStep, &magX, &magY, r_, p_, force);

		time += timeStep;
		timeCount += timeStep;

		if (timeCount >= timeOutput)
		{
			KineticEnergy(n, nLoc, &energKin, p_);
			PotentialEnergy(&energPot, magX, magY);
			energ = energKin + energPot;
			error = (energ - energ0) / energ0;
			error = fabs(error);
			// Colocar aqui um if para parar a simulação quandoo errofor grande 
			// Definir erro limite aceitavel
			//printf("%lf\t%1.2le\t%lf\t%lf\t%lf\n", time, error, npMean, var, statMoment4);
			timeCount = 0.0;
			if(world_rank == 0){			
				fprintf(enrg, "%lf\t%.12lf\t%.12lf\n", time, energKin, energPot);
				fprintf(fmag, "%lf %.12lf %.12lf %.12lf\n", time, magX, magY, sqrt(magX*magX + magY*magY));
				if (error > RMAX){
					printf("%lf\t%1.2le\n", time, error);
				}
			}
			//Aguardamos o rank 0 concluir
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	
	
	if(world_rank == 0){
		printf("Salvando os espacos de fase finais\n");
	}
	
	MPI_Gather(r_, nLoc, MPI_DOUBLE, r, nLoc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(p_, nLoc, MPI_DOUBLE, p, nLoc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	double rr = .0;
	if(world_rank == 0){
		for (long i = 0; i < n; i++)
		{
			rr = r[i];
			while (rr > dpi / 2.){
				rr -= dpi;
			}
			while (rr < -dpi / 2.){
				rr += dpi;
			}
			fprintf(finalSpace, "%lf\t%lf\n", r[i], p[i]);
		}
		free(r);
		free(p);
		fclose(enrg);
		fclose(finalSpace);
		fclose(fmag);
		fclose(init);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	free(r_);
	MPI_Barrier(MPI_COMM_WORLD);
	free(p_);
	MPI_Barrier(MPI_COMM_WORLD);
	free(force);
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(world_rank == 0){
		t_F = omp_get_wtime();
		printf("d_T = %lf\n", t_F - t_0);
	}
	
	MPI_Finalize();

	return 0;
}

void WaterBag(long n, long *idum, double p0, double r0, double *r, double *p)
{

	double aux = .0;
	for (long i = 0; i < n; i++)
	{
		r[i] = ((double)ran2(idum))*r0;
		p[i] = ((double)ran2(idum) - .5)*2.*p0;
		aux += p[i];
	}
	aux = aux / ((double)n);

	for (long i = 0; i < n; i++)
	{
		p[i] -= aux;
	}

	return;
}

void KineticEnergy(long n, long nLoc, double *energKin, double *p)
{
	double energLocal;
	*energKin = energLocal = .0;

	for (long i = 0; i < nLoc; i++)
	{
		energLocal += p[i] * p[i];
	}
	//MPI_Allreduce(&energLocal, energKin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Reduce(&energLocal, energKin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	*energKin = *energKin / ((double)2 * n);
	return;
}

void PotentialEnergy(double *energPot, double magX, double magY)
{
	*energPot = .0;
	*energPot = .5*(1. - magX*magX - magY*magY);

	return;
}

void Force(long n, long nLoc, double *force, double *r, double *magX, double *magY)
{
	double *as = (double *)malloc(n * sizeof(double));
	double *ac = (double *)malloc(n * sizeof(double));
	double aux1, aux2;
	double magX_, magY_;
	aux1 = .0;
	aux2 = .0;
	*magX = magX_ = .0;
	*magY = magY_ = .0;
	for (long i = 0; i < nLoc; i++)
	{
		aux1 = sin(r[i]);
		aux2 = cos(r[i]);
		magX_ += aux2;
		magY_ += aux1;
		as[i] = aux1;
		ac[i] = aux2;
	}
	//MPI_Allreduce(sendbuf,recvbuf,  count,datatype,  op,   comm)
	MPI_Allreduce(&magX_, magX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Reduce(&magX_, magX, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Bcast(magX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Allreduce(sendbuf,recvbuf,  count,datatype,  op,   comm)
	MPI_Allreduce(&magY_, magY, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Reduce(&magY_, magY, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Bcast(magY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	*magX = *magX / ((double)n);
	*magY = *magY / ((double)n);

	for (long i = 0; i < nLoc; i++)
	{
		force[i] = ac[i] * (*magY) - as[i] * (*magX);
	}

	free(as);
	free(ac);

	return;
}
					  
void Integration(long n, long nLoc, double dt, double *magX, double *magY, double *r, double *p, double *f)
{
	long i;
	double mx, my;

	mx = .0;
	my = .0;
	for (i = 0; i<nLoc; i++)
	{
		p[i] += B0*dt*f[i];
		r[i] += D0*dt*p[i];
	}

	Force(n, nLoc, f, r, &mx, &my);
	for (i = 0; i<nLoc; i++)
	{
		p[i] += B1*dt*f[i];
		r[i] += D1*dt*p[i];
	}

	Force(n, nLoc, f, r, &mx, &my);
	for (i = 0; i<nLoc; i++)
	{
		p[i] += B1*dt*f[i];
		r[i] += D0*dt*p[i];
	}

	Force(n, nLoc, f, r, &mx, &my);
	for (i = 0; i<nLoc; i++)
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
