EXECS= hmf_omp
CC=gcc

all: ${EXECS}

hmf_omp: main-openmp.c openmpFuncs.c openmpLibs.h
	${CC} -o hmf_omp main-openmp.c openmpFuncs.c -O3 -lm -fopenmp

clean:
	rm -f ${EXECS}
	rm -f *.dat
