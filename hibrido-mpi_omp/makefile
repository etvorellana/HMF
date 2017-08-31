EXECS=hmf_mpiomp
MPICC=mpicc
CC=gcc
all: ${EXECS}

hmf_mpiomp: main-mpiomp.c
	${MPICC} -o hmf_mpiomp main-mpiomp.c -O3 -lm -fopenmp

clean:
	rm -f ${EXECS}
	rm -f *.dat
