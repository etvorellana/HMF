EXECS=hmf hmf_mpi
MPICC=mpicc
CC=gcc
all: ${EXECS}
hmf: main.c
	${CC} -o hmf main.c -O3 -lm

hmf_mpi: main-openmpi.c
	${MPICC} -o hmf_mpi main-openmpi.c -O3 -lm

clean:
	rm -f ${EXECS}
	rm -f *.dat
