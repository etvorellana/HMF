#!/bin/bash
#SBATCH -J hmf_mpi2				# nome do job
#SBATCH --partition long		# fila onde o job será enfileirado, short ou long
#SBATCH --nodes 1				# quantidade de nós, 1 a 20
#SBATCH --ntasks 1				# quantidade de tarefas, 8 cpus por nó
#SBATCH --cpus-per-task 8		# quantidade de cpus por tarefa
##SBATCH --time 3:00:00			# tempo limite, se omitido usará o tempo total da fila
#SBATCH --exclusive             # se necessário, ativa uso exclusivo dos nós

source /usr/local/mpi/openmpi/2.1.0/openmpi_vars.sh

cp -f ../inputFiles/input01.in ./input.in
for i in 1 2 3
do
	mpirun -np 2 ./hmf_mpi
done
tar -czvf example_01_mpi2.tar.gz *.dat --remove-files

cp -f ../inputFiles/input01.2n ./input.in
for i in 1 2 3
do
	mpirun -np 2 ./hmf_mpi
done
tar -czvf example_02_mpi2.tar.gz *.dat --remove-files

cp -f ../inputFiles/input03.in ./input.in
for i in 1 2 3
do
	mpirun -np 2 ./hmf_mpi
done
tar -czvf example_03_mpi2.tar.gz *.dat --remove-files

cp -f ../inputFiles/input04.in ./input.in
for i in 1 2 3
do
	mpirun -np 2 ./hmf_mpi
done
tar -czvf example_04_mpi2.tar.gz *.dat --remove-files

cp -f ../inputFiles/input05.in ./input.in
for i in 1 2 3
do
	mpirun -np 2 ./hmf_mpi
done
tar -czvf example_05_mpi2.tar.gz *.dat --remove-files


