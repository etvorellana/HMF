#!/bin/bash
#SBATCH -J hmf_seq				# nome do job
#SBATCH --partition long		# fila onde o job será enfileirado, short ou long
#SBATCH --nodes 1				# quantidade de nós, 1 a 20
#SBATCH --ntasks 1				# quantidade de tarefas, 8 cpus por nó
#SBATCH --cpus-per-task 8		# quantidade de cpus por tarefa
##SBATCH --time 3:00:00			# tempo limite, se omitido usará o tempo total da fila
#SBATCH --exclusive             # se necessário, ativa uso exclusivo dos nós


cp -f ../inputFiles/input01.in ./input.in
for i in 1 2 3
do
	srun ./hmf
done
tar -czvf example_01.tar.gz *.dat --remove-files

cp -f ../inputFiles/input01.2n ./input.in
for i in 1 2 3
do
	srun ./hmf
done
tar -czvf example_02.tar.gz *.dat --remove-files

cp -f ../inputFiles/input03.in ./input.in
for i in 1 2 3
do
	srun ./hmf
done
tar -czvf example_03.tar.gz *.dat --remove-files

cp -f ../inputFiles/input04.in ./input.in
for i in 1 2 3
do
	srun ./hmf
done
tar -czvf example_04.tar.gz *.dat --remove-files

cp -f ../inputFiles/input05.in ./input.in
for i in 1 2 3
do
	srun ./hmf
done
tar -czvf example_05.tar.gz *.dat --remove-files

