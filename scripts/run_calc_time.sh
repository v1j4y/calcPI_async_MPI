#!/usr/bin/env bash

for i in {1..5..1}
do
	let nitermax=240000000
	let nsteps=50
	let ncpu=$(( 2**i + 1))
	echo " ----- $ncpu ------ "
	time mpirun -n $ncpu $PWD/lib/bin/mpiexample $nitermax $nsteps
done
