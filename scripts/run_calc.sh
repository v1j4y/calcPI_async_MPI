#!/usr/bin/env bash

for i in {5..22..1}
do
	let nitermax=$(( 2**i ))
	let nsteps=50000
	time mpirun -n 33 $PWD/lib/bin/mpiexample $nitermax $nsteps
done
