#!/bin/bash

for i in 1 2 3 4 5
do
	mkdir run$i

	cd run$i
	cp ../input/data.lmp data.lmp
	cp ../input/runMD runMD
	cp ../input/forcefield.data forcefield.data
	cp ../input/simulation.in simulation.in
	# Random Seed setting
	randomNumber=$(shuf -i 1-100 -n1)
	sed -i 's/VALUE/'$randomNumber'/' simulation.in
	sed -i 's/JOB_NAME/'$model'_run_'$i'/' runMD
	sbatch runMD
	cd ..

	echo "setting $model run $i completed, seed "$randomNumber"."
done
