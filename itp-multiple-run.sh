#!/bin/bash
#SBATCH --partition=itp
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
 
export OMP_NUM_THREADS=1
export MV2_ENABLE_AFFINITY=0


#running masses from 0.01 to 3
for mass in {1..300}
	
	do
		sleep 0.5
		#running r1 from 0 to 2
		for r1 in $(seq 0 20)

		do
			srun -n 1 python anicomp.py $r1 20 $mass 0.01 &
		done
	done

# Wait for all child processes to terminate.
wait