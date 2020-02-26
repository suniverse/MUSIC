#!/bin/sh

# Directives
#SBATCH -J MUSIC
#SBATCH -o MUSIC.o%j
#SBATCH -o MUSIC.e%j
#SBATCH -n 1            # total number of mpi tasks requested
#SBATCH --ntasks-per-node=28   # max 14 per node
#SBATCH -t 2:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=xinyu.li@columbia.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -p cca --qos=cca

./MUSIC enzo.conf >> est.out
python normalize.py >& est.out
#cp ./ic.enzo/Grid.* ../run/Cosmology/
# End of script
