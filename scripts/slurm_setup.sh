#!/bin/bash
#
# This script prepares slurm batch job scripts for those calculations performed in
# the workflow created by abc_workflow.sh. Scripts are written in each directory
# J=*/parity=*1/s_matrix/.
#
# Humberto Jr
# Jun, 2020

set -u
set -e

# Total angular momentum
J_min=0
J_max=10
J_step=1

# Executable
abc_exe="/home/humberto/balalab/H+D2/abc/bkmp2/exe/abc.out"

# Slurm batch configuration
wall_time="10:00:00"
max_memory="8G"
omp_threads=32
queue_name=""
modules=""

# Libraries to load by LD_LIBRARY_PATH (format: 'path_a:path_b:path_c:etc')
env_ld_path=""

# OpenMP configuration
env_omp_threads="OMP_NUM_THREADS=$omp_threads"

# GPU configuration (by generic resource scheduling)
gres="gpu:k80:4"

# Misc
input_filename="input.d"
slurm_filename="job.sh"

################################################################################

assert_file ()
{
	if [ ! -e $1 ]
	then
		echo
		echo "$0, error: $1 not found"
		echo
		exit 666
	fi
}

assert_file $abc_exe

if [ "$env_ld_path" != "" ]
then
	env_ld_path='export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:'$env_ld_path'"'
fi

for J in $(seq $J_min $J_step $J_max)
do
	work_dir="$PWD/J=$J"

	for parity in $(seq -1 2 +1)
	do
		if [ "$J" == "0" ] && [ "$parity" == "-1" ]
		then
			continue
		else
			if [ "$parity" == "-1" ]
			then
				parity_dir="$work_dir/parity=-1"
				job_name="J=-$J"
			else
				parity_dir="$work_dir/parity=+1"
				job_name="J=+$J"
			fi
		fi

		assert_file $parity_dir

		filename="$parity_dir/$slurm_filename"

		echo "#!/bin/sh"                                           > $filename
		echo '#SBATCH --job-name="'$job_name'"'                   >> $filename
		echo "#SBATCH --mem=$max_memory"                          >> $filename
		echo "#SBATCH --time=$wall_time"                          >> $filename
		echo "#SBATCH --nodes=1"                                  >> $filename
		echo "#SBATCH --ntasks=1"                                 >> $filename
		echo "#SBATCH --cpus-per-task=$omp_threads"               >> $filename

		if [ "$gres" != "" ]
		then
			echo "#SBATCH  -p GPU"                                 >> $filename
			echo "#SBATCH --gres=$gres"                            >> $filename
		fi

		echo ""                                                   >> $filename
		echo "# Job script generated at $(date) by $0"            >> $filename

		echo ""                                                   >> $filename
		echo "export $env_omp_threads"                            >> $filename

		if [ "$env_ld_path" != "" ]
		then
			echo $env_ld_path                                      >> $filename
		fi

		if [ "$modules" != "" ]
		then
			echo ""                                                >> $filename
			echo "module load $modules"                            >> $filename
		fi

		echo ""                                                   >> $filename
		echo 'abc_exe="'$abc_exe'"'                               >> $filename

		echo ""                                                   >> $filename
		echo 'input="'$input_filename'"'                          >> $filename

		echo ""                                                   >> $filename
		echo 'echo "# Calculation starting at $(date)" > abc.log' >> $filename

		echo 'echo "" >> abc.log'                                 >> $filename

		echo ""                                                   >> $filename

		echo '$abc_exe < $input >> abc.log'                       >> $filename

		echo ""                                                   >> $filename

		echo 'echo "" >> abc.log'                                 >> $filename

		echo 'echo "# Calculation ending at $(date)" >> abc.log'  >> $filename
	done

	echo "J=$J"
done
