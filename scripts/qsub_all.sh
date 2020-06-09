#!/bin/bash
#
# This script submits the PBS batch job scripts for every directory created by
# pbs_setup.sh, J=*/parity=*1/s_matrix/.
#
# Humberto Jr
# Jun, 2020

set -e
set -u

# Total angular momentum
J_min=0
J_max=1
J_step=1

# Misc
job_filename="job.pbs"
submit_cmd="qsub"

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

start_dir=$PWD

for J in $(seq $J_min $J_step $J_max)
do
	work_dir="$start_dir/J=$J"
	assert_file $work_dir

	for parity in $(seq -1 2 +1)
	do
		if [ "$parity" == "-1" ]
		then
			if [ "$J" == "0" ]
			then
				continue
			fi

			parity_dir="$work_dir/parity=-1"
		else
			parity_dir="$work_dir/parity=+1"
		fi

		assert_file $parity_dir

		batch_script="$parity_dir/$job_filename"
		assert_file $batch_script

		echo "#"
		echo "# $batch_script"

		cd $parity_dir

		$submit_cmd $batch_script

		cd $start_dir
	done
done
