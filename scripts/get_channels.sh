#!/bin/bash
#
# This script reads and sorts ABC's log file for those calculations performed in
# the workflow created by abc_workflow.sh. Each channel is listed for each J and
# parity.
#
# Humberto Jr
# Jun, 2020

set -u
set -e

# Total angular momentum
J_min=0
J_max=10
J_step=1

# Misc
abc_log="abc.log"

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

echo "#	J	p	a	v	j	k	E(eV)  	Ch."

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
			else
				parity_dir="$work_dir/parity=+1"
			fi
		fi

		input=$parity_dir/$abc_log
		assert_file $input

		declare -a a
		declare -a v
		declare -a j
		declare -a k
		declare -a e

		flag=""

		while IFS= read line
		do
			if [ "$flag" != "Channel" ]
			then
				flag=$(echo $line | cut -c -7)
				continue
			fi

			if [ "$flag" == "Channel" ]
			then
				if [ "$line" == "" ]
				then
					break
				fi

				if [ "${line:1:1}" == "-" ]
				then
					continue
				fi

				n=$(echo $line | awk '{print $1}')
				n=$(($n - 1))

				a[$n]=$(echo $line | awk '{print $2}')
				v[$n]=$(echo $line | awk '{print $3}')
				j[$n]=$(echo $line | awk '{print $4}')
				k[$n]=$(echo $line | awk '{print $5}')
				e[$n]=$(echo $line | awk '{print $6}')
			fi
		done < $input

		n_max=$n

		for n in $(seq 0 1 $n_max)
		do
			 echo "	$J	$parity	${a[$n]}	${v[$n]}	${j[$n]}	${k[$n]}	${e[$n]}	$((n + 1))"
		done

		unset a
		unset v
		unset j
		unset k
		unset e
	done
done
