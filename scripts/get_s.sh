#!/bin/bash
#
# This script reads and sorts ABC's log file for those calculations performed in
# the workflow created by abc_workflow.sh. Each pair of channels in the S-matrix
# as function of collision energy is written in a unique data file and stored in
# directories named J=*/parity=*1/s_matrix/.
#
# NOTE: the previous content of J=*/parity=*1/s_matrix/, if any, is deleted when
# running. If the directory does not exist it is created.
#
# Humberto Jr
# Jun, 2020

set -u
set -e

# Total angular momentum
J_min=0
J_max=1
J_step=1

# Misc
abc_log="abc.log"
buffer=".buffer.dat"
list=".list.dat"

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

build_dir ()
{
	if [ -d $1 ]
	then
		rm -rf $1/*
	else
		mkdir $1
	fi
}

header="#  Coll.                   a    v    j    k    a'   v'   j'   k'     Re(S)     Im(S)     |S|^2"

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

		s_dir="$parity_dir/s_matrix"
		build_dir $s_dir

		declare -a a
		declare -a v
		declare -a j
		declare -a k

		flag=""

		grep "E(eV) =" $input | cut -d "=" -f 2 | cut -d "S" -f 1 > $list

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
			fi
		done < $input

		n_max=$n

		for n in $(seq 0 1 $n_max)
		do
			if [ "${a[$n]}" == "2" ]
			then
				break
			fi

			ch_a="${a[$n]}    ${v[$n]}    ${j[$n]}    ${k[$n]}"

			for m in $(seq $n 1 $n_max)
			do
				ch_b="${a[$m]}    ${v[$m]}    ${j[$m]}    ${k[$m]}"

				grep "$ch_a    $ch_b" $input > $buffer

				if [ "$(wc -l < $buffer)" != "0" ]
				then
					filename="$s_dir/ch=$(($n + 1))-$(($m + 1)).dat"

					paste $list $buffer > $filename

					sed -i "1s/^/$header\n/" $filename
					sed -i "s/                       //g" $filename

					rm $buffer
					echo $filename
				fi
			done
		done

		rm $list

		unset a
		unset v
		unset j
		unset k
	done
done
