#!/bin/bash
set -u
set -e

# Isotopes
atom_a=1
atom_b=2
atom_c=2

# Total angular momentum
J_min=0
J_max=10
J_step=1

# Diatomic parity and rovib. states
jpar=1
jmax=8
kmax=$jmax
emax="3.5"

# Sectors
rmax="12.0"
mtr=200

# Collision energies
enrg="0.91916"
dnrg="0.000004308665"
nnrg=2000

# Output
nout=1000
jout=$jmax

# Misc
nmax=4000
input="input.d"
pes_extern_data=

################################################################################

build_dir ()
{
	if [ -d $1 ]
	then
		rm -rf $1/*
	else
		mkdir $1
	fi
}

for J in $(seq $J_min $J_step $J_max)
do
	work_dir="$PWD/J=$J"

	build_dir $work_dir

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

		build_dir $parity_dir

		filename="$parity_dir/$input"

		echo "&input"                              > $filename
		echo "  mass = $atom_a, $atom_b, $atom_c" >> $filename
		echo "  jtot = $J"                        >> $filename
		echo "  ipar = $parity"                   >> $filename
		echo "  jpar = $jpar"                     >> $filename
		echo "  jmax = $jmax"                     >> $filename
		echo "  kmax = $kmax"                     >> $filename
		echo "  rmax = $rmax"                     >> $filename
		echo "   mtr = $mtr"                      >> $filename
		echo "  emax = $emax"                     >> $filename
		echo "  enrg = $enrg"                     >> $filename
		echo "  dnrg = $dnrg"                     >> $filename
		echo "  nnrg = $nnrg"                     >> $filename
		echo "  nout = $nout"                     >> $filename
		echo "  jout = $jout"                     >> $filename
		echo "  nmax = $nmax"                     >> $filename
		echo "&end"                               >> $filename

		if [ ! "$pes_extern_data" == "" ]
		then
			cp -f $pes_extern_data $parity_dir/
		fi
	done

	echo "J=$J"
done
