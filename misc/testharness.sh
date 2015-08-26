#!/bin/bash

declare -a test_names
declare -a test_processes
declare -a test_config_file
declare -a test_num_timesteps
declare -a test_control_checkpoint

numberOfTests=0

runspecifictest() {
	echo "========================================================================================="
	echo "Running test "$4
	echo -e "=========================================================================================\n"
	mpiexec -np $1 ./monc --config=$2 --checkpoint_file="testharness.nc" --checkpoint_unique_per_dump=.false.	 --nn_timesteps=$6
	diff_lines=`diff <(ncdump testharness.nc) <(ncdump $5/$3) | wc -l`
	echo -e "\n========================================================================================="
	if [ $diff_lines -gt 23 ]
	then 
		echo -e "Failure for test '"$4"'\nComparing testharness.nc against "$3" resulted in "$diff_lines" differences">&2
		echo "========================================================================================="
		exit
	else
		echo "Success for test "$4
		rm testharness.nc
		echo "========================================================================================="
	fi
}

main() {
	if [ -z "$MONC_ROOT" ]; then echo "Variable MONC_ROOT is unset, this must point to the directory with the MONC executable"; exit; fi
	parseTests $1
	pushd $MONC_ROOT > /dev/null
	for (( i=0; i<$numberOfTests; i++ ))
	do
		runspecifictest "${test_processes[i]}" "${test_config_file[i]}" ${test_control_checkpoint[i]} "${test_names[i]}" $2 "${test_num_timesteps[i]}"
	done
	popd > /dev/null
	echo -e "\n\n#########################################################################"
	echo "Test harness complete, "$numberOfTests" tests completed successfully"
	echo "#########################################################################"
}

parseTests() {
	file=$1

	while IFS=: read col1 col2 col3 col4 col5
	do
		test_names[numberOfTests]="${col1}"
		test_processes[numberOfTests]="${col2}"
		test_config_file[numberOfTests]="${col3}"
		test_num_timesteps[numberOfTests]="${col4}"
		test_control_checkpoint[numberOfTests]="${col5}"
		((numberOfTests++))
	done < $file
}

if [ $# -ne 2 ]; then echo "Two command line arguments are needed; the testing file and path to the control NetCDF files" ; exit;  fi
main $1 $2