
determine_if_finished() {
	terminated_run=0
	local search_line=`grep "Model run complete due to model time" $1`
	local found_cont=`echo "$search_line" | wc -c`
	if [ $found_cont -gt 1 ]; then
		local mtime=`echo "$search_line" | awk '{ print $9 }'`
		echo "Terminating chain run as MONC has exceeded termination time, model time is $mtime seconds"
		terminated_run=1
	else 
		local search_line=`grep "messages file containing termination command" $1`
		local found_cont=`echo "$search_line" | wc -c`
		if [ $found_cont -gt 1 ]; then
			local mtime=`echo "$search_line" | awk '{ print $15 }'`
			echo "Terminating chain run as MONC was instructed to finish, model time is $mtime seconds"
			terminated_run=1
		else
			local search_line=`grep "timestep completion" $1`
    	local found_cont=`echo "$search_line" | wc -c`
			if [ $found_cont -gt 1 ]; then
				local mtime=`echo "$search_line" | awk '{ print $12 }'`
				echo "Terminating chain run as MONC exceeded timestep limit, model time is $mtime seconds"
				terminated_run=1
			fi
		fi
	fi
}

RUN_MONC_CONFIG=0
RUN_MONC_CP=0
outputid=0

run_monc() {
	if [ ! -f $MONC_EXEC ]; then
		echo "Error - executable $MONC_EXEC does not exist"
		exit
	fi
	if [ ! -z "$crun" ] && [ $crun -ge $MAX_CONTINUATION_RUNS ]; then
		echo "This has been run $crun times which exceeds your configured maximum number of runs"
		exit
	fi

	local output_filename=`ls -rt1 $STDOUT_DIR/output_$RUN_NAME* 2> /dev/null | tail -1`
	local checkpoint_filename=`ls -rt1 $CP_DIR/$RUN_NAME*.nc 2> /dev/null | tail -1`

	if [ ! -z "$output_filename" ] && [ ! -z "$checkpoint_filename" ]; then
		determine_if_finished $output_filename
		if [ $terminated_run -eq 0 ]; then
			outputid=`sed 's/.*_//' <<< "$output_filename"`
			if [ -z "$crun" ] || [ $crun -ne $outputid ]; then
				if [ -z "$cpfile" ] || [ "$cpfile" != "$checkpoint_filename" ]; then
					RUN_MONC_CP=1
				else
					echo "Not running MONC as the latest checkpoint file is the same that the previous run executed with"
					exit
				fi
			else
				echo "Not running MONC as there is no new output from the previous run, there was probably a submission error"
				exit			
			fi
		fi
        elif [ ! -z "$checkpoint_filename" ] && [ -z "$crun" ] && [ -z "$cpfile" ]; then
		RUN_MONC_CONFIG=2
	else
		if [ -z "$crun" ]; then
			RUN_MONC_CONFIG=1
		else
			echo "Error, this is configured as a continuation run but output and/or checkpoint file not found, check your script parameters"
			exit
		fi
	fi

	if [ $RUN_MONC_CONFIG -ge 1 ] || [ $RUN_MONC_CP -eq 1 ]; then
		export OMP_NUM_THREADS=1
		export MPICH_MAX_THREAD_SAFETY=multiple
                
		local submittedId=$(qsub -W depend=afterany:$PBS_JOBID -v crun=$outputid,cpfile=$checkpoint_filename $SUBMISSION_SCRIPT_NAME)

		((outputid++))
		local outputfn=$STDOUT_DIR"/output_"$RUN_NAME$outputid

                echo "This cycle is controlled by:$SUBMISSION_SCRIPT_NAME" > $outputfn
                echo "This cycle job:$PBS_JOBID:$PBS_JOBNAME" >> $outputfn
                echo "Next cycle job:$submittedId" >> $outputfn
                echo "" >> $outputfn

                echo ""

                # Cold start
		if [ $RUN_MONC_CONFIG -eq 1 ]; then
    		    echo "Start MONC with configuration file $TESTCASE"
		    eval 'aprun -B $MONC_EXEC --config=$TESTCASE >> $outputfn 2>&1'

                # Reconfiguration
                elif [ $RUN_MONC_CONFIG -eq 2 ]; then
                    echo "Reconfigure MONC with configuration file:"
                    echo "  $TESTCASE and its linked xml file,"
                    echo "  starting from checkpoint file:"
                    echo "     $checkpoint_filename"
                    eval 'aprun -B $MONC_EXEC --reconfig=$TESTCASE --checkpoint=$checkpoint_filename --retain_model_time=.true. >> $outputfn 2>&1'

                # Restart
		else
		    echo "Restarting MONC with checkpoint file $checkpoint_filename"
		    eval 'aprun -B $MONC_EXEC --checkpoint=$checkpoint_filename >> $outputfn 2>&1'
  	fi
fi
}
