# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error = job.err.$(jobid)
# @ output = job.out.$(jobid)
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 16
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 4:00:00
# @ notification = complete
# @ notify_user = derrien@fzu.cz
# @ queue

# run the program

module purge
module load intel impi
echo "/!\ NFFT is missing! Compile it !!"
#module load nfft/3.2.4
# LD_LIBRARY_PATH=/scratch/tjyderri/nfft/3.2.4/lib:$LD_LIBRARY_PATH
RunningTrunc=/ptmp/tjyderri/octopus.git/bin

# cd /ptmp/${USER}/
poe -n $NSLOTS $RunningTrunc/octopus >& $LOADL_JOB_NAME.log
