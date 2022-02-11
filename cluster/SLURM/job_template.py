# This string forms the first lines of every job file when using Sun Grid Engine scheduler. MEMORY is replaced by
# RosettaDockMemory or ProtacModelMemory. Instead of HEADER the content of the file which SCHEDULER_PARAMS points to
# is injected (line is removed if SCHEDULER_PARAMS is not set).
job_template = '#!/bin/bash\n' \
                '#SBATCH -J JOB_NAME\n' \
                '#SBATCH -o JOB_NAME.out\n' \
                '#SBTACH --cpus-per-task=CPUS_PER_TASK\n' \
                '#SBTACH --mem-per-cpu=MEM_PER_CPU\n' \
                'HEADER\n' \
                'source /etc/slurm.conf\n' \
                'echo `hostname`\n' \
                'cd $PBS_O_WORKDIR\n'