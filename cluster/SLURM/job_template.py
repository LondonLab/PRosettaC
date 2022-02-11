# This string forms the first lines of every job file when using Sun Grid Engine scheduler. MEMORY is replaced by
# RosettaDockMemory or ProtacModelMemory. Instead of HEADER the content of the file which SCHEDULER_PARAMS points to
# is injected (line is removed if SCHEDULER_PARAMS is not set).
job_template = '#!/bin/bash\n' \
                '#SBATCH -J PRosettaC_JOB_NAME\n' \
                '#SBATCH -o PRosettaC_JOB_NAME.out\n' \
                '#SBTACH --cpus-per-task=1\n' \
                'HEADER\n' \
                'source /etc/slurm.conf\n' \
                'echo `hostname`\n' \
                'cd $PBS_O_WORKDIR\n'