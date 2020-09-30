job_template = ['#!/bin/bash\n',
                '#$ -pe smp ',
                '\n',
                # TODO: memory
                'source /etc/sge.conf\n',
                'echo `hostname`\n',
                'cd $PBS_O_WORKDIR\n'
]
