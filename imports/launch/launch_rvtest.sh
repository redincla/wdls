#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=rvtests_HEV
#SBATCH -t 2-00:00
#SBATCH --mem-per-cpu=1G

module add Development/java/1.8.0_232

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar run /home/credin/scratch/WGS/wdls/imports/workflows/rv-tests.wdl -i /home/credin/scratch/WGS/wdls/imports/config-files/rv-tests.inputs.json -o /home/credin/scratch/WGS/wdls/imports/config-files/rv-tests.options.json
