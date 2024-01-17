# DE_pipeline


Steps to run:

pipeline_CIRIquant.py
-p  path to pipeline folder (Default=/home/thaight/projects/rrg-zovoilis/thaight/CIRIquant)
-i  path to project list
-t  organism type (rn6/rn7/mm10/tair10/hg38)
-m  mode to run (CIRI/DE/FIGURES)
-a  account to charge (Default=rrg-zovoilis)

In essence only -i, -t, -m are hard required to run

init_pipeline.sh
$1  list
$2  type
$3  mode

Sbatch this init_pipeline.sh file since certain modes can have high memory intensity

