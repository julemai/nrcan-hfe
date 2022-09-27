#!/bin/bash

# submit with:
#       sbatch submit-analyse-occurrence.sh

#SBATCH --account=rpp-julemai                      # your group
#SBATCH --mem-per-cpu=1G                           # memory; default unit is megabytes
#SBATCH --mail-user=juliane.mai@uwaterloo.ca       # email address for notifications
#SBATCH --mail-type=FAIL                           # email send only in case of failure
#SBATCH --time=0-12:00:00                          # time (DD-HH:MM:SS);
#SBATCH --job-name=analyse_occurence               # name of job in queque
#SBATCH --array=1-200


# job-id  :: ${SLURM_ARRAY_JOB_ID}
# task-id :: ${SLURM_ARRAY_TASK_ID}   # this will give you a counter from 1 to 10

module purge
module load StdEnv/2020 netcdf gcc/9.3.0 gdal/3.0.4
module load mpi4py proj
module load python/3.8

# load Python env
source /project/6070465/julemai/nrcan-hfe/env-3.8/bin/activate

# change to code
cd /project/6070465/julemai/nrcan-hfe/src/


# ---------------------
# run all
# ---------------------
nfeatures=1854
ntasks=200         # make sure this is the number of tasks set for array-job

features=( $( seq $(( ${SLURM_ARRAY_TASK_ID} - 1 )) ${ntasks} $(( ${nfeatures} -1 )) ) )  # (2,12,22,...)
features=$( printf "%s," "${features[@]}" )   # "2,12,22,"
features=$( echo ${features::-1} )            # "2,12,22"

# runs the script with set of features
python analyse_occurrence.py -i "${features}" --bbox_buffer 0.5 --dates_buffer 5.0,5.0 --tmpdir "/project/6070465/julemai/nrcan-hfe/data/output/"



# # ---------------------
# # run the ones that had an error
# # ---------------------
# features=( 261 461 484 561 661 684 687 732 761 861 884 887 932 961 1061 1084 1087 1132 1161 1261 1284 1287 1332 1339 1340 1361 1374 1461 1484 1487 1532 1539 1540 1561 1574 1661 1684 1687 1712 1732 1739 1740 1758 1761 1764 1774 )  # 46 total
# ifeature=$( echo ${features[$(( ${SLURM_ARRAY_TASK_ID} - 1 ))]} )   # This will give you one basin after the other

# # runs the script with the ith feature
# python analyse_occurrence.py -i "${ifeature}" --bbox_buffer 0.5 --dates_buffer 5.0,5.0 --tmpdir "/project/6070465/julemai/nrcan-hfe/data/output/"



# ntasks=200 --> 9 or 10 occurrences each

# JOBID
# 65508874  - testing     ::    6 occurrences;   3 tasks -->   2     occurrences per task --> took 1h
# 65513537  - final       :: 1854 occurrences; 200 tasks --> 9 or 10 occurrences per task --> took ~2h   (some geomet will fail because nodes do not have access to internet)
# 65531234  - redo failed ::   46 occurrences;  46 tasks -->   1     occurrences per task --> took ???   (some geomet will fail because nodes do not have access to internet)
# 65534294  - final       :: redo with --dates_buffer 5.0,5.0



# post-analyse results
#
# -------------------------
# get results sorted by accumulated_mm
# ---> some are really large --> event longer than >>1 year
# ---> some are really small --> there is just no event shortly before and after (you might want to post-process them again)
# -------------------------
# files=$( \ls ../data/output/analyse_occurrence_*/*.json )
# for ff in $files ; do mm=$( jq . $ff | grep "accumulated_mm" | cut -d : -f 2 | cut -d ',' -f 1 ) ; echo ${mm} ${ff} ; done | sort -n