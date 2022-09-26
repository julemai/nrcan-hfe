#!/bin/bash

# submit with:
#       sbatch submit-analyse-event.sh

#SBATCH --account=rpp-julemai                      # your group
#SBATCH --mem-per-cpu=16G                           # memory; default unit is megabytes
#SBATCH --mail-user=juliane.mai@uwaterloo.ca       # email address for notifications
#SBATCH --mail-type=FAIL                           # email send only in case of failure
#SBATCH --time=0-06:00:00                          # time (DD-HH:MM:SS);
#SBATCH --job-name=analyse_event                   # name of job in queque
#SBATCH --array=1-100


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
nfeatures=363
ntasks=100         # make sure this is the number of tasks set for array-job

features=( $( seq $(( ${SLURM_ARRAY_TASK_ID} - 1 )) ${ntasks} $(( ${nfeatures} -1 )) ) )  # (2,12,22,...)
features=$( printf "%s," "${features[@]}" )   # "2,12,22,"
ifeatures=$( echo ${features::-1} )            # "2,12,22"

# runs the script with set of features
python analyse_event.py -i "${ifeatures}" --bbox_buffer 0.5 --dates_buffer 5.0,5.0 --tmpdir "/project/6070465/julemai/nrcan-hfe/data/output/"



# # ---------------------
# # run the ones that had an error
# # ---------------------
# features=( 95 195 295 112 212 312 123 223 323 129 229 329 130 230 330 133 233 333 179 279 188 288 227 327  )  # 24 total
# features=( 133  )  # 1 total
# ifeatures=$( echo ${features[$(( ${SLURM_ARRAY_TASK_ID} - 1 ))]} )   # This will give you one basin after the other

# # runs the script with the ith feature
# python analyse_event.py --ifeatures "${ifeatures}" --bbox_buffer 0.5 --dates_buffer 5.0,5.0 --tmpdir "/project/6070465/julemai/nrcan-hfe/data/output/"



# ---------------------
# zip results (because they can contain quite a lot of PNGs)
# ---------------------
cd /project/6070465/julemai/nrcan-hfe/data/output/
ifeatures_list=$( tr -s ',' ' ' <<< "${ifeatures}" )    # "2 12 22"
for ifeature in ${ifeatures_list} ; do
    if [ -e "analyse_event_${ifeature}.zip" ] ; then rm "analyse_event_${ifeature}.zip" ; fi
    if [ -e "analyse_event_${ifeature}" ] ; then zip -r "analyse_event_${ifeature}.zip"  "analyse_event_${ifeature}" ; fi
    rm -r "analyse_event_${ifeature}"
done
cd -



# JOBID
# 65572013  - final       :: 363 events; 100 tasks --> 3 or 4 events per task --> took ~2h   (some geomet will fail because nodes do not have access to internet)
# 65573326  - redo        :: 24 events where requested memory of 4GB was not enough --> increase to 8GB
# 65576177  - redo        :: 363 events; 100 tasks --> 3 or 4 events per task --> took ~2h   (fixed that for some stations no event found; memory now set to 16GB)
#
