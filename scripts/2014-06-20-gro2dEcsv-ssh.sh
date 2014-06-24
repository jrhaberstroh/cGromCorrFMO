#!/bin/bash

# USEAGE:
# Edit script vairables for the location of the files, then pass in the numbers of the files to convert.
# Make sure that $num_frames is equal to the number of frames per file!

for i in $@; do
    JOB=$HOME/mass-storage/JobArchive/2014-06-20-FMO50ns-dEcsv
    PROJECT=$HOME/Code/photosynth/cGromCorrFMO
    SSH=jhabers@hopper.nersc.gov:/global/scratch2/sd/jhabers/2014-06-04-FMOeq/md/gro100
    
    fmo_top_path=$JOB/FMO_conf/4BCL_pp.top
    bcx_itp_path=$JOB/FMO_conf/amber_mod.ff/bcx.itp
    #output_path=$PROJECT/data/CsvOut/2014B_MDPost/md$i

    num_frames=10000

    fmo_ssh_path=$SSH/fmo50ns-$i.gro
    fmo_target_path=$JOB/fmo50ns-$i.gro
    output_path=$JOB/dEcsv$i
    
    echo "Using file $fmo_traj_path, and computing dEcsv with $num_frames frames"
    
    scp $fmo_ssh_path $fmo_target_path
    $PROJECT/src/traj2dEcsv $fmo_target_path $fmo_top_path $bcx_itp_path $output_path $num_frames
    rm $fmo_target_path
done
