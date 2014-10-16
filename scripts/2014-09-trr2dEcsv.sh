#!/bin/bash
f=$@

trr_files=$(printf "%s\n" ${f[@]} | sort -n -t _ -k 2)

for trr_file in $trr_files; do
    Tns=".1"
    dtfs=5

    FMO_conf=$HOME/Jobs/2014-08-15-4BCL/FMO_conf
    fmo_top_path=$FMO_conf/4BCL-pp.top
    bcx_itp_path=$FMO_conf/amber_mod.ff/bcx_cdc.itp
    traj2dEcsv=/home/jhaberstroh/Code/photosynth/cGromCorrFMO/src/traj2dEcsv

    Tfs=$(printf '%.0f' "$Tns"E6)
    Tps=$(echo "$Tfs / 1000" | bc)
    num_frames=$(echo "$Tfs / $dtfs" | bc)
    trr_dir=$(dirname $trr_file)
    trr_base=$trr_dir/$(basename $trr_file .trr)
    csv_out=$trr_base.csv
    echo $trr_base: "$Tns"ns, NumFrames=$num_frames
    #trjconv -f $trr_base.trr -o $trr_base.gro -s $trr_base.tpr -e $Tps -n $FMO_conf/index.ndx  <<< 22

    echo "Using file $trr_base.gro, and computing dEcsv with $num_frames frames"
    echo "Making csv file at $csv_out"
    #$traj2dEcsv $trr_base.gro $fmo_top_path $bcx_itp_path $csv_out $num_frames
    rm $trr_base.gro

    #csv_out=.csv
    python 2014-06-csv2hdf5.py $HOME/mass-storage/data/2014-08-4BCL/dETimeseries.hdf5 timeseries $csv_out.time -frames_per_file $num_frames -dt .005
done
