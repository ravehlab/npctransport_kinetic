#!/bin/bash 

RAN_factor=1

while getopts r: flag
do
    case "${flag}" in
        r) RAN_factor=${OPTARG};;
    esac
done

mkdir grid
for x_start in 5.0 10.0 15.0; do 
    for scale in 2.5 3.0 3.5 4.0 4.5 5.0; do 
        x_end=`echo $x_start $scale | awk '{printf "%.1f", $1*$2}'`
        tag=x_${x_start}_to_${x_end}_ran_factor_$RAN_factor
	echo TAG $tag
        cd $tag
        python3 ../../aggregate_results.py 100 $x_start $x_end
        cd ..
        cp ${tag}/stats_for_graph_${x_start}_${x_end}.pkl grid/${x_start}_${x_end}_$RAN_factor.pkl 
    done
done
