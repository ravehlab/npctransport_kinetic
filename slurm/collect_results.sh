#!/bin/bash 

mkdir grid
for x_start in 20 30 40; do 
    for scale in 1.5 2.5 5.0 10.0; do 
        x_end=$(echo $x_start $x_scale | awk '{printf "%4.1f\n", $1*$2}')
        tag=x_${x_start}_to_${x_end}
        cp ${tag}/stats_for_graph_${x_start}_${x_end}.pkl grid/${x_start}_${x_end}.pkl 
    done
done
