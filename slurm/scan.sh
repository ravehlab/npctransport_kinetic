#!/usr/local/bin/bash 

BASE_FOLDER="`pwd -P`"

#for x_start in 10.0 20.0 30.0 40.0 100.0; do 
#    for scale in 1.5 2.5 5.0 10.0 15.0 20.0; do 
for x_start in 5.0 10.0 15.0; do 
    for scale in 2.5 3.0 3.5 4.0 4.5 5.0; do 
        x_end=`echo $x_start $scale | awk '{printf "%.1f", $1*$2}'`
        tag=x_${x_start}_to_${x_end}
        echo $tag
	if [ ! -e "$tag" ] ; then
            mkdir "$tag"
	fi
        cd "$tag"
        bash ../run_all_mw.sh $x_start $x_end
	cd "$BASE_FOLDER"
    done
done
