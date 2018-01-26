#!/bin/bash
path=`pwd`
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
    printf -v i "%02d" $i
    cd $path
    echo "Starting $i"
    sed -r --in-place=.backup "s/atom_stress_([0-9]+).out/atom_stress_${i}.out/g" settings.py
    python run_sample.py
    echo "=============Finished Running for $i==============="
    cd $path/dump/results
    ls *.out | xargs -I {} mv {} hehe/${i}_{}
done