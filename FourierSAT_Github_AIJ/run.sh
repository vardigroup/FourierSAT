#!/bin/bash

function run(){
TL=60
NC=1

for file in `ls $1`
do
    if [ -d $1"/"$file ] 
    then
        run $1"/"$file $2 $3
    else
        filename=$(basename $1"/"$file)
        if [[ $filename != *.cnf ]]
        then 
            continue 
        fi
        filepath=$1/$file
        echo $file
	start=$(date +%s)
        outpath="res/"$3"_"$2".txt"
        if [[ $2 == 'cms' ]];
        then
            timeout $TL cryptominisat5 $filepath  > /dev/null 2>&1
	elif [[ $2 == 'walksat' ]];
        then
            timeout $TL solvers_binary/WalkSAT2013 $filepath 101 0.567
	elif [[ $2 == 'naps' ]];
	then
            timeout $TL ../binary_of_other_solvers/naps  $filepath 
	elif [[ $2 == 'roundingsat' ]];
	then
            timeout $TL ../binary_of_other_solvers/roundingsat < $filepath 
	elif [[ $2 == 'openwbo' ]];
	then
            timeout $TL ../open-wbo/open-wbo -formula=1 $filepath 
        elif [[ $2 == 'fouriersat' ]];
        then
            timeout $TL  python3 FourierSAT.py --timelimit 600 --cpus $NC --verbose 1 $filepath
        else 
            timeout $TL solvers_binary/$2 $filepath > /dev/null 2>&1
	fi
	status=$?
        if [ $status -eq 124 ]
        then
            echo "time $filename TO" >> $outpath
        else
           end=$(date +%s)
           time=$(( $end - $start ))
           echo "time $filename $time" >> $outpath
        fi
    fi
done 
}

#run benchmarks/cubic_vector_covering/ fouriersat vc



