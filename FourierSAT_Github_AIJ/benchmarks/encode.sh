function run(){
for file in `ls $1`
do
    if [ -d $1"/"$file ] 
    then
        run $1"/"$file $2 $3
    else
        filename=$(basename $1"/"$file)
	if [[ ($filename != *.cnf) && ($filename != *.opb)  ]]
        then 
            continue 
        fi
        filepath=$1/$file
        python3 encoding_pbo.py $filepath $2/$file.opb
    fi
done 
}

run cubic_vector_covering cubic_vector_covering_pb




