#!/bin/bash
echo "Building..."
make rspr


declare -a tests=(
    		  "trees2.txt"
		  "trees3.txt"
		  "trees4.txt"
		  "trees5.txt"
		  "trees6.txt"
		  "trees_100_9.txt"
		  "trees_100_17.txt"
		  #"trees_100_24.txt"
		  "multi_tree_basic_2.txt"
		  "multi_tree_basic_5.txt"
		  "multi_tree_7.1_test_00.txt"
		  "multi_tree_7.3_test_00.txt"
		  "multi_tree_7.4_test_00.txt"
		  "multi_tree_8.2_test_00.txt"
		  "multi_tree_8.4_test_00.txt"
		  "multi_tree_8.4_test_01.txt"
		  "multi_tree_8.5_test_00.txt"
		  "multi_tree_8.6_test_00.txt"
		  "multi_tree_8.7_test_00.txt"
		  "multi_tree_8.7_test_01.txt"
		  "cluster_test"
		  )
declare -a binary_tests_array=(
                  "trees2.txt"
		  "trees3.txt"
		  "trees4.txt"
		  "trees5.txt"
		  "trees6.txt"
		  "trees_100_9.txt"
		  "trees_100_17.txt"
		  "trees_100_24.txt"
		  "cluster_test"
		  "rho_test.txt"
		  "rho_test2.txt"
		  "cluster_1.txt"
		  "cluster_2.txt"
		  "cluster_3.txt"
		  "cluster_4.txt"
		  )

declare -a move_tests=(
    		  "trees2.txt"
		  "trees3.txt"
		  "multi_tree_7.1_test_00.txt"
		  "multi_tree_7.3_test_00.txt"
		  "multi_tree_7.4_test_00.txt"
		  "multi_tree_8.2_test_00.txt"
		  "multi_tree_8.4_test_00.txt"
		  "multi_tree_8.4_test_01.txt"
		  "multi_tree_8.5_test_00.txt"
		  "multi_tree_8.6_test_00.txt"
		  "multi_tree_8.7_test_00.txt"
		  "multi_tree_8.7_test_01.txt"
		  "cluster_test"
)

cluster="-bb"
binary_tests=false
reverse=""
leaf_reduction=""
show_moves=""
binary="-multifurcating"
incremental_tests=false
clean_incremental=false

for arg in "$@"
do
    if [ "$arg" = "-cluster" ] ;    then
	cluster="-cluster_test -leaf_reduction"
    elif [ "$arg" = "-binary_tests" ] ;    then
	binary_tests=true
    elif [ "$arg" = "-binary" ] ;    then
	binary=""
    elif [ "$arg" = "-reverse" ] ;    then
	reverse="-debug_reverse_trees"
    elif [ "$arg" = "-leaf_reduction" ] ;    then
	leaf_reduction=true
    elif [ "$arg" = "-show_moves" ] ;    then
	show_moves="-show_moves"
    elif [ "$arg" = "-incremental" ] ;    then
	incremental_tests=true
    elif [ "$arg" = "-clean_incremental_test_output" ] ;    then
	clean_incremental=true

    fi
done

elapsed=0

if [ "$clean_incremental" = true ] ; then
    rm ../../incremental_output/*.csv
    exit 1
fi

if [ "$incremental_tests" = true ] ; then
    TIMEFORMAT="%U"
    for d in `seq 0 20`
    do
	filename=../../incremental_output/100_leaf_"$d"_spr_timings.csv
	if [ ! -f $filename ]
	then
	    touch $filename
	fi
	for i in `seq 1 30`
	do
	    #python3 ../../randomMultifurcatingTree.py 1000 | ./rspr -multifurcating -random_spr 5 | time ./rspr -multifurcating -leaf_reduction
	    #time ls > /dev/null | tr -d . > ../../timeout
	    { time python3 ../../randomMultifurcatingTree.py 100 | ./rspr -multifurcating -random_spr "$d" | time ./rspr -multifurcating -leaf_reduction 2> ../../sleep.stderr > /dev/null; } 2>> $filename
	    #printf ',' >> ../../time.txt
	done
	echo "Finished dSPR test $d"
    done
    exit 1
fi

if [ "$binary_tests" = true ] ; then
   for i in ${binary_tests_array[@]}
   do
       echo "\n\n////////////////////////////////////////////////////////"
       echo $i
       echo "////////////////////////////////////////////////////////\n\n"
       time ./rspr $show_moves $binary $reverse $cluster< test_trees/$i
   done
else
   for i in ${tests[@]}
   do
       echo "\n\n////////////////////////////////////////////////////////"
       echo $i
       echo "////////////////////////////////////////////////////////\n\n"
       time ./rspr -multifurcating $cluster $reverse < test_trees/$i    
   done
fi
echo "Finished"
