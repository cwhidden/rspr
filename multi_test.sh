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
		  "cluster_1.txt"
		  "cluster_2.txt"
		  "cluster_3.txt"
		  #"cluster_4.txt"
		  "cluster_5.txt"
		  "cluster_6.txt"
		  "cluster_7.txt"
		  "cluster_8.txt"
		  "cluster_9.txt"
		  "cluster_a.txt"
		  "cluster_b.txt"
		  "cluster_c.txt"
		  "cluster_d.txt"
		  "cluster_e.txt"
		  #"cluster_f.txt"
		  "cluster_10.txt"
		  "cluster_11.txt"

		  "cluster_12.txt"
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
		  "cluster_5.txt"
		  "cluster_6.txt"
		  "cluster_7.txt"
		  "cluster_8.txt"
		  "cluster_9.txt"
		  "cluster_a.txt"
		  "cluster_b.txt"
		  "cluster_c.txt"
		  "cluster_d.txt"
		  "cluster_e.txt"
		  "cluster_f.txt"
		  "cluster_10.txt"
		  "cluster_11.txt"
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

cluster=""
binary_tests=false
reverse=""
leaf_reduction=""
show_moves=false
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
	show_moves=true
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


if [ "$binary_tests" = true ] ; then
   for i in ${binary_tests_array[@]}
   do
       echo "\n\n////////////////////////////////////////////////////////"
       echo $i
       echo "////////////////////////////////////////////////////////\n\n"
       time ./rspr  $binary $reverse $cluster< test_trees/$i
   done
else
   for i in ${tests[@]}
   do
       echo "\n\n////////////////////////////////////////////////////////"
       echo $i
       echo "////////////////////////////////////////////////////////\n\n"
       time ./rspr -multi_4_branch $cluster $reverse  < test_trees/$i    
   done
fi
echo "Finished"
