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
		  "trees_100_24.txt"
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
		  
for i in ${tests[@]}
do
    echo "\n\n////////////////////////////////////////////////////////"
    echo $i
    echo "////////////////////////////////////////////////////////\n\n"
    if [ -n "$1" ]
       then
	   if [ $1 == "-cluster" ]
	   then
	       time ./rspr -multifurcating -cluster_test -leaf_reduction2 < test_trees/$i
	   fi
    else
	time ./rspr -multifurcating < test_trees/$i    
    fi
done

echo "Finished"
