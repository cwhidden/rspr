#!/bin/bash
echo "Building..."
make
#if (! test -d "test_output" ) then
 #  mkdir test_output
#fi

declare -a tests=(
    		  "trees2.txt"
		  "trees3.txt"
		  "trees4.txt"
		  "trees5.txt"
		  "trees6.txt"
		  "trees_100_9.txt"
		  "trees_100_17.txt"
		  "multi_tree_basic_2.txt"
		  "multi_tree_basic_5.txt"
		  "multi_tree_7.1_test_00.txt"
		  "multi_tree_7.3_test_00.txt"
		  "multi_tree_7.4_test_00.txt"
		  "multi_tree_8.2_test_00.txt"
		  "multi_tree_8.4_test_00.txt"
		  "multi_tree_8.5_test_00.txt"
		  )
		  
for i in ${tests[@]}
do
    echo "\n\n////////////////////////////////////////////////////////"
    echo $i
    echo "////////////////////////////////////////////////////////\n\n"

    ./rspr -multifurcating < test_trees/$i
done

echo "Finished"
