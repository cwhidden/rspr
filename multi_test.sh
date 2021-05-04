#!/bin/bash
make debug
#if (! test -d "test_output" ) then
 #  mkdir test_output
#fi
./rspr -multi_approx < test_trees/multi_tree_basic.txt
echo ""
./rspr -multi_approx < test_trees/multi_tree_basic_5.txt
echo "Finished"
