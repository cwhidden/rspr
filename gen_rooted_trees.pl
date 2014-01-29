#!/usr/bin/perl

################################################################################
# gen_rooted_trees.pl
#
# Output to stdout all rootings of an unrooted binary newick tree
# from stdin
#
# Example
# STDIN:
# (1, 2, 3);
# STDOUT:
# (1,(2,3));
# (2,(1,3));
# (3,(1,2));
#
# Works by finding 1,2,3 subtrees:
# assume 1 is a leaf
# 1st tree = 1 , (2,3)
# recurse on left and right subtrees as unrooted, ie
# 1=(1,2) 2 = left of 3, 3 = right of 3
# 1=(1,3) 2 = left of 2, 3 = right of 2
#
# Copyright 2009 Chris Whidden
# cwhidden@dal.ca
# http://kiwi.cs.dal.ca/Software/RSPR
# September 7, 2009
# 
# This file is part of rspr.
# 
# rspr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rspr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rspr.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

use strict;

sub unrooted_move($$);
sub get_subtrees($);


while ($#ARGV >= 0) {
if ($ARGV[0] eq "-h" || $ARGV[0] eq "--help") {
    print "Usage: gen_rooted_trees.pl\n";
    print "Output to stdout all rootings of an unrooted binary newick tree\n";
	print "from stdin.\n";
	print "\n";
	print "Example\n";
	print "STDIN:\n";
	print "(1, 2, 3);\n";
	print "STDOUT:\n";
	print "(1,(2,3));\n";
	print "(2,(1,3));\n";
	print "(3,(1,2));\n";

    exit(0);
}
}

while(<>) {
	my $unrooted_tree = $_;
	exit if ($unrooted_tree =~ "NULL");
	my @subtrees = get_subtrees($unrooted_tree);
	unrooted_move($unrooted_tree, 0);
	unrooted_move($unrooted_tree, 1);
	unrooted_move($unrooted_tree, 2);
			
};

# Subroutines
sub unrooted_move($$) {
	my ($unrooted_tree, $dir) = @_;
	#print "u=$unrooted_tree  d=$dir\n";
	my @subtrees = get_subtrees($unrooted_tree);

	# generate a rooted tree based on the move
	my $rooted_tree = "(";
	$rooted_tree .= $subtrees[$dir];
	$rooted_tree .= ",(";
	my @new_dir = (0 .. $dir-1, $dir+1 .. $#subtrees);
	$rooted_tree .= $subtrees[$new_dir[0]];
	$rooted_tree .= ",";
	$rooted_tree .= $subtrees[$new_dir[1]];
	$rooted_tree .= "));";
	print "$rooted_tree\n";

	my @new_subtrees = get_subtrees($subtrees[$dir]);
	return if ($#new_subtrees <= 0);
	my $new_unrooted_tree = "((";
	$new_unrooted_tree .= $subtrees[$new_dir[0]];
	$new_unrooted_tree .= ",";
	$new_unrooted_tree .= $subtrees[$new_dir[1]];
	$new_unrooted_tree .= "),";
	$new_unrooted_tree .= $new_subtrees[0];
	$new_unrooted_tree .= ",";
	$new_unrooted_tree .= $new_subtrees[1];
	$new_unrooted_tree .= ");";

	unrooted_move($new_unrooted_tree, 1);
	unrooted_move($new_unrooted_tree, 2);
}

sub get_subtrees($) {
	my $unrooted_tree = shift;
	my @subtrees = ();
	my $subtree = "";
	my $parity = 0;
	while($unrooted_tree =~ /./g) {
		my $char = $&;
		if ($char eq "," && $parity == 1) {
			push(@subtrees, $subtree);
			$subtree = "";
		}
		else {
			if ($char eq ")") {
				$parity--;
			}
			if ($parity > 0) {
				$subtree .= $char;
			}
			if ($char eq "(") {
				$parity++;
			}
		}
		#print "c=$char  p=$parity s=$subtree\n";
	}
	if ($subtree ne "") {
		push(@subtrees, $subtree);
	}
	return @subtrees;
}


