#! /usr/bin/perl
# This runs 2 perl script to find non-covalent interaction by every amino acid side and main chain. 
# Usage potential.pl protein_name coordinate(.gro) topology(.top) mdpfile(.mdp) trajectory(.trr) starting_time(in picosecond) 
use strict;
use File::Path qw(make_path remove_tree);
#Call perl scripts.
if (system ("perl rerun.pl $ARGV[4] $ARGV[0]/$ARGV[0]\_enrggrps.txt $ARGV[0]/$ARGV[0] $ARGV[2] $ARGV[1] $ARGV[3]")){
die "Unable to rerun gromacs. Make sure that gromacs is in path or bring gromacs in path. Exiting\n";
}
wait;
if (system ("perl calculate_energy.pl $ARGV[0]/$ARGV[0] $ARGV[5]")){
die "Unable to calculate energy.Exiting\n";
}
