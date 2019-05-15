# run main chain and side chain separately 
#! /usr/bin/perl 
# Programe to rerun GROMACS to to calculate potential energy of each index group.
# Usage rerun.pl input_file_name(.xtc/.trr) energygroup_filename(.txt) index_filename topology(.top) structure(.gro) mdp(.mdp) 
use strict;
use warnings; 
no strict 'refs';
use File::Copy qw(move);
my @line; 
my @unlink;
my $command;
my $engrgrp;
open (READ, "$ARGV[1]");
while (<READ>){
@line = split();
if ($#line){
	print "reruning for protein $ARGV[0] and amino acid @line  $line[1]\n";
# Modify mdp file	
	$engrgrp=0;
	$ARGV[5]=~s/.mdp//;
	my $old= "$ARGV[5]\_1.mdp";
	my $new= "$ARGV[5]\.mdp";
	open INPUT,  '<',  "$new" or die "Can't read old file: $!";
	open OUTPUT, '>', "$old" or die "Can't write new file: $!";
	while( <INPUT> )
    {
	if ( m/energygrps /) {
	$_="energygrps               =	$line[1]\_main	$line[1]\_main\_com";
	$engrgrp=1;
	}
	print OUTPUT "$_";
	}
	close INPUT;
	close OUTPUT;
	move $old, $new;
	unless ($engrgrp){
	open (OUT, ">>$ARGV[5]\.mdp");
	print OUT "energygrps               =	$line[1]\_main	$line[1]\_main\_com \n";
	close OUT;	
	}
# Creat .tpr
	$command = "gmx grompp -f $ARGV[5]\.mdp -n $ARGV[2]\_$line[0].ndx -c $ARGV[4] -p $ARGV[3] -o $ARGV[2]\_$line[0]\_m.tpr";
	system ($command);
	wait;
# Call mdrun rerun
	$command = "gmx mdrun -deffnm $ARGV[2]\_$line[0]\_m -rerun $ARGV[0]";
	system ($command);
	wait;
	@unlink=("$ARGV[2]\_$line[0]\_m\.log", "$ARGV[2]\_$line[0]\_m\.trr", "mdout.mdp");
	unlink @unlink;
#side chain energy
	unless ($line[1]=~m/GLY/){
	open INPUT,  '<',  "$new" or die "Can't read old file: $!";
	open OUTPUT, '>', "$old" or die "Can't write new file: $!";
	while( <INPUT> )
    {
	if ( m/energygrps /) {
	$_="energygrps               =	$line[1]\_side	$line[1]\_side\_com";
	}
	print OUTPUT "$_";
	}
	close INPUT;
	close OUTPUT;
	move $old, $new;
# Creat .tpr
	$command = "gmx grompp -f $ARGV[5]\.mdp -n $ARGV[2]\_$line[0].ndx -c $ARGV[4] -p $ARGV[3] -o $ARGV[2]\_$line[0]\_s.tpr";
	system ($command);
	wait;
# Call mdrun rerun
	$command = "gmx mdrun -deffnm $ARGV[2]\_$line[0]\_s -rerun $ARGV[0]";
	system ($command);
	wait;
	@unlink=("$ARGV[2]\_$line[0]\_s\.log", "$ARGV[2]\_$line[0]\_s\.trr", "mdout.mdp");
	unlink @unlink;
	}
}
}
close READ;
