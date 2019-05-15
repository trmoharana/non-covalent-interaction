#! /usr/bin/perl
# This runs all other perl script to find non-covalent interaction by every amino acid main and side chains.  
# Creates one index group for each amino acids and file containing name of energy groups.  
# Usage master.pl -p protein_name -s coordinate(.gro) -t topology(.top) -m mdpfile(.mdp) -r trajectory(.trr) -b starting_time(in picosecond) -i run_input(.tpr)
use strict;
use File::Path qw(make_path remove_tree);
my $protein_name="Unnamedprotein";
my $coordinate;
my $topology;
my $mdpfile;
my $trajectory;
my $time;
my $runinput;
my $i;
#Check for inputs
if ("@ARGV" !~ /.gro/){
die "Coordinate file(.gro) not specified";
}
if ("@ARGV" !~ /.top/){
die "Topology file(.top) not specified";
}
if ("@ARGV" !~ /.mdp/){
die "mdp file(.mdp) not specified";
}
if ("@ARGV" !~ /.trr/){
die "Trajectory file(.xtc) not specified";
}
if ("@ARGV" !~ /.tpr/){
die "Run input file (.tpr) not specified";
}
#Assign file name
for (my $i=0; $i<$#ARGV; $i++){
if ($ARGV[$i]eq"-p"){
$protein_name=$ARGV[$i+1];
}elsif ($ARGV[$i]eq"-s"){
$coordinate=$ARGV[$i+1];
}elsif ($ARGV[$i]eq"-t"){
$topology=$ARGV[$i+1];
}elsif ($ARGV[$i]eq"-m"){
$mdpfile=$ARGV[$i+1];
}elsif ($ARGV[$i]eq"-r"){
$trajectory=$ARGV[$i+1];
}elsif ($ARGV[$i]eq"-b"){
$time=$ARGV[$i+1];
}elsif ($ARGV[$i]eq"-i"){
$runinput=$ARGV[$i+1];
}}
#Creat separate directory
unless(make_path $protein_name){
die "Unable to create $protein_name. Either the directory already exist or you don't have writting permision\n";
}
#Check for gromacs versions. 
system("gmx check");
while(<out>){
if ( m/VERSION/){
	open(VER, ">$protein_name\_analysis.txt");
	print VER "# $_ used for analysis \n";
	close VER;
	last; 
	}
}
system("gmx check -f $trajectory");
while(<out>){
if ( m/VERSION/){
	open(VER, ">>$protein_name\_analysis.txt");
	print VER "# $_ used for production run \n";
	close VER;
	last; 
	}
}
#Print input variables
print "Protein name : $protein_name \n topology: $topology \n structure : $coordinate \n trajectory : $trajectory \n mdp: $mdpfile\n Start time: $time\n";
open(FINAL, ">$protein_name/$protein_name\_energy.txt");
print FINAL "# Protein name : $protein_name \n# topology: $topology \n# structure : $coordinate \n# trajectory : $trajectory \n# mdp: $mdpfile \n# first frame: $time\n";
close FINAL;
#Creat error file.
open (ERROR, ">$protein_name/error.txt");
close ERROR;
#Call perl scripts.
if (system ("perl prepare_indexgrp.pl $coordinate $protein_name/$protein_name")){
remove_tree $protein_name;
die "Unable to prepare index group. Exiting\n$!\n";
}
wait;
if (system ("perl potential.pl $protein_name $coordinate $topology $mdpfile $trajectory $time")){
die "Unable to calculate potential energy.Exiting\n $!\n";
}
wait;
print "Succesfuly calculated ennthalpy contribution of each amino acid side chain. Results are stored in $protein_name/$protein_name\_energy.txt and errors are stored in $protein_name/error.txt\n"; 



