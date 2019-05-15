#! /usr/bin/perl
# Programe to extract energy contribution of different types of amino acids 
# This program takes output of master.pl as one of its input
# Usage list.pl input_file_name output_file_name number_of_atom end_of_protein
use strict;
use warnings;
#check whether all input file exist
for (my $i=0;$i<$#ARGV;$i++){
if (-e "$ARGV[$i]\_energy.txt"){
print"$ARGV[$i]";
}else{die "\nInput file $ARGV[$i] doesn't exist. Run terminated.";}
}
my $a;
$a=Extract("ALA");
$a=Extract("ARG");
$a=Extract("ASN");
$a=Extract("ASP");
$a=Extract("CYS");
$a=Extract("GLN");
$a=Extract("GLU");
$a=Extract("HIS");
$a=Extract("ILE");
$a=Extract("LEU");
$a=Extract("LYS");
$a=Extract("MET");
$a=Extract("PHE");
$a=Extract("PRO");
$a=Extract("SER");
$a=Extract("THR");
$a=Extract("TRP");
$a=Extract("TYR");
$a=Extract("VAL");
sub Extract {
my $temp=join("", @_);
open(WRITE, ">$temp.txt") || die "Couldn,t creat file $temp.txt";
print WRITE "Protein_name \t Aminoacid_number \t Aminoacid_name \t internal_coulombic \t internal_lj \t external_coulombic \t external_lj \t total_internal \t total_external \n";
close WRITE;
print "Calculating for $temp... \n";
for (my $j=0; $j<$#ARGV+1;$j++){
open(READ, "$ARGV[$j]\_energy.txt");
print "Reading $ARGV[$j]\_energy.txt \n";
while (<READ>){
if ( m/$temp/){
open(APPEND, ">>$temp.txt") || die "Couldn,t creat file $temp.txt";
print APPEND "$ARGV[$j] \t $_";
close APPEND;
}}
close READ;
}
}
