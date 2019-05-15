#! /usr/bin/perl
# Programe to calculate number of mutations surrounding a given amino acid
# This program takes output of find_neighbour1.1.pl as input
# Usage count_mutations.pl input.txt output.txt -ignore (side or main) -positions 1 2 3 
use warnings;
my $neighbour; my $contri; my $number; my $ig; my $index; my $total_num; my $total_con; my $total_rat; my $con_rat; my @neighbour; my @positions; my @temp; my %neighbour; my %count;
# Obtain mutation position
for (my $i=0;$i<=$#ARGV;$i++){
if ($ARGV[$i] eq "-positions"){
$index=$i+1;
}
}
for (my $i=$index;$i<=$#ARGV;$i++){
push @positions, $ARGV[$i];
}
# Read input and write output
open (WRITE, ">$ARGV[1]");
print WRITE "# @ARGV \n# Amino_acid 	 Neighbours(occupancy) 	 Expected_energy 	 Actual_energy 	 Difference(observed-expected) total_neighbours mutant_neighbour total_contribution mutant_contribution number_rat contribution_ratio \n";  
open (READ, "$ARGV[0]");
while (<READ>){
if (/^#/){
next;
}
if (/^\s/){
next;
}
@temp=split();
$neighbour=$temp[1];
chomp();
print WRITE;
$neighbour=~s/\(/start/g;
$neighbour=~s/\)/end/g;
@neighbour=split('end', $neighbour);
# Remove ignored neighbours
for (my $i=0;$i<=$#ARGV;$i++){
if ($ARGV[$i] eq "-ignore"){
$index=$i;
}
}
$ig=$ARGV[$index+1];
@neighbour = grep !/$ig/, @neighbour;
foreach my $a (@neighbour){
@temp=split('start', $a); 
# Remove amino acid name and keep only positional information
$temp[0]=~ s/[^0-9]//g;
$neighbour{$temp[0]}=$temp[1];
}
@temp=keys%neighbour;
$contri=0; $number=0; $total_num=0; $total_con=0;
foreach my $b (@temp){
$total_num++; 
$total_con=$total_con+$neighbour{$b};
}
foreach my $c (@positions){
if (grep {$_ eq $c} @temp){
$number++;
$contri=$contri+$neighbour{$c};
}
}
$total_rat=$number/$total_num;
$con_rat=$contri/$total_con;
print WRITE "\t $total_num \t $number \t $total_con \t $contri \t $total_rat \t $con_rat \n";
undef $total_rat; 
undef $con_rat;
undef $neighbour;
undef @neighbour;
undef @temp;
undef %neighbour;
}
close READ;
close WRITE;
