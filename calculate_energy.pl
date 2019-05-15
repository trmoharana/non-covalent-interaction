#! /usr/bin/perl
# Programe to find energy of each index group
# Usage calculate_energy.pl protein_name/protein_name
use strict;
use warnings;
use IPC::Open2;
use IO::Handle;  # so we can call methods on filehandles
my %energy;
my @residue;
my @energy;
my @write;
my $side_int_coul;
my $side_int_lj;
my $side_int_total;
my $side_ext_coul;
my $side_ext_lj;
my $side_ext_total;
my $side_total;
my $main_int_coul;
my $main_int_lj;
my $main_int_total;
my $main_ext_coul;
my $main_ext_lj;
my $main_ext_total;
my $main_total;
my $command;

# Creat output file.
open(GIV, ">>$ARGV[0]\_potential.txt")|| die "you don't have permision to write in this path";
print GIV "# Aminoacid_number \t Aminoacid_name \t mainchain_internal_coulombic \t mainchain_internal_lj \t mainchain_internal_total \t mainchain_external_coulombic \t mainchain_external_lj \t mainchain_external_total \t mainchain_total \t sidechain_internal_coulombic \t sidechain_internal_lj \t sidechain_internal_total \t sidechain_external_coulombic \t sidechain_external_lj \t sidechain_external_total \t sidechain_total\n";
close GIV;
open(ENRG, "$ARGV[0]\_enrggrps.txt")|| die "$ARGV[0]\_enrggrps.txt doesn't exist. Run rerun.pl with appropriate arguments";

# Read energygroup file.
while (<ENRG>)
{
@residue=split();
if($#residue==1){

# Main chain atom potential energy contributions.
$command="gmx energy -f $ARGV[0]\_$residue[0]\_m.edr -s $ARGV[0]\_$residue[0]\_m.tpr -b $ARGV[1] -o dummy.xvg \n";
open2( my $out, my $in, $command ) or die "Can't open $command: $!";

# Set both filehandles to print immediately and not wait for a newline. Just a good idea to prevent hanging. 
print $in "Coul-SR:$residue[1]\_main-$residue[1]\_main \n Coul-SR:$residue[1]\_main-$residue[1]\_main\_com \n Coul-14:$residue[1]\_main-$residue[1]\_main \n Coul-14:$residue[1]\_main-$residue[1]\_main\_com \n LJ-SR:$residue[1]\_main-$residue[1]\_main \n LJ-SR:$residue[1]\_main-$residue[1]\_main\_com \n LJ-14:$residue[1]\_main-$residue[1]\_main \n LJ-14:$residue[1]\_main-$residue[1]\_main\_com \n \n";
close $in;
while (<$out>){

# Extract energy terms.
if ( m/kJ\/mol/){
@energy=split();
$energy[0]=~s/://;
$energy{$energy[0]}=$energy[1];
undef(@energy);
$main_int_coul=$energy{"Coul-SR$residue[1]\_main-$residue[1]\_main"}+$energy{"Coul-14$residue[1]\_main-$residue[1]\_main"};
$main_int_lj=$energy{"LJ-SR$residue[1]\_main-$residue[1]\_main"}+$energy{"LJ-14$residue[1]\_main-$residue[1]\_main"};
$main_int_total=$main_int_coul+$main_int_lj;
$main_ext_coul=$energy{"Coul-SR$residue[1]\_main-$residue[1]\_main\_com"}+$energy{"Coul-14$residue[1]\_main-$residue[1]\_main\_com"};
$main_ext_lj=$energy{"LJ-SR$residue[1]\_main-$residue[1]\_main\_com"}+$energy{"LJ-14$residue[1]\_main-$residue[1]\_main\_com"};
$main_ext_total=$main_ext_coul+$main_ext_lj;
$main_total=$main_int_total+$main_ext_total;
}}

# Write energy terms to output file.
unlink "dummy.xvg";
open(OUTPUTFILE, ">>$ARGV[0]\_potential.txt");
print OUTPUTFILE " $residue[0] \t $residue[1] \t $main_int_coul \t $main_int_lj \t $main_int_total \t $main_ext_coul \t $main_ext_lj \t $main_ext_total \t $main_total ";
close OUTPUTFILE;
undef $main_int_coul;
undef $main_int_lj;
undef $main_int_total;
undef $main_ext_coul;
undef $main_ext_lj;
undef $main_ext_total;
undef $main_total;
# Glycine doesnot have side chain. Adding end line. 
if(/GLY/){
open(OUTPUTFILE, ">>$ARGV[0]\_potential.txt");
print OUTPUTFILE " \n ";
close OUTPUTFILE;
}else{

#Side chain atom potential energy contributions
$command="gmx energy -f $ARGV[0]\_$residue[0]\_s.edr -s $ARGV[0]\_$residue[0]\_s.tpr -b $ARGV[1] -o dummy.xvg \n";
open2( my $out, my $in, $command ) or die "Can't open $command: $!";
print $in "Coul-SR:$residue[1]\_side-$residue[1]\_side \n Coul-SR:$residue[1]\_side-$residue[1]\_side\_com \n Coul-14:$residue[1]\_side-$residue[1]\_side \n Coul-14:$residue[1]\_side-$residue[1]\_side\_com \n LJ-SR:$residue[1]\_side-$residue[1]\_side \n LJ-SR:$residue[1]\_side-$residue[1]\_side\_com \n LJ-14:$residue[1]\_side-$residue[1]\_side \n LJ-14:$residue[1]\_side-$residue[1]\_side\_com \n \n";
close $in;
while (<$out>){
if ( m/kJ\/mol/){
@energy=split();
$energy[0]=~s/://;
$energy{$energy[0]}=$energy[1];
undef(@energy);
$side_int_coul=$energy{"Coul-SR$residue[1]\_side-$residue[1]\_side"}+$energy{"Coul-14$residue[1]\_side-$residue[1]\_side"};
$side_int_lj=$energy{"LJ-SR$residue[1]\_side-$residue[1]\_side"}+$energy{"LJ-14$residue[1]\_side-$residue[1]\_side"};
$side_int_total=$side_int_coul+$side_int_lj;
$side_ext_coul=$energy{"Coul-SR$residue[1]\_side-$residue[1]\_side\_com"}+$energy{"Coul-14$residue[1]\_side-$residue[1]\_side\_com"};
$side_ext_lj=$energy{"LJ-SR$residue[1]\_side-$residue[1]\_side\_com"}+$energy{"LJ-14$residue[1]\_side-$residue[1]\_side\_com"};
$side_ext_total=$side_ext_coul+$side_ext_lj;
$side_total=$side_int_total+$side_ext_total;
}}
unlink "dummy.xvg";
open(OUTPUTFILE, ">>$ARGV[0]\_potential.txt");
print OUTPUTFILE " \t $side_int_coul \t $side_int_lj \t $side_int_total \t $side_ext_coul \t $side_ext_lj \t $side_ext_total \t $side_total \n";
close OUTPUTFILE;
undef $side_int_coul;
undef $side_int_lj;
undef $side_int_total;
undef $side_ext_coul;
undef $side_ext_lj;
undef $side_ext_total;
undef $side_total;
}
}
}
close ENRG;
