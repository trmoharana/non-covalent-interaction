#! /usr/bin/perl
# Programe to find atoms within given distance from given atoms
# Usage find_neighbour.pl 0trajectory(.trr) 1runinput(.tpr) 2starting_time(pico second) 3Energygroup_filename 4distance(nm) 5output_file(without .txt) 6protein_name 7secondary_structure 8-ignoress(if) 9-bypasatomsel(if) 10-expected expected_energy 11-ignore aminoacids 
use strict;
use warnings;
my @ignore;
my @coil;
my @helix;
my @sheet;
my %atmfrqnc;
my %grpfrqnc;
my %expected;
my %actual;
my $firstaa;
my $lastaa;
my $out;
my $expected;
for (my $i=0;$i<=$#ARGV;$i++){
if ($ARGV[$i] eq "-ignore"){
for (my $j=$i+1;$j<=$#ARGV;$j++){
@ignore=(@ignore, "$ARGV[$j]");
}}}
for (my $i=0;$i<=$#ARGV;$i++){
if ($ARGV[$i] eq "-expected"){
$expected=$ARGV[$i+1];
}}
$out="$ARGV[6]/$ARGV[5]\_$expected\_report";
if (grep {$_ eq "-ignoress"} @ARGV){
$out="$ARGV[6]/$ARGV[5]\_$expected\_ignoress";
}
if(grep {$_ eq "-ignore"} @ARGV){
my $temp=join '_', @ignore;
$out="$out\_$temp";
}
$out="$out\.txt";
# Read expected energy contribution by different amino acid main and side chain.
open (DATA, "$expected\.txt");
while (<DATA>){
my @read = split();
$expected{"$read[0]"}="$read[1]";
}
close DATA;
open (REPORT, ">$out");
print REPORT "# @ARGV \n # Amino_acid \t Neighbours(occupancy) \t Expected_energy \t Actual_energy \t Difference\(observed-expected\) \n";
close REPORT;
# Read observed energy contribution by different amino acid main and side chain.
open (READ, "$ARGV[6]/$ARGV[6]\_potential.txt");
while (<READ>){
if (/^#/){
next;
}
my @observ = split();
# Find first and last amino acids.
unless ($lastaa){
$firstaa=$observ[1];
}
$lastaa=$observ[1];
$actual{"$observ[1]\_main"}=$observ[8];
unless ( m/GLY/){
$actual{"$observ[1]\_side"}=$observ[15];
}}
close READ;
# Obtain secondary structure information (list of amino acids in the helix, sheet and coil confwermation. See example_ss.txt).
unless(grep {$_ eq "-ignoress"} @ARGV){ #?
open (SS, "$ARGV[7]");
while (<SS>){
my @temp=split;
my $compare=shift @temp;
if ($compare eq "Helix" || $compare eq "helix"){
@helix=@temp;
}
if ($compare eq "Sheet" || $compare eq "sheet"){
@sheet=@temp;
}
if ($compare eq "Coil" || $compare eq "coil"){
@coil=@temp;
}}
close SS;
# Check for duplicate assignment
for (my $i=0; $i<=$#helix; $i++){
if (grep {$_ eq $helix[$i]} @coil){
die "You have assign $helix[$i] both as helix and coil";
}
if (grep {$_ eq $helix[$i]} @sheet){
die "You have assign $helix[$i] both as helix and sheet";
}
}
for (my $i=0; $i<=$#sheet; $i++){
if (grep {$_ eq $sheet[$i]} @coil){
die "You have assign $sheet[$i] both as sheet and coil";
}}}
# Find atoms within neighbourhood.
my $grpname;
open (INPUT, "$ARGV[3]") || die "Couldn,t open file $ARGV[3]";
while (<INPUT>){
if( m/GLY/){
next;
}
$grpname=(split)[1];
unless (grep {$_ eq "-bypasatomsel"} @ARGV){ #correct it.
my $cmd="gmx select -f $ARGV[0] -s $ARGV[1] -pbc -rmpbc -b $ARGV[2] -n $ARGV[6]/$ARGV[6]\_aminoacids.ndx -selrpos atom -select \'within $ARGV[4] of group \"$grpname\_side\"\' -on $ARGV[6]/$ARGV[5]\_$grpname.ndx\n |";
system ($cmd ); 
wait;
}
# calculate frequency of different atoms
my $frame_number=0;
open (INDEX, "$ARGV[6]/$ARGV[5]\_$grpname.ndx");
while (<INDEX>){
if( m/within/){
$frame_number++;
next;
}else{
my @temp=split();
for (my $i=0;$i<=$#temp; $i++){
$atmfrqnc{"$temp[$i]"}++;
}}}
close INDEX;
#Convert atom frequency to occupancy
my @temp=keys %atmfrqnc;
@temp= sort { $a <=> $b }@temp;
my $grpnm;
open (TRNST, "$ARGV[6]/$ARGV[6]\_aminoacids.ndx");
while (<TRNST>){
if ( m/main/ || m/side/ ){
$grpnm = (split)[1];
}else{
my @aaa = split();
if ($#aaa<2&&$grpnm=~/TRP/){
next;
undef @aaa;
}
#TRP contains more atoms
if ($grpnm=~/TRP/){
my $next=$aaa[-1]+1;
push @aaa, $next;
}
if ($temp[0]<=$aaa[-1]&&$aaa[0]<=$temp[-1]){
for (my $i=0;$i<=$#aaa; $i++){
$grpfrqnc{$grpnm}+=$atmfrqnc{$aaa[$i]};
}
$grpfrqnc{$grpnm}/=scalar(@aaa);
$grpfrqnc{$grpnm}/=$frame_number;
}}}
my $exp_enrg=0;
my $obs_enrg=0;
my $difference=0;
my @neighbours= keys %grpfrqnc;
# Remove first and last amino acids from calculation
@neighbours = grep {$_ ne "$firstaa\_side"} @neighbours;
@neighbours = grep {$_ ne "$lastaa\_side"} @neighbours;
@neighbours = grep {$_ ne "$firstaa\_main"} @neighbours;
@neighbours = grep {$_ ne "$lastaa\_main"} @neighbours;
# Remove ignored from calculation
foreach my $i (@ignore){
@neighbours = grep !/$i/, @neighbours;
}
# Calculate $exp_enrg and $obs_enrg
open (REPORT, ">>$out");
print REPORT "$grpname \t ";
for (my $i=0;$i<=$#neighbours; $i++){
my $occupancy=$grpfrqnc{"$neighbours[$i]"};
if ($occupancy){
print REPORT "$neighbours[$i]\($occupancy\)";
}
my $aat=$neighbours[$i];
if ($neighbours[$i]=~/side/){
$aat =~ s/\d//g;
$exp_enrg+= $occupancy*$expected{$aat};
}else{
$aat =~ s/\_main//;
if (grep {$_ eq "-ignoress"} @ARGV){
$exp_enrg+= $occupancy*$expected{average};
}elsif (grep {$_ eq $aat} @helix){
$exp_enrg+= $occupancy*$expected{helix};
}elsif (grep {$_ eq $aat} @sheet){
$exp_enrg+= $occupancy*$expected{sheet};
}elsif (grep {$_ eq $aat} @coil){
$exp_enrg+= $occupancy*$expected{coil};
}else{
die "Kindly assign the secondary structure for $aat";
}}
$obs_enrg+= $occupancy*$actual{"$neighbours[$i]"};
}
$difference=$obs_enrg-$exp_enrg;
print REPORT "\t $exp_enrg \t $obs_enrg \t $difference \n";
close REPORT;
undef %atmfrqnc;
undef %grpfrqnc;
undef $grpname;
}
close INPUT;
