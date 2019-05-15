#! /usr/bin/perl
# Programe to creat separate index group and compliment index group for each amino acid main and side chain from .gro file 
# Creates one index group for each amino acids and file containing name of energy groups  
# Usage list.pl input_file_name output_file_name 
use strict;
my %grpnm;
my %atmno;
my @correct;
my @protein;
my @nonprotein;
my @side;
my @sidec;
my @system;
my @main;
my @mainc;
my $lpa;
my $totatm;
my $ln=1;
my $sol=0;
my $an=2;
my $gn=1;
my $us;
#Read structure file (.gro) to find total number of atoms and last protein atom number
open (READ, "$ARGV[0]") || die "Couldn,t open file $ARGV[0]";
while (<READ>){
if ($ln==2){
$totatm=$_;
}
if ( m/SOL\s/){
$lpa=(split)[2];
last;
}
$ln++;
}
open (WRITE, ">$ARGV[1]\_enrggrps.txt") || die "Couldn,t creat file $ARGV[1]";
close WRITE;
#Detecting first and last atom of an aminoacid
open (INPUT, "$ARGV[0]") || die "Couldn,t open file $ARGV[0]";
while (<INPUT> ){
	if ( m/\sO1\s/){
	last; 
	}
	if ( m/\sCA\s/){
	$atmno{$an}= (split)[2];
	$an++;
	$grpnm{$gn}= (split)[0];
	$gn++;
	}
	if ( m/\sO\s/){
	$atmno{$an}= (split)[2];
	$an++;
	}
	}
close INPUT;
#Creating text file containing name of energy groups (amino acid number and 3 digit name) 
for (my $i=1; $i<$gn; $i++){
	open (WRITE, ">>$ARGV[1]\_enrggrps.txt");
	print WRITE "$i \t $grpnm{$i} \n";
	close WRITE;
}
#Generating atom numbers for each index groups
for (my $i=1; $i<$gn+1; $i++){
	@protein=(1 .. $lpa-1);
	@nonprotein=($lpa .. $totatm);
	@side=($atmno{$i*2}+1 .. $atmno{$i*2 +1}-2);
	@sidec=(1 .. $atmno{$i*2}, $atmno{$i*2 +1}-1 .. $totatm);
	if ($grpnm{$i}=~ /PRO/){
		@mainc=(1 .. $atmno{$i*2}-2, $atmno{$i*2}+1 .. $atmno{$i*2 +1}-2, $atmno{$i*2 +1}+1 .. $totatm);
	}else{
		@mainc=(1 .. $atmno{$i*2}-3, $atmno{$i*2}+1 .. $atmno{$i*2 +1}-2, $atmno{$i*2 +1}+1 .. $totatm);
	}
	@system=(1..$totatm);
	if ($grpnm{$i}=~ /PRO/){
		@main=($atmno{$i*2}-1, $atmno{$i*2}, $atmno{$i*2 +1}-1, $atmno{$i*2 +1});
	}else{	
		@main=($atmno{$i*2}-2, $atmno{$i*2}-1, $atmno{$i*2}, $atmno{$i*2 +1}-1, $atmno{$i*2 +1});
	}
#Arrenging atom numbers in a group of 15
	@system=modify(@system);
	@protein=modify(@protein);
	@nonprotein=modify(@nonprotein);
	@side=modify(@side);
	@main=modify(@main);
	@sidec=modify(@sidec);
	@mainc=modify(@mainc);
#Creating index groups 
	open (WRITE, ">$ARGV[1]\_$i.ndx") || die "Couldn,t creat file $ARGV[1]";
	if ($grpnm{$i}=~ /GLY/){
	print WRITE "[ System ] \n @system \n [ Protein ] \n @protein\n [ Non-protein ] \n @nonprotein \n [ $grpnm{$i}\_main ] \n @main \n [ $grpnm{$i}\_main_com ] \n @mainc \n";
	}else{
	print WRITE "[ System ] \n @system \n [ Protein ] \n @protein\n [ Non-protein ] \n @nonprotein \n [ $grpnm{$i}\_main ] \n @main \n [ $grpnm{$i}\_side ] \n @side \n [ $grpnm{$i}\_main_com ] \n @mainc \n [ $grpnm{$i}\_side_com ] \n @sidec \n";
	}
	close WRITE;
#Creating index groups for each amino acid side chain and main chain.
	open (AMNACD, ">>$ARGV[1]\_aminoacids.ndx");
	print AMNACD "[ $grpnm{$i}_main ] \n @main \n [ $grpnm{$i}_side ] \n @side \n  ";
	close AMNACD;
		if ($i==$gn+1){
		open (CMPLT, ">>$ARGV[1]\_aminoacids.ndx");
		print CMPLT "[ System ] \n @system \n [ Protein ] \n @protein\n [ Non-protein ] \n @nonprotein \n"; 
		close CMPLT;
		}
}
#subroutine modify to add new line after 15 numbers
sub modify  {
	my $k=15; 
	@correct=@_;
	if ($#_>14){
	for (my $i=1; $i<=$#_/15; $i++){
		splice(@correct, $k, 0, "\n",); 
		$k=$k+16;
		}
	}
	return @correct;
}

