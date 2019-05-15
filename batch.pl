#! /usr/bin/perl
# Programe to run serise of command listed in a text file 
# Usage batch.pl input.txt
use strict;
use warnings;
open (READ, "$ARGV[0]") || die "Couldn,t open file $ARGV[0]";
while (<READ>){
my $cmd=$_;
system ($cmd );
wait;
}
close READ;
print "Completed";
