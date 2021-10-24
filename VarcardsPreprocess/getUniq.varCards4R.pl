#!/usr/bin/perl -w
use strict;

#what this script does is to remove the duplicate entries in Varcard annotation output
my $in = shift;
open O, ">$in.uniq";
my %met;
open I, $in;
while(<I>){
	chomp;
	my @F = split /\t/;
	$F[0] = "$F[0]_$F[1]_$F[2]";
	$_ = join "\t", @F;
	if(!exists $met{$F[0]}){
		print O "$_\n"; 
	}
	$met{$F[0]} = $.;
}
close I;
close O;

for my $k(keys %met){
	#print ",$k,$met{$k},\n";
}
