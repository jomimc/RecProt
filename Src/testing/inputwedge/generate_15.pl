#!/usr/bin/perl
use strict;
{
    open(my $IN,"<ligand3.dat")|| die;
    my @ligs=();
    while(my $line=<$IN>){
	chomp($line);
	$line=~s/1/0/g;
	$line=~s/2/1/g;
	$line=~ s/ //g;
	push(@ligs,$line);

    }

    for my $i (0..8191){
	for my $k (0..$#ligs){
	    my $width = 14;

	    printf "%${width}b", $i+8192;

	    print "$ligs[$k]\n";
	}
    }
}
