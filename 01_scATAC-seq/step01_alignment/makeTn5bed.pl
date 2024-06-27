#!/usr/bin/perl
use strict;
use warnings;
use perlSAM qw(:all);
use Sort::Naturally;

my %flag = (
	0  => 'read_paired',
	1  => 'read_properly',
	2  => 'read_umapped',
	3  => 'mate_unmapped',
	4  => 'read_reverse',
	5  => 'mate_reverse',
	6  => 'first_pair',
	7  => 'second_pair',
	8  => 'secondary',
	9  => 'fail_qc',
	10 => 'duplicate',
	11 => 'sup_align'
);
	
my $counts = 0;
open F, "samtools view $ARGV[0] | " or die;
while(<F>){
	chomp;
	$counts++;
	if(($counts % 1000000) == 0){
		print STDERR " - iterated over $counts reads ... \n";
	}
	my @col = split("\t",$_);
	my $bc;
	foreach(@col[11..$#col]){
		if($_ =~ /CB:Z:/){
			$bc = $_;
			last;
		}
	}
	my $bin = sprintf ("%.12b", $col[1]);
	my @bins = split("",$bin);
	my @flags;

	# load flags
	for(my $i = 0; $i < @bins; $i++){
		if($bins[$i] == 1){
			push(@flags, $flag{$i});
		}
	}

	# get start and end pos of read
	my $chr = $col[2];
	my $pos1 = $col[3];
	my $cigar = parseCIGAR($col[5]);
	my %cig = %$cigar;
	my $netdif = 1;
	my @act = nsort keys %cig;
	for (my $i = 0; $i < @act; $i++){
		my @action = split("_",$act[$i]);
		my $times = $cig{$act[$i]};
		if($action[1] eq 'M' || $action[1] eq 'S'){
			$netdif = $netdif + $times;
		}elsif($action[1] eq 'D'){
			$netdif = $netdif + $times;
		}
	}
	my $pos2 = $netdif + $pos1;
	if(grep {$_ =~ /read_reverse/} @flags){
		my $end = $pos2 - 4;
		my $start = $end - 1;
		print "$chr\t$start\t$end\t$bc\t-\n";
	}else{
		my $start = $pos1 + 5;
		my $end = $start + 1;
		print "$chr\t$start\t$end\t$bc\t+\n";
	}
}
close F;

# subroutines

