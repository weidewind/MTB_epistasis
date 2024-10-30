#!/usr/bin/env perl
#This script aggregates site pairs of highly linked sites into clusters
#Usage: options <branch_pair_table>
#   Note: the first two columns in the <branch_pair_table> are expected to be site IDs.
#   Note: the table have to contain columns in the following order: 
#        <number of mutations in branch pairs in the first site> - #1
#        <number of mutations in branch pairs in the second site>
#        <total number of mutations in the first site>
#        <total number of mutations in the second site>
#   Note: F_min - the ratio of mutation not linked (not on the same branches) with mutations in other site in the site with the minimal number of mutations in a pair
#         F_max - the ratio of mutation not linked (not on the same branches) with mutations in other site in the site with the maximal number of mutations in a pair
#options:
#	-f ufloat - the threshold for F_min statistic, site pairs with F_min below threshold would be aggregated
#	[-F] ufloat - the threshold for F_max statistic, site pairs with F_max upper threshold would be aggregated
#		default=0.0
#	-c uint - the index of column #1 with the number of mutations in branch pairs in the first site.


use strict;
use Getopt::Std;
my %args;
if(!getopts('f:F:c:',\%args)){
	die "\nError in option string!";
}
my $fmin_upper_threshold;
if(defined $args{f}){
	$fmin_upper_threshold=$args{f} if $args{f}=~/^0\.\d+/;
}
die "\nThe F_min threshold is not defined! Use -f option." unless defined $fmin_upper_threshold;
my $fmax_lower_threshold=0;
if(defined $args{F}){
	$fmax_lower_threshold=$args{F} if $args{F}=~/^0\.\d+/;
}
my ($n1_bp_coln,$n2_bp_coln,$n1_total_coln,$n2_total_coln);
if(defined $args{c}){
	if($args{c}=~/^\d+/){
		$n1_bp_coln=$args{c};
		$n1_bp_coln--;
		$n2_bp_coln=$n1_bp_coln+1;
		$n1_total_coln=$n1_bp_coln+2;
		$n2_total_coln=$n1_bp_coln+3;
	}
}
die "\nThe column with the number of branch pair mutations is not specified! Use -c option." unless defined $n1_bp_coln;

my @branch_pairs;
my %clusters;
my %in_cluster;
open INPF, "<$ARGV[0]" or die "\Unable to open input file $ARGV[0]!";
my $header;
my $i=0;
while(<INPF>){
	if(/\S+/){
		if(defined $header){
			my @line=split "\t",$_,-1;
			push @branch_pairs,$_;
			die "\nFirst two columns must be site IDs!" unless ($line[0]=~/^\d+$/)&&($line[1]=~/^\d+$/);
			die "\nThe number of columns in the input table is less than expected" unless $n2_total_coln<@line;
			my ($fmin,$fmax);
			my $n1=$line[$n1_bp_coln];
			my $n2=$line[$n2_bp_coln];
			my $n1t=$line[$n1_total_coln];
			my $n2t=$line[$n2_total_coln];
			my $index_site;
			if($n1t<$n2t){
				$fmin=1.0-$n1/$n1t;
				$fmax=1.0-$n2/$n2t;
				$index_site=$line[1];
			}else{
				$fmin=1.0-$n2/$n2t;
				$fmax=1.0-$n1/$n1t;
				$index_site=$line[0];
			}
			
			if(($fmin<$fmin_upper_threshold)&&($fmax>$fmax_lower_threshold)){
				$clusters{$index_site}=[] unless defined $clusters{$index_site};
				push @{$clusters{$index_site}},$i;
				$in_cluster{$i}=1;
			}
			$i++;
		}else{
			$header=$_;
		}
	}
}
close INPF;
#print clusters
my @site_keys=sort {scalar(@{$clusters{$b}})<=>scalar(@{$clusters{$a}})} keys %clusters;
print "branch_pair_clstr_id\t$header";
my $k=1;
foreach my $index_site(@site_keys){
	foreach my $i(@{$clusters{$index_site}}){
		print "$k\t$branch_pairs[$i]";
	}
	$k++;
}
for(my $i=0;$i<@branch_pairs;$i++){
	unless($in_cluster{$i}){
		print "$k\t$branch_pairs[$i]" ;
	}
}
