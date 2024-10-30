#!/usr/bin/env perl
#This script selects subset of sites in xparr file
#Usage: options <xparr_fn>
#options:
#	-f <site_id_list_fn> - Transfer selected sites only
use strict;
use Getopt::Std;
use Bio::Phylo::IO;

my %args;
if(!getopts('f:',\%args)){
	die "\nError in option string!";
}
my $gene_idx=4;
my $site_filter_fn=$args{f};
my %allowed_sites;
if(defined $site_filter_fn){
	open INPF, "<$site_filter_fn" or die "\nUnable to open input file $site_filter_fn!";
	while(<INPF>){
		chomp;
		s/^\s+//;
		s/\s+$//;
		if(/\S+/){
			if(/(\d+)/){
				my $id=$1;
				$allowed_sites{$id}=1;
			}
		}
	}
	close INPF;
}else{
	die "\nThe set of allowed sites is not specified! Use '-f' option.";
}
my $xpar_fn=$ARGV[0];

open INPF,"<$ARGV[0]" or die "\nUnable to open input file $ARGV[0]!";
my $header;
while(<INPF>){
	chomp;
	if(/\S+/){
		if(defined $header){
			my @line=split "\t",$_,-1;
			my $name=$line[0];
			my $mut_str=$line[$gene_idx];
			if($mut_str ne ""){
				my @mutations=split ";",$mut_str;
				my @allowed;
				foreach my $mut(@mutations){
					if($mut=~/(\d+)/){
						push @allowed, $mut if(defined $allowed_sites{$1});
					}
				}
				$mut_str=join ";",@allowed;
			}
			my $str=join("\t",@line[0 .. 3])."\t".$mut_str."\t".$mut_str;
			if($#line>5){
				$str.="\t".join("\t",@line[6 .. $#line]);
			}
			print "\n$str";
		}else{
			$header=$_;
			print $header;
		}
	}
}
close INPF;