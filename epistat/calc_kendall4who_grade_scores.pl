#!/usr/bin/env perl
#This script calculates for each pair the number of its sites having WHO grades below 3 (associated with resistance) and then calculates Kendall's correlation with the specified table column.
#Usage: options <table_fn>
#options:
#	-g uint(,uint)* - list of columns containing the WHO grades
#	-s uint - the column with statistics wich is used for calculation of Kendall's correlation
#	[-o] STR - a suffix for output file name
use strict;
use Getopt::Std;
use File::Basename;
my $cor_test_cmd="/home/common/fedonin.gg/devel/epistat.10/stat/cor_test.R";

my %args;
if(!getopts('g:s:o:',\%args)){
	die "\nError in option string!";
}
my $out_fn=$args{o};
my @column_ids;
if(defined $args{g}){
	@column_ids=split ",",$args{g};
	for(my $i=0;$i<@column_ids;$i++){
		$column_ids[$i]--;
	}
}else{
	die "\nUnable to automatically determine columns with WHO grades! Use -g option!";
}
my $stat_col;
if(defined $args{s}){
	if($args{s}=~/^(\d+)/){
		$stat_col=$1-1;
	}
}
die "\nUnable to determine the column with statistics to test!" unless defined $stat_col;
if(defined $out_fn){
	my ($basename,$dir,$ext) = fileparse($ARGV[0],'\.[^\.]*$');
	$out_fn=$dir.$basename.$out_fn;
	open OUTF,">$out_fn" or die "\nUnable to open output file $out_fn!";
}else{
	*OUTF=*STDOUT;
}
my $header_str;
my $stat_col_name;
my @stats;
my @pair_scores;
open INPF,"<$ARGV[0]" or die "\nUnable to open input file $ARGV[0]!";
while(<INPF>){
	chomp;
	if(/\S+/){
		my @line=split "\t",$_,-1;
		if(defined $header_str){
			my $str;
			my @site_score=(0,0);
			my $j=0;
			foreach my $i(@column_ids){
				my $grades=$line[$i];
				my $grade=10;
				if(defined $grades){
					my @drugs=split ";", $grades;
					foreach my $drug(@drugs){
						if($drug=~/\:(\d+)/){
							$grade=$1 if($1<$grade);
						}
					}
					$site_score[$j]=1 if $grade<3; #accociated with resistance to something
				}
				$j++;
			}
			my $pair_score=$site_score[0]+$site_score[1];
			$str=join("\t",@line)."\t$pair_score\n";
			print OUTF $str;
			push @pair_scores,$pair_score;
			push @stats,$line[$stat_col];
		}else{
			$header_str=join("\t",@line)."\tdrug_resist_score\n";
			print OUTF $header_str;
			$stat_col_name=$line[$stat_col];
		}
	}
}
close INPF;

sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
}

sub print_child_termination_status{
	if ($? == -1) {
		print "failed to execute: $!\n";
		exit 1;
	}elsif ($? & 127) {
		printf "\tchild died with signal %d, %s coredump\n",
					($? & 127),  ($? & 128) ? 'with' : 'without';
	}else {
		if($? >> 8!=0){
			printf "\tchild exited with value %d\n", $? >> 8;
		};
	}
}

my $tmp_fn=gen_tempname(6).".tab";
open OPF,">$tmp_fn" or die "\nUnable to open output file $tmp_fn!";
print OPF "drug_resist_score\t$stat_col_name";
for(my $i=0;$i<@stats;$i++){
	print OPF "\n$pair_scores[$i]\t$stats[$i]";
}
close OPF;
my $str=$cor_test_cmd." $tmp_fn kendall";
$str.=">>$out_fn" if defined $out_fn;
print STDERR "\n$str";
system($str);
print_child_termination_status();
unlink $tmp_fn;