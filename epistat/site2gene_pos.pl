#!/usr/bin/env perl
#The script converts site id to the gene and site position
#Usage: options <table>
#options:
#	-c uint(,uint)* - the list of indices of columns in the table file containing site IDs
#	-a STR - a name of file with annotations for site IDs
#	[-g] STR - a name of file with genome annotation in gff.v4 format

use strict;
use Getopt::Std;
my %args;
if(!getopts('c:a:g:',\%args)){
	die "\nError in option string!";
}
my $site_annotation_fn=$args{a};
die "\nThe file with site annotations is required. Please specify '-a ' option!" unless defined $args{a};
die "\nUnable to find columns with site IDs. Please specify '-c ' option!" unless defined $args{c};
my $gff_fn=$args{g};
my @column_ids=split(',', $args{c});
for(my $i=0;$i<@column_ids;$i++){
	$column_ids[$i]--;
}
my $table_fn=$ARGV[0];
my %gene_type;
if(defined $gff_fn){
	open INPF, "<$gff_fn" or die "\nUnable to open gff file: $gff_fn!";
	my $i=1;
	while(<INPF>){
		chomp;
		if(/\S+/){
			my @line= split "\t",$_,-1;
			if($line[8]=~/Name=(.+?);/){
				$gene_type{$1}=$line[2];
			}else{
				die "\nWrong gff format in file $gff_fn line $i!";
			}
			$i++;
		}
	}
	close INPF;
}

my %site2gene_pos;
open INPF, "<$site_annotation_fn" or die "\nUnable to open input file $site_annotation_fn!";
while(<INPF>){
	chomp;
	if(/^\d+\t/){
		my @line=split "\t",$_,-1;
		die "\nWrong format of the file with site id to gene conversion table $ARGV[1]!" if @line <3;
		my $str=join "\t",$line[1],$line[2];
		if(defined $gff_fn){
			my $mtype=$line[-1];
			my $stype="";
			if(($mtype eq "snp")||($mtype eq "del")||($mtype eq "ins")){
				if($gene_type{$line[1]} eq "CDS"){
					if($line[2]<0){
						$stype="nuc";
					}else{
						$stype="aac";
					}
				}else{
					$stype="nuc";
				}
			}
			my $n;
			if($mtype eq "del"){
				if(($gene_type{$line[1]} eq "CDS")&&($stype eq "aac")){
					$n=$line[-2];
				}else{
					$n=length($line[-3])-length($line[-2]);
				}
			}elsif($mtype eq "ins"){
				$n=length($line[-2]);
				$n-=length($line[-3]) unless ($gene_type{$line[1]} eq "CDS")&&($stype eq "aac");
			}
			$str.="\t$mtype.$n.$stype";
		}
		$site2gene_pos{$line[0]}=$str;
	}
}
close INPF;
open INPF, "<$ARGV[0]" or die "\nUnable to open input file $ARGV[0]!";
my $header_str;
while(<INPF>){
	chomp;
	if(/\S+/){
		my @line=split "\t",$_,-1;
		my $begin=0;
		if(defined($header_str)){
			my $str;
			foreach my $i(@column_ids){
				my $site=$line[$i];
				my $tmp=$site2gene_pos{$site};
				die "\nUnknown site id found in the file $ARGV[0]!" unless defined $tmp;
				$str.=join("\t",@line[$begin .. $i])."\t$tmp\t";
				$begin=$i+1;
			}
			if($begin<=$#line){
				$str.=join("\t",@line[$begin .. $#line]);
			}else{
				$str=~s/\t$//;
			}
			$str.="\n";
			print $str;
		}else{
			foreach my $i(@column_ids){
				$header_str.=join("\t",@line[$begin .. $i])."\tgene\tpos\t";
				$header_str.="mtype\t" if(defined $gff_fn);
				$begin=$i+1;
			}
			if($begin<=$#line){
				$header_str.=join("\t",@line[$begin .. $#line]);
			}else{
				$header_str=~s/\t$//;
			}
			print "$header_str\n";
		}
	}
}
close INPF;