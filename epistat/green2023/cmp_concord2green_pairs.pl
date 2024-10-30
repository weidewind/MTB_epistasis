#!/usr/bin/env perl
#This script compares the set of concordant site pairs with site pairs from the Green 2023 doi: 10.1093/molbev/msad131
#The script expects the list of concordantly evolving site pairs with the site annotation
#Usage: options <snp_analyzed> <sequential_pairs> <simultaneous_pairs>
#options:
#	-o STR - the name of output file
#	-p STR - the name of a file with a list of concordant pairs
#	-a STR - the name of annotation file in the gff format
#	-g uint,uint - numbers of columns containing gene names for sites in a pairs

use strict;
use Getopt::Std;
use Bio::Phylo::IO;
use IO::File;
use Class::Struct;

struct MutationInfo => {
	gene_id => '$',
	gene_name => '$',
	gene_site => '$',
	position => '$',
	alleles => '@',
};

struct GeneInfo => {
	type => '$',
	strand => '$',
	locus => '@',
};

struct MutPairInfo => {
	pos1 => '$',
	pos2 => '$',
	rank => '$',
	pvalue => '$',
};

my %args;
if(!getopts('o:g:p:a:',\%args)){
	die "\nError in option string!";
}

my $gff_fn;
my $site_pairs_fn;
if(defined $args{p}){
	$site_pairs_fn=$args{p};
}else{
	die "\nThe epistat site pairs was not specified! Use -p option!";
}
my @column_ids;
if(defined $args{g}){
	if($args{g}=~/(\d+),(\d+)/){
		@column_ids=($1,$2);
		$column_ids[0]--;
		$column_ids[1]--;
	}else{
		die "\nWrong value for -g option, 'uint,uint' - is expected!";
	}
}else{
	die "\nThe column identifiers are not specified! Use -g option!";
}
if(defined $args{a}){
	$gff_fn=$args{a};
}else{
	die "\nThe annotation file is required! Use -a option to specify appropriate GFF file!";
}
my $out_fn=$args{o};
my $ofh;
if(defined $out_fn){
	$ofh= IO::File->new(">$out_fn") or die "\nUnable to open output file $out_fn: $!";
}else{
	$ofh=*STDOUT;
}

my @mutations;
my %gene_id2name;
my %gene_name2id;
my %gene2mutation_id;
my %position2mutation_id;
my $mid=0;
open INPF,"<$ARGV[0]" or die "\nUnable to open input file $ARGV[0]!";
while(<INPF>){
	chomp;
	if(/^\d+\t/){
		my @line=split "\t",$_,-1;
		my $si=MutationInfo->new();
		my $pos=$line[1];
		$si->position($pos);
		die "\nError: Unable to parse genome position for a mutation in $ARGV[0]:\n$_!" unless defined $pos;
		$si->gene_id($line[2]);
		my $gname=$line[3];
		$si->gene_name($gname);
		my $gsite;
		if($line[4]=~/([A-Z*])(\d+)([A-Z*])/){
			$gsite=$2;
			$si->gene_site($gsite);
			$si->alleles([$1,$3]);
		}elsif($line[4] ne "None"){
			print STDERR "\nUnable to parse mutation in $ARGV[0]:\n$_";
		}
		$mutations[$mid]=$si;
		unless($si->gene_id eq "intergenic"){
			$gene_id2name{$si->gene_id}=$gname;
		}
		unless($gname eq "None"){
			$gene_name2id{$gname}=$si->gene_id;
			if(defined $gsite){
				$gene2mutation_id{$gname}={} unless defined $gene2mutation_id{$gname};
				$gene2mutation_id{$gname}->{$gsite}=[] unless defined $gene2mutation_id{$gname}->{$gsite};
				push @{$gene2mutation_id{$gname}->{$gsite}},$mid;
			}
		}
		$position2mutation_id{$pos}=[] unless defined $position2mutation_id{$pos};
		push @{$position2mutation_id{$pos}},$mid;
		$mid++;
	}
}
close INPF;

my %genes;
open INPF, "<$gff_fn" or die "\nUnable to open the input file $gff_fn!";
while(<INPF>){
	chomp;
	if(/\t/){
		my @line=split "\t",$_,-1;
		my $str=$line[8];
		if($str=~/Locus=([^;]+);/){
			my $gene_id=$1;
			if(exists $gene_id2name{$gene_id}){
				my $gi=GeneInfo->new();
				$gi->type($line[2]);
				$gi->strand($line[6]);
				$gi->locus([$line[3],$line[4]]);
				$genes{$gene_id}=$gi;
			}
		}else{
			die "\nUnable to parse line in the GFF file $gff_fn\n$_!";
		}
	}
}
close INPF;
		
sub get_mutations_ids{
	my ($rh_gname2mut,$rh_pos2mut,$rh_gname2gid,$rh_gid2locus,$gname,$site,$mtype)=@_;
	return () unless $mtype=~/^snp/;
	if($mtype eq "snp..aac"){
		my @mids;
		my $rh_sites;
		if(defined $rh_gname2mut->{$gname}){
			$rh_sites=$rh_gname2mut->{$gname};
		}else{
			my $gene_id=$rh_gname2gid->{$gname};
			$rh_sites=$rh_gname2mut->{$gene_id} if defined $gene_id;
		}
		return () unless defined $rh_sites;
		if(defined $rh_sites->{$site}){
			return @{$rh_sites->{$site}};
		}
		return ();
	}
	#nucleotide snp
	if($gname eq "-"){
		if(defined $rh_pos2mut->{$site}){
			return @{$rh_pos2mut->{$site}};
		}
		return ();
	}elsif($site<0){
		my $gene_id=$rh_gname2gid->{$gname};
		if(defined $gene_id){
			my $ginf=$rh_gid2locus->{$gene_id};
			if(defined $ginf){
				my $pos;
				if($ginf->strand eq "+"){
					$pos=$ginf->locus->[0]+$site;
				}else{
					$pos=$ginf->locus->[1]-$site;
				}
				if(defined $rh_pos2mut->{$pos}){
					return @{$rh_pos2mut->{$pos}};
				}
				return ();
			}
		}
		return ();
	}
	return ();
}

my @sequential_pair_info;
my %sequential_pairs;
open INPF,"<$ARGV[1]" or die "\nUnable to open input file $ARGV[1]!";
my $nseq=0;
while(<INPF>){
	chomp;
	if(/^\d+\t/){
		my @line=split "\t",$_,-1;
		my @si=@{$position2mutation_id{$line[0]}};
		my @sj=@{$position2mutation_id{$line[2]}};
		die "\nUnable convert genome positions for the mutation pair from $ARGV[1]:\n$_" unless @si>0&&@sj>0;
		my $pinf=MutPairInfo->new();
		$pinf->pos1($line[0]);
		$pinf->pos1($line[2]);
		$pinf->pvalue($line[12]);
		$pinf->rank($nseq);
		push @sequential_pair_info,$pinf;
		foreach my $i(@si){
			foreach my $j(@sj){
				$sequential_pairs{"$i,$j"}=[] unless defined $sequential_pairs{"$i,$j"};
				push @{$sequential_pairs{"$i,$j"}},$nseq;
			}
		}
		$nseq++;
	}
}
close INPF;

my @simultaneous_pair_info;
my %simultaneous_pairs;
open INPF,"<$ARGV[2]" or die "\nUnable to open input file $ARGV[2]!";
my $nsim=0;
while(<INPF>){
	chomp;
	if(/^\d+\t/){
		my @line=split "\t",$_,-1;
		my @si=@{$position2mutation_id{$line[0]}};
		my @sj=@{$position2mutation_id{$line[2]}};
		die "\nUnable convert genome positions for the mutation pair from $ARGV[2]:\n$_" unless @si>0&&@sj>0;
		my $pinf=MutPairInfo->new();
		$pinf->pos1($line[0]);
		$pinf->pos1($line[2]);
		$pinf->pvalue($line[12]);
		$pinf->rank($nsim);
		push @simultaneous_pair_info,$pinf;
		foreach my $i(@si){
			foreach my $j(@sj){
				$simultaneous_pairs{"$i,$j"}=[] unless defined $simultaneous_pairs{"$i,$j"};
				push @{$simultaneous_pairs{"$i,$j"}},$nsim;
			}
		}
		$nsim++;
	}
}
close INPF;
open INPF, "<$site_pairs_fn" or die "\nUnable to open input file $site_pairs_fn!";
my $header;
while(<INPF>){
	chomp;
	if(/\t/){
		my @line=split "\t",$_,-1;
		if(defined $header){
			my ($gene1,$site1,$mtype1)=($line[$column_ids[0]],$line[$column_ids[0]+1],$line[$column_ids[0]+2]);
			my ($gene2,$site2,$mtype2)=($line[$column_ids[1]],$line[$column_ids[1]+1],$line[$column_ids[1]+2]);
			my @site1_ids=get_mutations_ids(\%gene2mutation_id,\%position2mutation_id,\%gene_name2id,\%genes,$gene1,$site1,$mtype1);
			my @site2_ids=get_mutations_ids(\%gene2mutation_id,\%position2mutation_id,\%gene_name2id,\%genes,$gene2,$site2,$mtype2);
			my @seqpairs;
			if(@site1_ids>0&&@site2_ids>0){
				foreach my $sid1(@site1_ids){
					foreach my $sid2(@site2_ids){
						push @seqpairs,@{$sequential_pairs{"$sid1,$sid2"}} if defined $sequential_pairs{"$sid1,$sid2"};
						push @seqpairs,@{$sequential_pairs{"$sid2,$sid1"}} if defined $sequential_pairs{"$sid2,$sid1"};
					}
				}
			}
			my @simpairs;
			if(@site1_ids>0&&@site2_ids>0){
				foreach my $sid1(@site1_ids){
					foreach my $sid2(@site2_ids){
						push @simpairs,@{$simultaneous_pairs{"$sid1,$sid2"}} if defined $simultaneous_pairs{"$sid1,$sid2"};
						push @simpairs,@{$simultaneous_pairs{"$sid2,$sid1"}} if defined $simultaneous_pairs{"$sid2,$sid1"};
					}
				}
			}
			print $ofh "$_\t";
			print $ofh @site1_ids>0?"+":"-";
			print $ofh @site2_ids>0?"+":"-";
			print $ofh "\t";
			if(@seqpairs){
				my @ranks;
				my @pvals;
				foreach my $pid(@seqpairs){
					my $pinf=$sequential_pair_info[$pid];
					push @ranks,$pinf->rank;
					push @pvals,$pinf->pvalue;
				}
				my $str=join ";",@ranks;
				print $ofh "$str\t";
				$str=join ";",@pvals;
				print $ofh "$str\t";
			}else{
				print $ofh "\t\t";
			}
			if(@simpairs){
				my @ranks;
				my @pvals;
				foreach my $pid(@simpairs){
					my $pinf=$simultaneous_pair_info[$pid];
					push @ranks,$pinf->rank;
					push @pvals,$pinf->pvalue;
				}
				my $str=join ";",@ranks;
				print $ofh "$str\t";
				$str=join ";",@pvals;
				print $ofh "$str\n";
			}else{
				print $ofh "\t\n";
			}
		}else{
			$header=$_."\tf_sites\tseq_ranks\tseq_pvals\tsim_ranks\tsim_pvals";
			print $ofh "$header\n";
		}
	}
}
close INPF;
