#!/usr/bin/env perl
#This script assignes a status for an annotated site in the MTB genome according to the presence or absence of mutations conferring drug resistance in this site in the WHO catalog 2021.
#Usage:	options <table_fn> <catalog_genome_indices_fn>
#options:
#	-c uint(,uint)* - the list of indices of columns in the table file containing annotations of MTB genome sites
#		The site annotation in the table expected within two adjacent colums: 'gene_name' and 'gene position'.
#		The indices of the first column ('gene_name') in each pair is expected.

use strict;
use Getopt::Std;
use Class::Struct;

struct Mutation => {
	variant => '$',
	annotation => '$',
	position => '$',
	effect => '$',
	drug_conf_grades => '@'
};

my %args;
if(!getopts('c:',\%args)){
	die "\nError in option string!";
}
die "\nUnable to find columns to assign catalog status. Please specify '-c ' option!" unless defined $args{c};
my @column_ids=split(',', $args{c});
for(my $i=0;$i<@column_ids;$i++){
	$column_ids[$i]--;
}
my $table_fn=$ARGV[0];
my $catalog_fn=$ARGV[1];
my %catalog;
my $header_str;
my %catalog_header;
my @drugs=qw(INH RIF STM EMB PZA LEV MXF AMI CAP KAN ETH BDQ LZD CFZ DLM);
my $ndrugs=@drugs;
my $n=0;
open INPF, "<$catalog_fn" or die "Unable to open catalog file $catalog_fn!";
while(<INPF>){
	chomp;
	$n++;
	if(/\S+/){
		my @line=split "\t",$_,-1;
		if(defined $header_str){
			my $gene_name=$line[$catalog_header{'gene_name'}];
			my $variant=$line[$catalog_header{'variant'}];
			$catalog{$gene_name}={(gene_product => {}, upstream=>{}, lof => {})} unless defined $catalog{$gene_name};
			my $mut=Mutation->new();
			my $effect=$line[$catalog_header{'effect'}];
			$mut->effect($effect);
			$mut->variant($variant);
			my $ann_gene=$line[$catalog_header{'gene_annotation'}];
			my $pos;
			my $annotation;
			if($effect=~/upstream_gene_variant/){
				$annotation=$line[$catalog_header{'nuc_annotation'}];
				if($annotation=~/^[cn]\.(-\d+)/){
					$pos=$1;
					$mut->position($pos);
					$mut->annotation($annotation);
					$catalog{$gene_name}->{upstream}->{$pos}=[] unless defined $catalog{$gene_name}->{upstream}->{$pos};
					push @{$catalog{$gene_name}->{upstream}->{$pos}},$mut;
				}
			}elsif($effect=~/non_coding_transcript_exon_variant/){
				$annotation=$line[$catalog_header{'nuc_annotation'}];
				if(($ann_gene eq $gene_name)&&($annotation=~/^[cn]\.(-?\d+)/)){
					$pos=$1;
					$mut->position($pos);
					$mut->annotation($annotation);
					if($pos<0){
						$catalog{$gene_name}->{upstream}->{$pos}=[] unless defined $catalog{$gene_name}->{upstream}->{$pos};
						push @{$catalog{$gene_name}->{upstream}->{$pos}},$mut;
					}else{
						$catalog{$gene_name}->{gene_product}->{$pos}=[] unless defined $catalog{$gene_name}->{gene_product}->{$pos};
						push @{$catalog{$gene_name}->{gene_product}->{$pos}},$mut;
					}
				}
			}else{
				$annotation=$line[$catalog_header{'prot_annotation'}];
				if($annotation=~/^p\.[a-z]+(\d+)/i){
					$pos=$1;
					$mut->position($pos);
					$mut->annotation($annotation);
				}
				if(defined $pos){
					if($effect=~/frameshift_variant/||
						$effect=~/stop_gained/||
						$effect=~/stop_lost/||
						$effect=~/start_lost/){
						$catalog{$gene_name}->{lof}->{$pos}=[] unless defined $catalog{$gene_name}->{lof}->{$pos};
						push @{$catalog{$gene_name}->{lof}->{$pos}},$mut;
					}
					if($effect=~/missense_variant/||
						$effect=~/_inframe/){
						$catalog{$gene_name}->{gene_product}->{$pos}=[] unless defined $catalog{$gene_name}->{gene_product}->{$pos};
						push @{$catalog{$gene_name}->{gene_product}->{$pos}},$mut;
					}
				}
			}
			if(defined $pos){
				#read drug confidence grading
				for(my $i=0;$i<$ndrugs;$i++){
					my $dr=$drugs[$i];
					my $str=$line[$catalog_header{$dr."_Conf_Grade"}];
					if($str=~/^([1-5])\)/){
						$mut->drug_conf_grades->[$i]=$1;
					}else{
						$mut->drug_conf_grades->[$i]=0;
					}
				}
			}else{
				warn "Parser error line $n: variant=$variant, effect=$effect, annotation=$annotation, annotation_gene=$ann_gene";
			}
		}else{
			$header_str=$_;
			die "\nError: wrong format of the WHO catalog file $catalog_fn!" unless $header_str=~/^gene_name\tgene_locus\tvariant\tcodon_number\tgenome_index/;
			for(my $i=0;$i<@line;$i++){
				if($line[$i] eq "gene_name"||
				$line[$i] eq  "variant"||
				$line[$i]=~/_Conf_Grade/){
					$catalog_header{$line[$i]}=$i;
				}elsif($line[$i] eq "final_annotation.Effect"){
					$catalog_header{'effect'}=$i;
				}elsif($line[$i] eq "final_annotation.TentativeHGVSNucleotidicAnnotation"){
					$catalog_header{'nuc_annotation'}=$i;
				}elsif($line[$i] eq "final_annotation.TentativeHGVSProteicAnnotation"){
					$catalog_header{'prot_annotation'}=$i;
				}elsif($line[$i] eq "final_annotation.Gene"){
					$catalog_header{'gene_annotation'}=$i;
				}
			}
		}
	}
}
close INPF;

sub _set_grades{
	my ($ra_mutations,$ra_grades)=@_;
	my @grades=();
	@grades=@{$ra_grades} if defined $ra_grades;
	if(defined $ra_mutations){
		foreach my $mut(@{$ra_mutations}){
			for(my $i=0;$i<$ndrugs;$i++){
				my $gr=$mut->drug_conf_grades->[$i];
				if($gr){
					if(!defined($grades[$i])){
						$grades[$i]=$gr;
					}elsif($grades[$i]==3){
						$grades[$i]=$gr if $gr!=3;
					}else{
						if($gr<3){
							$grades[$i]=$gr if $gr<$grades[$i];
						}elsif($gr>3){
							if($grades[$i]>3){
								$grades[$i]=$gr if $gr>$grades[$i];
							}
						}
					}
				}
			}
		}
	}
	return @grades;
}

sub calc_site_conf_grades{
	my ($gene_name,$pos)=@_;
	return () if(!defined($gene_name)||$gene_name eq "-");
	return () unless defined $catalog{$gene_name};
	my @grades;
	my $ra_muts;
	if($pos eq "broken"){
		foreach my $pos(keys %{$catalog{$gene_name}->{lof}}){
			$ra_muts=$catalog{$gene_name}->{lof}->{$pos};
			@grades=_set_grades($ra_muts,\@grades);
		}
	}else{
		if($pos<0){
			$ra_muts=$catalog{$gene_name}->{upstream}->{$pos};
		}else{
			$ra_muts=$catalog{$gene_name}->{gene_product}->{$pos};
		}
		@grades=_set_grades($ra_muts);
	}
	return @grades;
}

#start script
$header_str=undef;
open INPF,"<$table_fn" or die "\nUnable to open input file $table_fn!";
while(<INPF>){
	chomp;
	if(/\S+/){
		my @line=split "\t",$_,-1;
		my $begin=0;
		if(defined $header_str){
			my $str;
			foreach my $i(@column_ids){
				my @grades=calc_site_conf_grades($line[$i],$line[$i+1]);
				my $grade_str;
				if(@grades){
					for(my $j=0;$j<$ndrugs;$j++){
						$grade_str.=$drugs[$j].":".$grades[$j].";" if defined $grades[$j];
					}
					$grade_str=~s/;$//;
				}
				$str.=join("\t",@line[$begin .. $i+1])."\t$grade_str\t";
				$begin=$i+2;
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
				$header_str.=join("\t",@line[$begin .. $i+1])."\tWHO_grade\t";
				$begin=$i+2;
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