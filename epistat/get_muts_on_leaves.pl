#!/usr/bin/env perl
#For each site in the provided site pairs the script finds all mutations in the tree and a list of leafs carrying each mutation.
#Usage: options <xparr_fn> 
#options:
#	-o STR - the name of output file
#	-m 1|2 - the xparr slot where mutations mapped on the tree branches are
#		1 - background
#		2 - foreground
#	[-p] STR - the site pairs file name. Processes sites from the specified list of pairs
#	[-s] STR - the list of sites file name. Processes sites from the specified list.

use strict;
use Getopt::Std;
use Bio::Phylo::IO;
use IO::File;

my %args;
if(!getopts('o:p:m:s:',\%args)){
	die "\nError in option string!";
}

my $f_sites=0;
my $site_list_fn;
my $site_pairs_fn;

$site_pairs_fn=$args{p} if(defined $args{p});
$site_list_fn=$args{s} if(defined $args{s});
my $out_fn=$args{o} if(defined $args{o});
my $ofh;
if(defined $out_fn){
	$ofh= IO::File->new(">$out_fn") or die "\nUnable to open output file $out_fn: $!";
}else{
	$ofh=*STDOUT;
}
my $mut_col_idx;
if(defined $args{m}){
	if($args{m}=~/^([1|2])$/){
		if($args{m}==1){
			$mut_col_idx=5;
		}else{
			$mut_col_idx=4;
		}
	}else{
		die "\nWrong parameter type for the '-m' option, the uint is expected!";
	}
}else{
	die "\nUnable to identify xparr mutations slot, use '-m' option to specify!";
}
my %sites;
if(defined $site_pairs_fn){
	open INPF,"<$site_pairs_fn" or die "\nUnable to open input file $site_pairs_fn!";
	while(<INPF>){
		chomp;
		if(/^(\d+)\t(\d+)/){
			$sites{$1}={} unless defined $sites{$1};
			$sites{$2}={} unless defined $sites{$2};
		}
	}
	$f_sites=1;
	close INPF;
}
if(defined $site_list_fn){
	open INPF,"<$site_list_fn" or die "\nUnable to open input file $site_list_fn!";
	while(<INPF>){
		chomp;
		if(/^(\d+)/){
			$sites{$1}={} unless defined $sites{$1};
		}
	}
	$f_sites=1;
	close INPF;
}

my $xpar_fn=$ARGV[0];
my $tree;

$tree=Bio::Phylo::IO->parse(
	-file => $xpar_fn,
	-format => 'adjacency');
#!!!The return value of Adjacency parser not corresponded to interface of the Bio::Phylo::IO->parse
$tree=$tree->[0]->first;
open INPF,"<$ARGV[0]" or die "\nUnable to open input file $ARGV[0]!";
my $header;
my %subst_map;
while(<INPF>){
	chomp;
	if(/\S+/){
		if(defined $header){
			my @line=split "\t",$_,-1;
			my $name=$line[0];
			if($line[$mut_col_idx] ne ""){
				my @muts=split ";",$line[$mut_col_idx];
				foreach my $mut(@muts){
					if($mut=~/([A-Z])(\d+)([A-Z])/){
						my $site_id=$2;
						unless($f_sites){
							$sites{$site_id}={} unless defined $sites{$site_id};
						}
						if(defined $sites{$site_id}){
							$subst_map{$name}={} unless defined $subst_map{$name};
							$subst_map{$name}->{$site_id}=[($1,$3)];
							$sites{$site_id}->{$name}=1
						}
					}
				}
			}
		}else{
			$header=$_;
		}
	}
}
close INPF;

sub print_leafs_with_mutation{
	my ($fh,$branch_name,$site_id,$ra_alleles,$ra_leafs)=@_;
	print $fh "\nBranch=\"$branch_name\"\tmutation=\"".$ra_alleles->[0].$site_id.$ra_alleles->[1]."\"";
	foreach my $name(@{$ra_leafs}){
		print $fh "\n$name";
	}
	print $fh "\n";
}

sub delete_descendent_leafs{
	my ($node,$rh_recent_mutations,$rh_descendent_leafs)=@_;
	return if $node->is_terminal;
	my $name=$node->get_name;
	foreach my $mbranch_name(@{$rh_recent_mutations->{$name}}){
		foreach my $term(keys %{$rh_descendent_leafs->{$mbranch_name}}){
			delete $rh_descendent_leafs->{$mbranch_name}->{$term};
		}
		$rh_descendent_leafs->{$mbranch_name}=undef;
		delete $rh_descendent_leafs->{$mbranch_name};
	}
}

#calculate node genotypes
foreach my $site_id(keys %sites){
	print STDERR "\n$site_id";
	my %recent_mutations;
	my %descendent_leafs;
	$tree->visit_depth_first(
		-in   => sub {
			my $node=shift;
			my $name=$node->get_name;
			$descendent_leafs{$name}={};
			unless($node->is_terminal()){
				$recent_mutations{$name}=[];
				foreach my $chnode(@{$node->get_children}){
					my $chname=$chnode->get_name;
					foreach my $leaf(keys %{$descendent_leafs{$chname}}){
						$descendent_leafs{$name}->{$leaf}=1;
					}
					unless($chnode->is_terminal){
						unless(defined($subst_map{$chname})&&defined($subst_map{$chname}->{$site_id})){
							push @{$recent_mutations{$name}},@{$recent_mutations{$chname}};
						}
					}
					if(defined $subst_map{$chname}->{$site_id}){
						push @{$recent_mutations{$name}},$chname;
						die "\nUnexpected XPARR parser error!" unless defined $sites{$site_id}->{$chname};
						#get list of descendent terminals
						my @desc_terms;
						if($chnode->is_terminal){
							push @desc_terms,$chname;
						}else{
							foreach my $term(keys %{$descendent_leafs{$chname}}){
								my $t=1;
								foreach my $mbranch_name(@{$recent_mutations{$chname}}){
									die "\nUnable to find leaves for $mbranch_name!" unless defined $descendent_leafs{$mbranch_name};
									if(defined $descendent_leafs{$mbranch_name}->{$term}){
										$t=0;
										last;
									}
								}
								push @desc_terms,$term if $t;
							}
						}
						print_leafs_with_mutation($ofh,$chname,$site_id,$subst_map{$chname}->{$site_id},\@desc_terms);
						delete_descendent_leafs($chnode,\%recent_mutations,\%descendent_leafs);
					}
				}
			}else{
				$descendent_leafs{$name}->{$name}=1;
			}
		},
		-post   => sub {
			my $node=shift;
			my $name=$node->get_name;
			#free allocated memory
			unless($node->is_terminal){
				foreach my $chnode(@{$node->get_children}){
					my $chname=$chnode->get_name;
					unless(defined $sites{$site_id}->{$chname}){
						foreach my $term(keys %{$descendent_leafs{$chname}}){
							delete $descendent_leafs{$chname}->{$term};
						}
						$descendent_leafs{$chname}=undef;
						delete $descendent_leafs{$chname};
					}
					$recent_mutations{$chname}=undef;
					delete $recent_mutations{$chname};
				}
			}
		}
	);
}
