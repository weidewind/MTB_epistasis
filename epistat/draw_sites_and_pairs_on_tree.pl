#!/usr/bin/env perl
#This script draws mutations on tree branches
#options:
#	-x <FN> - input file for EpiStat application in XPAR format [Kryazhimsky11]. Contains tree in the adjacent file format
#	-p <epistat_prm> - input file with parameters for running epistat.pl script
#	[-u] unordered site pairs. Print both ordered site pairs on one tree
#	[-s] <FN> - a list of other sites which have to be drawn additionally to each pair
#	[-a] <FN> - a list of all site pairs which is used to find sites with consecutive mutations with sites from the <site_pairs_list>
#		This option is mandatory if the '-s' option has been specified
#	[-c] uint - an id of the "npair_subst" column in the table with all site pairs
#		This option is mandatory if the '-a' option has been specified
#Params:
#<site_pairs_list> - a file with list of site pairs in TAB delimited format. Mutations for each pair are drawn on separate tree

use strict;
use Getopt::Std;
use File::Basename;

my $epistat_cmd="$ENV{EPISTAT_HOME}/stat/draw_site_pairs_on_tree.pl";

my %args;
if(!getopts('x:p:us:a:c:',\%args)){
	die "\nError in option string!";
}

my $xparr_fn=$args{x};
die "\nNo default parameter for -x option!" unless defined $xparr_fn;
my $epistat_prm=$args{p};
die "\nNo default parameter for -p option!" unless defined $epistat_prm;
my $site_pairs_fn=$ARGV[0];
die "\nThe file with a list of site pairs to draw is not accounted!" unless defined $site_pairs_fn;
my $other_sites_fn=$args{s};
my $all_site_pars_fn=$args{a};
die "\nIncomplete information about other sites has been provided! Use both '-a', '-s' and '-c' options." unless defined($other_sites_fn)&&defined($all_site_pars_fn)&&defined($args{c});
my $npairs_cid;
my %other_sites;
if(defined $other_sites_fn){
	$npairs_cid=$args{c} if $args{c}=~/^\d+$/ or die "\nWrong type of the argument of the '-c' option.";
	$npairs_cid--;
	open  INPF, "<$other_sites_fn" or die "\nUnable to open input file $other_sites_fn!";
	while(<INPF>){
		chomp;
		if(/^(\d+)/){
			$other_sites{$1}=1;
		}
	}
	close INPF;
}

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
		}
	}
}

my @tmp_fn;
open INPF, "<$site_pairs_fn" or die "\nUnable to open input file: $site_pairs_fn!";
my @site_pairs;
my %sites_in_pairs;
while(<INPF>){
	my @line=split '\t';
	$line[0]=~s/\s+//;
	$line[1]=~s/\s+//;
	if($line[0]=~/\d+/ || $line[1]=~/\d+/){
		my ($bgs,$tgs)=($line[0],$line[1]);
		$sites_in_pairs{$tgs}=1;
		$sites_in_pairs{$bgs}=1;
		my $rp=[];
		if(defined $args{u}){
			if($bgs>$tgs){
				@{$rp}=($tgs,$bgs);
			}else{
			@{$rp}=($bgs,$tgs);
			}
		}else{
			@{$rp}=($bgs,$tgs);
		}
		push @site_pairs,$rp;
	}
}
close INPF;

my %site2others;
if(defined $all_site_pars_fn){
	open INPF, "<$all_site_pars_fn" or die "\nUnable to open input file $all_site_pars_fn!";
	while(<INPF>){
		chomp;
		if(/^\d+\t\d+\t/){
			my @line=split "\t",$_,-1;
			my @pair=@line[0..1];
			my $other_site;
			my $pair_site;
			if($line[$npairs_cid]){
				for(my $i=0;$i<2;$i++){
					my $j=($i+1)%2;
					if(defined($sites_in_pairs{$pair[$i]})&&defined($other_sites{$pair[$j]})){
						$site2others{$pair[$i]}=[] unless defined $site2others{$pair[$i]};
						push @{$site2others{$pair[$i]}},$pair[$j];
					}
				}
			}
		}
	}
	close INPF;
}

foreach my $sp(@site_pairs){
	my $tmp_fn="pair_".join("_",$sp->[0],$sp->[1]).".tab";
	open OPF, ">$tmp_fn" or die "\nUnable to open output file: $tmp_fn!";
	print OPF join("\t",$sp->[0],$sp->[1]);
	if(defined $all_site_pars_fn){
		my %others;
		for(my $i=0;$i<2;$i++){
			if(defined $site2others{$sp->[$i]}){
				foreach my $site(@{$site2others{$sp->[$i]}}){
					$others{$site}++;
				}
			}
		}
		foreach my $site (keys %others){
			if($others{$site}==2){
				print OPF "\n$site\t".$sp->[0];
				print OPF "\n$site\t".$sp->[1];
			}
		}
	}
	close OPF;
	my $str="$epistat_cmd -x $xparr_fn -p $epistat_prm";
	$str.=" -u" if defined $args{u};
	$str.=" $tmp_fn";
	print STDERR "\n\nPrinting trees!\n$str";
	system($str);
	print_child_termination_status();
}
