#!/usr/bin/perl

use strict;
use FindBin qw( $RealBin );
use lib $RealBin;
use DistanceMap;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Getopt::ArgvFile;
use Getopt::Long;
use Data::Dumper;
use Parsers qw(parse_pairs);

my $protein = "toy";
my $state = "nsyn";
my $verbose;
my $fdr;
my $fakenum;
my $phenotype;
my $site_pairs_file; #contains pairs of sites, one tab-delimited pair per line (e.g. 123	144\n)
my $data_folder;
my $fdr_folder;
my $newick_file;
my $site_pairs_file;
my $site2pheno_file;
my $output_folder;
my $config;



GetOptions (	
		'verbose' => \$verbose,
		'fdr' => \$fdr,
		'fakenum=i' => \$fakenum,
		'phenotype=s' => \$phenotype,
		'config=s' => \$config,
);

$| = 1;


open INFILE, $config or die "\nUnable to open config file $config!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "tree"){
			$protein = $value;
		}elsif($key eq "state"){
			$state=$value;
		}elsif($key eq "output_folder"){
			$output_folder=$value;
		}elsif($key eq "site_pairs_file"){
			$site_pairs_file = $value
		}elsif($key eq "fdr_folder"){
			$fdr_folder=$value;
		}elsif($key eq "data_folder"){
			$data_folder=$value;
		}elsif($key eq "newick_file"){
			$newick_file=$value;
		}elsif($key eq "site_pairs_file"){
			$site_pairs_file=$value;
		}elsif($key eq "site2pheno_file"){
			$site2pheno_file=$value;
		}
	}
}
close INFILE;

my $args = {protein => $protein, state => $state, 
			newick_file=> $newick_file, site_pairs_file => $site_pairs_file, site2pheno_file => $site2pheno_file,
			output_folder => $output_folder, 
			data_folder => $data_folder, fdr_folder => $fdr_folder, fdr => $fdr, 
			fakenum => $fakenum, phenotype => $phenotype};
print Dumper $args;
my $mutmap = DistanceMap->new($args);

print STDOUT "DistanceMap created\n";

$mutmap->all_nonsequential_distances();




