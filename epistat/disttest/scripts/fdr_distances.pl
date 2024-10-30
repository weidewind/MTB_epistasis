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
use Parsers;
use Parallel::ForkManager;
use File::Path qw(make_path remove_tree);

my $protein = "toy";
my $state = "nsyn";
my $phenotype;
my $verbose;
# my $switch;
my $site_pairs_file;
my $maxfakenum;
my $fdr_folder;
my $data_folder;
my $newick_file;
my $site_pairs_file;
my $site2pheno_file;
my $output_folder;
my $config;

GetOptions (	
		'verbose' => \$verbose,
		'phenotype=s' => \$phenotype,
		'maxfakenum=i' => \$maxfakenum,
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
			data_folder => $data_folder, fdr_folder => $fdr_folder, fdr => 1, phenotype => $phenotype};
my $folder = DistanceMap::get_input_base($args);
make_path(DistanceMap::get_output_base($args));
my $log = File::Spec->catfile(DistanceMap::get_output_base($args), "$protein.log");

open LOG, ">$log" or die "Cannot create log file $log $! \n";
print "Printing some data to logfile $log\n";

my @files = dirfiles($folder);
my @statfiles = grep {/^[0-9]+\.stat\.ident/} @files;
print LOG Dumper (\@statfiles );
my @fakenums;
for my $filename(@statfiles){
		push @fakenums, $filename =~ /([0-9]+)/;
		if (scalar @fakenums == $maxfakenum){last;}
}
print LOG Dumper (\@fakenums);
##
## Preparing commands for Forkmanager
my @commands;



foreach my $fakenum(@fakenums){	
	#check if this fake is already processed
	my $outfile = File::Spec->catfile(DistanceMap::get_output_base($args), $protein."_".$fakenum."_all_nonsequential");
	if (-e $outfile){
		print LOG "$outfile already exists. Skipping\n";
		next;
	}
	my $command = mycomm($fakenum, $site_pairs_file);
	push @commands, $command;
	print LOG $command."\n";	
}

##
## Forkmanager setup
my $manager = new Parallel::ForkManager(14);

$manager->run_on_start( 
	  sub {
		my $pid = shift;
		print LOG "Starting child processes under process id $pid\n";
	  }
	);
$manager->run_on_finish( 
	  sub {
		 my ( $pid, $exit_code, $signal, $core ) = @_;
		 if ( $core ) {
			print LOG "Process (pid: $pid) core dumped.\n";
		 } else {
			print LOG "Process (pid: $pid) exited with code $exit_code and signal $signal.\n";
		 }
	  }
   );
$manager->run_on_wait( 
	  sub {
		 print LOG "Waiting for all children to terminate... \n";
	  },
	  180 # time interval between checks
   ); 
##
## Launching a series of new iteration_gulps if needed 

foreach my $command (@commands) {
	  $manager->start and next;
	  system( $command );
	  $manager->finish;
   }
$manager->wait_all_children;

close LOG;

##
sub mycomm {
	my $fakenum = shift;
	my $site_pairs_file = shift;
	my $perlocation = "perl";
	my $exports = "";
	# if ($switch) {
		# $perlocation = "~/perl5/perlbrew/perls/perl-5.22.1/bin/perl";
	 	# $exports = "export PERL5LIB=/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/.perl/lib/perl5/5.22.1/x86_64-linux:/export/home/popova/.perl/lib/perl5/5.22.1:/export/home/popova/.perl/lib/perl5/x86_64-linux:/export/home/popova/.perl/lib/perl5:/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/workspace/evopoisson:$PERL5LIB; ";
	# }
	my $command = $exports.$perlocation." test_distances.pl --input $input --protein $protein --state $state --output_folder $output_folder --fdr --fakenum $fakenum";
	if ($site_pairs_file){
		$command .= " --site_pairs_file $site_pairs_file";
	}
	if ($phenotype){
		$command .= " --phenotype $phenotype";
	}
	return $command;
	
}

sub dirfiles{
        my $dirname = shift;
        opendir(DH, $dirname);
        my @files = readdir(DH);
        closedir(DH);
        if (scalar @files == 0){
                print LOG "No files found in $dirname!\n";
				print "No files found in $dirname!\n";
        }
        return @files
}
