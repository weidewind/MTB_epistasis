#!/usr/bin/env bash

config=${1}

perl test_distances.pl --config $config >disttest_pheno_mtub.log 2>disttest_pheno.errlog
parallel perl fdr_distances.pl --config $config --phenotype {1} --maxfakenum 500 >fdrdist_pheno_mtub_{1}.out 2>fdrdist_pheno_{1}.err ::: 1 2 4 6 8 9 12 13 14 # use no more than 500 fakes
perl grep_pheno_distances_fdr_stepwise.pl --config $config --sites --nopdb --numfakes 500 --numchunks 10 2>fdrdist_pheno_grepper.err
 
