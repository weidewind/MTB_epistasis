#!/usr/bin/env bash

config=${1}

perl test_distances.pl --config $config >${outfolder}/disttest_mtub.log 2>disttest.errlog
perl fdr_distances.pl --config $config --maxfakenum 500 >fdrdist.out 2>fdrdist.err # use no more than 500 fakes 
perl grep_distances_fdr.pl --config $config --sites --nopdb