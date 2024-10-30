#!/usr/bin/env bash
$EPISTAT_HOME/run_epistat.pl -r -x 9drugs.filtered.xparr -m 2 -p 10000 epistat.phen.prm run_epistat.prm
#estimating FDR for pvalues:
cat estimate_fdr.txt|parallel $EPISTAT_HOME/estimate_fdr.pl "{}.prm>{}"
$EPISTAT_HOME/mk_coevolmtx.pl -m Z-score 9drugs.mk_coevolmtx.prm
$EPISTAT_HOME/minvert/cor2pcor.R -f 9drugs.block.mtx -n 0 -l 0.9 -s .n0.l90.cor2pcor.R.out
$EPISTAT_HOME/mk_coevolmtx.pl -m Z-score 9drugs.all.mk_coevolmtx.prm