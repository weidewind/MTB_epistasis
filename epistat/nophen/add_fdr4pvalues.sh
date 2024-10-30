#!/usr/bin/env bash
$EPISTAT_HOME/mk_summary.pl -u -a 9drugs.sites2.all.mk_summary.cfg -s 1 9drugs.all.mtx>9drugs.z-scores.all.sites2.tab
$EPISTAT_HOME/mk_summary.pl -u -a 9drugs.pairs2.all.mk_summary.cfg -s 1 9drugs.all.mtx>9drugs.z-scores.all.pairs2.tab
$EPISTAT_HOME/utils/add_fdr4pvalues.pl -f 9drugs.upper.pvalue.unord_pairs.pairs2.unlinked.fdr -p 4 9drugs.z-scores.all.pairs2.tab>9drugs.z-scores.all.pairs2+upval_FDR.tab
$EPISTAT_HOME/utils/add_fdr4pvalues.pl -f 9drugs.upper.pvalue.unord_pairs.sites2.unlinked.fdr -p 5 9drugs.z-scores.all.sites2.tab>9drugs.z-scores.all.sites2+upval_FDR.tab
$EPISTAT_HOME/utils/add_fdr4pvalues.pl -f 9drugs.lower.pvalue.unord_pairs.sites2.unlinked.fdr -p 4 9drugs.z-scores.all.sites2.tab>9drugs.z-scores.all.sites2+lpval_FDR.tab