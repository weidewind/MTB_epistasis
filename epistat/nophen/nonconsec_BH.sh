#!/usr/bin/env bash
$EPISTAT_HOME/stat/mk_BH_fdr_correction.pl -s 15 -a 0.1 9drugs.pairs2.pos.FDR10+noncosec.tab
$EPISTAT_HOME/stat/mk_BH_fdr_correction.pl -s 13 -a 0.1 9drugs.sites2.neg.FDR10+noncosec.tab