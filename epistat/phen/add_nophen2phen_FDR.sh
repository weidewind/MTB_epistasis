#!/usr/bin/env bash
nophen_dir=../nophen

$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 4,5 -p nophen. 9drugs.pairs2.pos.FDR10+noncosec.tab $nophen_dir/9drugs.z-scores.all.pairs2+upval_FDR.tab>9drugs.pairs2.pos.FDR10+noncosec+nophen_FDR.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 4,5 -p nophen. 9drugs.sites2.neg.FDR10+noncosec.tab $nophen_dir/9drugs.z-scores.all.sites2+lpval_FDR.tab>9drugs.sites2.neg.FDR10+noncosec+nophen_FDR.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 1,2 -p nophen. 9drugs.pairs2.pos.FDR10+noncosec.BH_greater_FDR10.tab $nophen_dir/9drugs.pairs2.pos.FDR10+noncosec.BH_greater_FDR10.tab>9drugs.pairs2.FDR10.pos_epistatic.join_nophen.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 1,2 -p nophen. 9drugs.sites2.neg.FDR10+noncosec.BH_greater_FDR10.tab $nophen_dir/9drugs.sites2.neg.FDR10+noncosec.BH_greater_FDR10.tab>9drugs.sites2.FDR10.neg_epistatic.join_nophen.tab