#!/usr/bin/env bash
##############################
disttest_dir=../disttest/nsyn/nopheno
epistat_dir=..
#obtaining the list of discordantly evolving pairs for 10% FDR pvalue threshold
# using the prepared configuration file 9drugs.sites2.pos.FDR10.mk_summary.cfg (need to be edited by user!!!):
$EPISTAT_HOME/mk_summary.pl -u -a 9drugs.sites2.neg.FDR10.mk_summary.cfg -s 1 9drugs.all.mtx>9drugs.sites2.neg.FDR10.tab
#obtaining the list of concordantly evolving pairs for 10% FDR pvalue threshold
# using the prepared configuration file 9drugs.pairs2.pos.FDR10.mk_summary.cfg (need to be edited by user!!!):
$EPISTAT_HOME/mk_summary.pl -u -a 9drugs.pairs2.pos.FDR10.mk_summary.cfg 9drugs.n0.l90.cor2pcor.R.out>9drugs.pairs2.pos.FDR10.tab
$epistat_dir/site2gene_pos.pl -c 1,2 -a $epistat_dir/mut_codes_rc2.txt -g $epistat_dir/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff 9drugs.pairs2.pos.FDR10.tab>9drugs.pairs2.pos.FDR10+w_gene_names.tab
$epistat_dir/site2gene_pos.pl -c 1,2 -a $epistat_dir/mut_codes_rc2.txt -g $epistat_dir/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff 9drugs.sites2.neg.FDR10.tab>9drugs.sites2.neg.FDR10+w_gene_names.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 8,9,10 -p nconsec. 9drugs.pairs2.pos.FDR10.tab $disttest_dir/9drugs.filtered_dist.fdr_results>9drugs.pairs2.pos.FDR10+noncosec.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 8,9,10 -p nconsec. 9drugs.sites2.neg.FDR10.tab $disttest_dir/9drugs.filtered_dist.fdr_results>9drugs.sites2.neg.FDR10+noncosec.tab
