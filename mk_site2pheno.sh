#!/bin/bash
cp merged_GWAS/gwas.site_pairs merged_min_GWAS_DDSS
$EPISTAT_HOME/utils/two_vec_min.pl merged_GWAS/gwas.lower.pvalue merged_DDSS/ddss.lower.pvalue>merged_min_GWAS_DDSS/min.lower.pvalue
$EPISTAT_HOME/utils/two_vec_min.pl merged_GWAS/gwas.upper.pvalue merged_DDSS/ddss.upper.pvalue>merged_min_GWAS_DDSS/min.upper.pvalue
cd merged_min_GWAS_DDSS
$EPISTAT_HOME/conv_epistat_prj2site2pheno.pl -p gwas.site_pairs -P 1 -l min.lower.pvalue -u min.upper.pvalue -r 6,13,12,4,14,1,2,9,8>min_RR_SR.site2pheno
