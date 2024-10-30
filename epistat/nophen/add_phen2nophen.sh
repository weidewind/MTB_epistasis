#!/usr/bin/env bash
phen_dir=../phen
epistat_dir=..

$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 4,5 -p phen. 9drugs.pairs2.pos.FDR10+noncosec.tab $phen_dir/9drugs.z-scores.all.pairs2+upval_FDR.tab>9drugs.pairs2.pos.FDR10+noncosec+phen_FDR.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 4,5 -p phen. 9drugs.sites2.neg.FDR10+noncosec.tab $phen_dir/9drugs.z-scores.all.sites2+lpval_FDR.tab>9drugs.sites2.neg.FDR10+noncosec+phen_FDR.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 4,5 -p phen. 9drugs.pairs2.pos.FDR10+noncosec.BH_greater_FDR10.tab $phen_dir/9drugs.z-scores.all.pairs2+upval_FDR.tab>9drugs.pairs2.FDR10.pos_epistatic+noncosec+phen_FDR.tab
$epistat_dir/site2gene_pos.pl -c 1,2 -a $epistat_dir/mut_codes_rc2.txt -g $epistat_dir/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff 9drugs.pairs2.FDR10.pos_epistatic+noncosec+phen_FDR.tab>9drugs.pairs2.FDR10.pos_epistatic+noncosec+w_gene_names+phen_FDR.tab
$epistat_dir/who/add_who_mutation_catalog_info.pl -c 2,6 9drugs.pairs2.FDR10.pos_epistatic+noncosec+w_gene_names+phen_FDR.tab $epistat_dir/who/WHO-UCN-GTB-PCI-2021.7-GenomeIndices.tab>9drugs.pairs2.FDR10.pos_epistatic+noncosec+WHO+phen_FDR.tab
$epistat_dir/calc_kendall4who_grade_scores.pl -g 4,9 -s 26 -o .phen_FDR2WHO_grade.cor.out 9drugs.pairs2.FDR10.pos_epistatic+noncosec+WHO+phen_FDR.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 1,2 -p phen. 9drugs.pairs2.pos.FDR10+noncosec.BH_greater_FDR10.tab $phen_dir/9drugs.pairs2.pos.FDR10+noncosec.BH_greater_FDR10.tab>9drugs.pairs2.FDR10.pos_epistatic.join_phen.tab
$EPISTAT_HOME/utils/join_site_pair_tables.pl -s 1,2 -p phen. 9drugs.sites2.neg.FDR10+noncosec.BH_greater_FDR10.tab $phen_dir/9drugs.sites2.neg.FDR10+noncosec.BH_greater_FDR10.tab>9drugs.sites2.FDR10.neg_epistatic.join_phen.tab

$epistat_dir/site2gene_pos.pl -c 1,2 -a $epistat_dir/mut_codes_rc2.txt -g $epistat_dir/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff 9drugs.pairs2.pos.FDR10+noncosec+phen_FDR.tab>9drugs.pairs2.pos.FDR10+noncosec+w_gene_names+phen_FDR.tab
$epistat_dir/who/add_who_mutation_catalog_info.pl -c 2,6 9drugs.pairs2.pos.FDR10+noncosec+w_gene_names+phen_FDR.tab $epistat_dir/who/WHO-UCN-GTB-PCI-2021.7-GenomeIndices.tab>9drugs.pairs2.pos.FDR10+noncosec+WHO+phen_FDR.tab
$epistat_dir/calc_kendall4who_grade_scores.pl -g 4,9 -s 25 -o .phen_FDR2WHO_grade.cor.out 9drugs.pairs2.pos.FDR10+noncosec+WHO+phen_FDR.tab