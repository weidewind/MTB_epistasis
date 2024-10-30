#!/usr/bin/env bash
epistat_dir=..

$epistat_dir/site2gene_pos.pl -c 1,2 -a $epistat_dir/mut_codes_rc2.txt -g $epistat_dir/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff 9drugs.pairs2.pos.FDR10+noncosec.tab>9drugs.pairs2.pos.FDR10+noncosec+w_gene_names.tab
$epistat_dir/site2gene_pos.pl -c 1,2 -a $epistat_dir/mut_codes_rc2.txt -g $epistat_dir/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff 9drugs.sites2.neg.FDR10+noncosec.tab>9drugs.sites2.neg.FDR10+noncosec+w_gene_names.tab
$epistat_dir/site2gene_pos.pl -c 1,2 -a $epistat_dir/mut_codes_rc2.txt -g $epistat_dir/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff 9drugs.pairs2.pos.FDR10+noncosec.BH_greater_FDR10.tab>9drugs.pairs2.FDR10.pos_epistatic+noncosec+w_gene_names.tab
$epistat_dir/site2gene_pos.pl -c 1,2 -a $epistat_dir/mut_codes_rc2.txt -g $epistat_dir/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff 9drugs.sites2.neg.FDR10+noncosec.BH_greater_FDR10.tab>9drugs.sites2.FDR10.neg_epistatic+noncosec+w_gene_names.tab
$epistat_dir/who/add_who_mutation_catalog_info.pl -c 2,6 9drugs.pairs2.FDR10.pos_epistatic+noncosec+w_gene_names.tab $epistat_dir/who/WHO-UCN-GTB-PCI-2021.7-GenomeIndices.tab>9drugs.pairs2.FDR10.pos_epistatic+noncosec+WHO.tab
$epistat_dir/who/add_who_mutation_catalog_info.pl -c 2,6 9drugs.sites2.FDR10.neg_epistatic+noncosec+w_gene_names.tab $epistat_dir/who/WHO-UCN-GTB-PCI-2021.7-GenomeIndices.tab>9drugs.sites2.FDR10.neg_epistatic+noncosec+WHO.tab
