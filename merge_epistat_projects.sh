#!/bin/bash
$EPISTAT_HOME/union_epistat_projects.pl -p gwas.prj_dir.list -o merged_GWAS GWAS/Amikacin/Amikacin.epistat.prm GWAS/run_epistat.prm
$EPISTAT_HOME/union_epistat_projects.pl -p ddss.prj_dir.list -o merged_DDSS DDSS/Amikacin/Amikacin.epistat.prm DDSS/run_epistat.prm