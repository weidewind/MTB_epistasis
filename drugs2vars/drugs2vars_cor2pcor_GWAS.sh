#!/bin/bash
$EPISTAT_HOME/mk_coevolmtx.pl -m Z-score drugs2vars.gwas.mk_coevolmtx.3.prm
$EPISTAT_HOME/minvert/cor2pcor.R -f mtb.block.mtx -l 0.9 -s .l90.APC.cor2pcor.R.out