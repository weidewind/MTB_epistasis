#!/bin/bash
$EPISTAT_HOME/mk_coevolmtx.pl -m Z-score drugs2vars.ddss.all.mk_coevolmtx.3.prm
$EPISTAT_HOME/minvert/cor2pcor.R -f mtb_RR.all.block.mtx -t 1 -l 0.9 -s .abs.l90.APC.cor2pcor.R.out