#!/usr/bin/env bash
EPISTAT7_HOME=path/to/installed/EpiStat_v7.1
fdr_dir="9drugs.filtered/fdr"
find $fdr_dir/*.xpar|parallel $EPISTAT7_HOME/epistat.pl -p -x {} cumdist.epistat.prm

$EPISTAT7_HOME/epistat.pl -p -x 9drugs.xparr cumdist.epistat.prm