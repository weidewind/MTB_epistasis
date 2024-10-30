#!/usr/bin/env bash
EPISTAT7_HOME=/path/to/the/epistat.7.1
pheno="./9drugs/pheno"
drugIdList=(1 2 4 6 8 9 12 13 14)
for d in ${drugIdList[@]}; do
   seq 500|parallel $EPISTAT7_HOME/epistat.pl -p -x $pheno/$d/samples/L0I0/{}.xpar cumdist.epistat.prm
done

$EPISTAT7_HOME/epistat.pl -p -x 9drugs.xparr cumdist.epistat.prm