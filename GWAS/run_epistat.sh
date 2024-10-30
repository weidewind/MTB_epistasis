#!/usr/bin/bash
# check first argument is an existing regular file
if [[ ! -f $1 ]]
then
    echo "$1 is not a regular file" 1>&2
    exit 1
fi

drugs_fn=$1
while read drug; do
  printf "\n"
  echo "$EPISTAT_HOME/run_epistat.pl -r -P $drug/${drug}_gnu_parallel.opts -x $drug/${drug}.xparr -m 2 -p 10000 $drug/${drug}.epistat.prm run_epistat.prm"
  $EPISTAT_HOME/run_epistat.pl -r -P $drug/${drug}_gnu_parallel.opts -x $drug/${drug}.xparr -m 2 -p 10000 $drug/${drug}.epistat.prm run_epistat.prm
  echo "$EPISTAT_HOME/estimate_fdr.pl $drug/${drug}.lower.pvalue.fdr.prm > $drug/${drug}.lower.pvalue.fdr"
  $EPISTAT_HOME/estimate_fdr.pl $drug/${drug}.lower.pvalue.fdr.prm > $drug/${drug}.lower.pvalue.fdr
  echo "$EPISTAT_HOME/estimate_fdr.pl $drug/${drug}.upper.pvalue.fdr.prm > $drug/${drug}.upper.pvalue.fdr"
  $EPISTAT_HOME/estimate_fdr.pl $drug/${drug}.upper.pvalue.fdr.prm > $drug/${drug}.upper.pvalue.fdr
done <$drugs_fn