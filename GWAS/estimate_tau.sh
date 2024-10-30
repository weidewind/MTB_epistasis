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
  echo "$EPISTAT_HOME/estimate_tau.pl -b 0.95 $drug/${drug}.xparr>$drug/${drug}.xparr.tau"
  $EPISTAT_HOME/estimate_tau.pl -b 0.95 $drug/${drug}.xparr>$drug/${drug}.xparr.tau
done <$drugs_fn