#!/bin/bash
set -euo pipefail
awk '{if (substr($0, 0, 1) == "#") print ""; else printf("%s ", $4)}END{print "";}' stem.txt | grep -v '^\s*$' > stem_0.txt
NUM=expr `grep "\-\-Z" stem.txt | wc -l`

if [ $NUM -gt 0 ]
then
seq 1 $NUM | xargs -I{} echo "tail -n +{} stem_0.txt | head -n 1 > stem_{}.txt" | bash
awk '{if (substr($0, 0, 1) == "#") print ""; else printf("%s ", $4)}END{print "";}' acc.txt | grep -v '^\s*$' > acc_0.txt
NUM=`grep "\-\-Z" acc.txt | wc -l`
seq 1 $NUM | xargs -I{} echo "tail -n +{} acc_0.txt | head -n 1 > acc_{}.txt" | bash
R --vanilla --slave --args $NUM < plot.R
fi