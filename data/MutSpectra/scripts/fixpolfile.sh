#!/usr/bin/env bash

source .env
set -euo pipefail

rawFN=${PATH2ANCESTOR}/homo_sapiens_ancestor_

fastaFiles=""
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22;
do
    echo "${rawFN}${i}.fa"
    fastaFiles="${fastaFiles} ${rawFN}${i}.fa"
done

path2out="${PATH2ANCESTOR}/hs_ancestor.fa"

# gsed is gnu-sed
cat ${fastaFiles} | \
    sed 's/>ANCESTOR_for_chromosome:GRCh37:/>chr/' | \
    sed 's/:.*$//' >> ${path2out}
