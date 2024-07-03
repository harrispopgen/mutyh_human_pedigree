#!/usr/bin/env bash

module load modules modules-init
module load python/3.7.7
module load IGV/2.9.4
module load igv-reports/1.7.0

eqFile="TMP.nobackup/src/reference/equivalence.csv"

sample1=$(egrep ${snakemake_wildcards[sample]} ${eqFile} | awk 'BEGIN {FS=","} {print($2)}')
parent1=$(egrep ${snakemake_wildcards[trio1]} ${eqFile} | awk 'BEGIN {FS=","} {print($2)}')
parent2=$(egrep ${snakemake_wildcards[trio2]} ${eqFile} | awk 'BEGIN {FS=","} {print($2)}')
parent3=$(egrep ${snakemake_wildcards[trio3]} ${eqFile} | awk 'BEGIN {FS=","} {print($2)}')

create_report ${snakemake_input[filtvcf]} \
		      --genome hg38 \
                      --flanking 10000 \
                      --info-columns AC QD hiConfDeNovo \
                      --samples $sample1 $parent1 $parent2 $parent3 \
                      --sample-columns AD GQ GT PL \
                      --tracks ${snakemake_input[bam]} \
                               ${snakemake_input[bam_nodups]} \
                               ${snakemake_input[bam_realign]} \
                               ${snakemake_input[bamtrio1]} \
                               ${snakemake_input[bamtrio1_nodups]} \
                               ${snakemake_input[bamtrio1_realign]} \
                               ${snakemake_input[bamtrio2]} \
                               ${snakemake_input[bamtrio2_nodups]} \
                               ${snakemake_input[bamtrio2_realign]} \
                               ${snakemake_input[bamtrio3]} \
                               ${snakemake_input[bamtrio3_nodups]} \
                               ${snakemake_input[bamtrio3_realign]} \
		      --output ${snakemake_output[html]}
