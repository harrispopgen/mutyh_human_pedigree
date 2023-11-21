#!/bin/bash
#$ -N coverage_job
#$ -o /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/errors
#$ -e /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/errors
#$ -m bea
#$ -t 1-345
#$ -cwd

### run samtools depth command to compute coverage scores for individuals across 10K bp chunks; add zeros to any empty rows in output ###

module load modules modules-init
module load samtools/1.9

BAM_PATH=/net/harris/vol1/sharing/MUTYH_pedigree/bamFiles/
OUTPUT_DIR=/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/output
INDIVIDUALS=(C11.1079256 C12.1079257 C21.1079259 C22.1079260 C23.1079261 C31.1079264 C32.1079265 C41.1079267 C42.1079268 P1.1079254 P2.1079258 P3.1079262 P4.1079266 S1.1079255 S3.1079263)
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
CHUNK_SIZE=10000

INDEX=$(expr $SGE_TASK_ID - 1)
INDIVIDUAL_INDEX=$((INDEX / ${#CHROMOSOMES[@]}))
CHROMOSOME_INDEX=$((INDEX % ${#CHROMOSOMES[@]}))

INDIVIDUAL=${INDIVIDUALS[${INDIVIDUAL_INDEX}]}
CHROMOSOME=${CHROMOSOMES[${CHROMOSOME_INDEX}]}

BAM_FILE="${BAM_PATH}/${INDIVIDUAL}.bam"
OUT_FILE="${OUTPUT_DIR}/${INDIVIDUAL}_${CHROMOSOME}.txt"
JOB_NAME="${INDIVIDUAL}_${CHROMOSOME}_coverage_job"

# get chromosome length
CHR_LENGTH=$(samtools view -H ${BAM_FILE} | grep "^@SQ" | grep -P "${CHROMOSOME}\t" | cut -f 3 -d ':' | cut -f 2)

# divide chromosome into chunks and get the coverage for each chunk, and run samtools depth with -Q 20 for minimum mapq
for ((START_POS=1; START_POS<=${CHR_LENGTH}; START_POS+=${CHUNK_SIZE})); do
  END_POS=$((START_POS+CHUNK_SIZE-1))
  if ((END_POS > CHR_LENGTH)); then
    END_POS=${CHR_LENGTH}
  fi
  DEPTH=$(samtools depth -r ${CHROMOSOME}:${START_POS}-${END_POS} -Q 20 ${BAM_FILE} | awk '{sum+=$3} END {print sum/NR}')
  echo -e "${CHROMOSOME}\t${START_POS}\t${END_POS}\t${DEPTH}" >> ${OUT_FILE}
done
