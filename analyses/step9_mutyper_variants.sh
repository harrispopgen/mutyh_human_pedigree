#! /bin/bash
#$ -l h_rt=100:00:00,h_data=10G
#$ -o /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/mutyper_3mers/errors
#$ -e /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/mutyper_3mers/errors
#$ -m bea

module load modules modules-init modules-gs
module load python/3.7.7
module load samtools/1.9
module load htslib/1.9 bcftools/1.9

set -o pipefail

vcfdirectory=/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/mutyper_3mers/vcfs
refgenome=/net/harris/vol1/data/hg38/GRCh38_full_analysis_set_plus_decoy_hla_singlespaced.fa

sep="\s"
kmersize=3

outdir=/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/mutyper_3mers/output
variantsdir=$outdir/mutyper_variant_files
spectrumdir=$outdir/mutyper_spectrum_files

mkdir -p $outdir
mkdir -p $variantsdir
mkdir -p $spectrumdir

# Select SAMPLE based on task ID
SAMPLE_LIST=($(ls $vcfdirectory/*.vcf.gz | xargs -n 1 basename | sed 's/\.vcf\.gz//'))

if [ ${#SAMPLE_LIST[@]} -eq 0 ]; then
    echo "No samples found."
    exit 1
fi

# Loop through all samples
for SAMPLE in "${SAMPLE_LIST[@]}"; do
    # Extract the "pure" sample name (without appended chromosome name)
    PURE_SAMPLE=$(echo $SAMPLE | awk -F '_' '{print $1}') 

    vcffilename=$vcfdirectory/${SAMPLE}.vcf.gz

    # Check if index file already exists before indexing
    if [ ! -f "${vcffilename}.csi" ]; then
        bcftools index $vcffilename || { echo "Indexing failed for $vcffilename"; continue; }
    else
        echo "Index file already exists for $vcffilename. Skipping indexing."
    fi

    # Get list of unique chromosomes present in the sample's VCF file
    CHROM_LIST=$(bcftools view -H $vcffilename | cut -f 1 | sort -u)

    # Loop through all chromosomes
    for CHROM in $CHROM_LIST; do
        chromvcffile=$vcfdirectory/${PURE_SAMPLE}_$CHROM.vcf.gz

        # Check if the chromosome-specific VCF file already exists before creating a new one
        if [ ! -f "$chromvcffile" ]; then
            bcftools view -r $CHROM $vcffilename -Oz -o $chromvcffile || { echo "bcftools view failed"; continue; }
        else
            echo "Chromosome-specific VCF file already exists for $SAMPLE on $CHROM. Skipping extraction."
        fi

        # Check if index file already exists before indexing
        if [ ! -f "${chromvcffile}.csi" ]; then
            bcftools index $chromvcffile || { echo "Indexing failed for $chromvcffile"; continue; }
        else
            echo "Index file already exists for $chromvcffile. Skipping indexing."
        fi

        hetvariantsoutfile=$variantsdir/${PURE_SAMPLE}_${CHROM}.mutyper.variants.mutationTypes.IGVfiltered.NOSTRICT.vcf.gz
        zcat $chromvcffile | grep -v "\.\/\." | bcftools view -c 1:minor | mutyper variants --k $kmersize --sep $sep $refgenome - | bgzip -cf > ${hetvariantsoutfile}
        exitVal=$?

        if [ ! -s "$hetvariantsoutfile" ]; then
            echo "Intermediate file $hetvariantsoutfile is empty."
            continue
        fi

        if [ ${exitVal} -ne 0 ]; then
                echo "error in mutyper variants for $SAMPLE on $CHROM"
                continue
        else
                echo "finished variant calling for $SAMPLE on $CHROM"
        fi

        hetspectrumoutfile=$spectrumdir/${PURE_SAMPLE}_${CHROM}.mutyper.spectrum.IGVfiltered.NORANDOMIZE.NOSTRICT.txt
        mutyper spectra $hetvariantsoutfile > $hetspectrumoutfile
        exitVal=$?

        if [ ${exitVal} -ne 0 ]; then
                echo "error in mutyper spectra for $SAMPLE on $CHROM"
                continue
        else
                echo "finished spectra calling for $SAMPLE on $CHROM"
        fi

    done
done






