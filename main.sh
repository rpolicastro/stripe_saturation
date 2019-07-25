#!/bin/bash

## Downsampling BAMs
## ----------

## Get total read number.

READS=$(samtools flagstat $BAM | awk 'NR==1' | cut -d" " -f1)

## Divide number by 2 if paired.

if [ $PAIRED = TRUE ]; then
        READS=$(Rscript --slave -e 'cat(${READS}/2)')
fi

## Get sampling percentages.

SAMPLES=($(seq $FROM $BY $TO))

SAMPLE_FRACS=()
for SAMPLE in ${SAMPLES[@]}; do
        SAMPLE_FRACS+=($(R --slave -e 'cat(${SAMPLE}/${READS})'))
done

## Sample BAMs and remove PCR duplicates.

mkdir -p ${OUTDIR}/sampled_bams

N_SAMPLES=${#SAMPLE_FRACS[@]}
for N in $(seq 0 1 $(bc <<< $N_SAMPLES-1)); do
        samtools view -bs 0${SAMPLE_FRACS[$N]} $BAM | \
        samtools sort -n -@ $CORES - | \
        samtools fixmate -m - - | \
        samtools sort -@ $CORES - | \
        samtools markdup - - | \
        samtools view -F 3852 -f 3 -O BAM -@ $CORES -o ${OUTDIR}/sampled_bams/${SAMPLES[$N]}_sampled_$(basename $BAM)

        samtools index ${OUTDIR}/sampled_bams/${SAMPLES[$N]}_sampled_$(basename $BAM)
done

## Saturation Analysis
## ----------

Rscript R/stripe-saturation.R \
--outdir $OUTDIR \
--paired $PAIRED \
--cores $CORES \
--gtf $GTF
