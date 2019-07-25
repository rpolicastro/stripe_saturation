#!/bin/bash

##########################################
## Saturation Analysis of STRIPE-seq Data
##########################################

## Load Settings
## ----------

source settings.conf

## Downsampling BAMs
## ----------

## Get total read number.

READS=$(samtools flagstat $BAM | awk 'NR==1' | cut -d" " -f1)

## Divide number by 2 if paired.

[[ $PAIRED = TRUE ]] && READS=$(python -c "print $READS/2")
