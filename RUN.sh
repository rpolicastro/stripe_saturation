#!/bin/bash

CONTAINER='library://rpolicastro/default/stripeseq_saturation:1.0.0'

##########################################
## Saturation Analysis of STRIPE-seq Data
##########################################

## Load Settings
## ----------

source settings.conf

## Prepare Singularity Container
## ----------

## Make sure the container directory exists.

mkdir -p ${OUTDIR}/container

## Download container if it doesn't exist.

CONTAINER_NAME=$(basename $CONTAINER | tr ":" "_").sif

if [ ! -f ${OUTDIR}/container/$CONTAINER_NAME ]; then
	singularity pull -U $CONTAINER
	mv $CONTAINER_NAME ${OUTDIR}/container
fi

## Run Workflow
## ----------

SINGULARITYENV_OUTDIR=$OUTDIR
SINGUALRITYENV_CORES=$CORES
SINGUALRITYENV_PAIRED=$PAIRED
SINGUALRITYENV_GTF=$GTF
SINGUALRITYENV_FROM=$FROM
SINGUALRITYENV_TO=$TO
SINGUALRITYENV_BY=$BY

singularity exec \
-eCB $PWD,$OUTDIR,$(dirname ${GTF}) \
-H $PWD \
main.sh
