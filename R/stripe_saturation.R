#!/usr/bin/env Rscript

library("tidyverse")
library("getopt")
library("ChIPseeker")
library("GenomicFeatures")
library("GenomicRanges")
library("TSRchitect")

##########################################
## Saturation Analysis of STRIPE-seq Data
##########################################

## Load Settings
## ----------

options <- matrix(c(
	"outdir", "o", 1, "character", "output directory",
	"cores", "c", 1, "character", "number of cores",
	"gtf", "g", 1, "character", "gtf annotation file",
	"paired", "p", 1, "character", "whether the bams are paired or not"
), byrow=TRUE, ncol=5)

opt <- getopt(options)

## Make Sample Sheet
## ----------

## Get sample files.

sample.files <- list.files(file.path(opt$outdir, "sampled_bams"), pattern=".*\\.bam$")

## Get sample names.

sample.names <- str_replace(sample.files, "_Aligned.out.bam", "")

## Make sample sheet.

samples <- tibble(
	"SAMPLE"=sample.names,
	"ReplicateID"=1,
	"FILE"=file.path(opt$outdir, "sampled_bams", sample.files)
)

## Export sample sheet.

write.table(
	samples, file.path(opt$outdir, "sampled_bams", "samples.tsv"),
	sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

## Find TSSs
## ----------

## Load TSS object.

tss.obj <- loadTSSobj(
	experimentTitle="STRIPE-seq Saturation Analysis",
	inputDir=file.path(opt$outdir, "sampled_bams"),
	n.cores=as.numeric(opt$cores),
	isPairedBAM=as.logical(opt$paired),
	isPairedBED=FALSE,
	sampleSheet="samples.tsv",
	sampleNames="",
	replicateIDs=""
)

## Get fragment 5` ends.

tss.obj <- inputToTSS(tss.obj)

## Find TSSs.

tss.obj <- processTSS(
	tss.obj,
	n.cores=as.numeric(opt$cores),
	tssSet="all",
	writeTable=FALSE
)

## Grab TSSs from object.

TSSs <- tss.obj@tssCountData %>%
	setNames(tss.obj@sampleNames)

## Convert TSSs to GRanges object.

TSSs.granges <- map(
	TSSs,
	~ GRanges(
		seqnames=.$seq,
		ranges=IRanges(start=.$TSS, width=1),
		strand=.$strand,
		score=.$nTAGs
	)
)

## Export TSS Granges.

dir.create(file.path(opt$outdir, "Rdata"))

saveRDS(TSSs.granges, file.path(opt$outdir, "Rdata", "TSSs.RDS"))

## Annotate TSSs
## ----------

## Load genome annotation.

annotation <- makeTxDbFromGFF(opt$gtf, format="gtf")

## Annotate TSSs.

TSSs.annotated <- map(
	TSSs.granges,
	~ annotatePeak(.,
		tssRegion=c(-1000,100),
		TxDb=annotation,
		sameStrand=TRUE
	) %>%
	as_tibble
)

## Saturation Analysis
## ----------

## Saturation analysis function.

saturation.analysis <- function(x) {
	n.genes <- x %>%
		filter(score >= 3) %>%
		pull(geneId) %>%
		unique %>%
		length

	promoter.proximal.frac <- x %>%
		filter(score >= 3) %>%
		count(annotation) %>%
		mutate(frac=n/sum(n)) %>%
		filter(annotation == "Promoter") %>%
		pull(frac)

	return(tibble(n.genes, promoter.proximal.frac))
}

## Get saturation.

TSS.saturation <- map(TSSs.annotated, ~saturation.analysis(.))

## Plot Saturation
## ----------

## Create outdir directory.

dir.create(file.path(opt$outdir, "saturation_analysis_plots"))

## Format data for plotting.

plotting.data <- TSS.saturation %>%
	bind_rows(.id="sample") %>%
	mutate(sample_value=str_remove(sample, pattern="_sampled.+")) %>%
	mutate(sample_value=as.numeric(sample_value)) %>%
	arrange(sample_value)

## Plot data.

p <- ggplot(plotting.data, aes(x=sample_value, y=n.genes)) +
	geom_line() +
	geom_point(aes(size=promoter.proximal.frac, color=promoter.proximal.frac)) +
	scale_color_viridis_c() +
	theme_bw() +
	labs(
		x="Number of Sampled Reads",
		y="Genes with Promoter Proximal TSS"
	)

ggsave(
	file.path(opt$outdir, "saturation_analysis_plots", "saturation_plot.pdf"),
	plot=p, device="pdf", height=4, width=7
)
