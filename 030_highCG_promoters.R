# Libraries ---------------------------------------------------------------

library(SummarizedExperiment)
library(org.Bt.eg.db)

# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'
genomeDir <- "bostaurus"
plotDir <- "collapsedGroups"
highCGdir <- file.path(plotDir, "highCG")

# Create folders ----

if (!dir.exists(highCGdir)){
  dir.create(highCGdir)
}

# Import promoter information ----

promoters.5 <- readRDS(file.path(outdir, "promoters.5.rds"))

highPromoters <- promoters.5[promoters.5$CG.promoter > 100,]

range(highPromoters$CG.promoter)  # all have 3 digits

# Import smoothed methylation collapsed by group ----

BS.infection.smooth <- readRDS(file.path(outdir, "BS.infection.smooth.rds"))

# Save line colour in phenoData slot
BS.infection.smooth$Infection <- sampleNames(BS.infection.smooth)
BS.infection.smooth$col <- ifelse(
  BS.infection.smooth$Infection == "Control", "blue", "red"
)
pData(BS.infection.smooth)

# Plot regions ----

genes.gr <- readRDS(file.path(genomeDir, "biomart_ensembl_genes.rds"))

exons.df <- readRDS(file.path(genomeDir, "biomart_ensembl_exons.rds"))
exons.gr <- GRanges(
  seqnames = exons.df$chromosome_name,
  ranges = IRanges(
    start = exons.df$exon_chrom_start,
    end = exons.df$exon_chrom_end,
    names = exons.df$ensembl_gene_id
  ),
  strand = exons.df$strand,
  external_gene_name = exons.df$external_gene_name
)
exons.grl <- split(exons.gr, names(exons.gr))


ensGenes2plot <- names(highPromoters)

for (ensGene in ensGenes2plot){
  geneGR <- genes.gr[ensGene]
  geneName <- geneGR$external_gene_name
  promoterGR <- promoters(geneGR, 1500, 500)
  annotGR <- list(
    Promoter = promoterGR,
    Gene = geneGR,
    Exons = exons.grl[[ensGene]]
  )
  CGcount <- mcols(highPromoters[ensGene])[,"CG.promoter"]
  MethPromoter <- format(
    mcols(highPromoters[ensGene])[,"Meth.promoter"], digits = 2)
  names(annotGR)[2] <- ifelse(geneName == "", ensGene, geneName)
  pdf(
    file = file.path(highCGdir, sprintf(
      "%i-%s-%s-%s_gene_exons_promoter.pdf",
      CGcount, MethPromoter, geneName, ensGene)),
    width = 9, height = 7
  )
  plotRegion(
    BS.infection.smooth, promoterGR, extend = 2E3,
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
}
