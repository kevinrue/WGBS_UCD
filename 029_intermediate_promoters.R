
# Libraries ---------------------------------------------------------------

library(SummarizedExperiment)
library(org.Bt.eg.db)

# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'
genomeDir <- "bostaurus"
plotDir <- "collapsedGroups"
immuneDir <- file.path(plotDir, "immune_defense")
imprintedDir <- file.path(plotDir, "imprinted")

# Create folders ----

if (!dir.exists(immuneDir)){
  dir.create(immuneDir)
}
if (!dir.exists(imprintedDir)){
  dir.create(imprintedDir)
}

# Import promoter information ----

promoters.5 <- readRDS(file.path(outdir, "promoters.5.rds"))

# Identify high-CG promoters ----

head(promoters.5)
length(promoters.5)
names(mcols(promoters.5))

# reorder by highest CG content
promoters.5 <- promoters.5[order(promoters.5$CG.promoter, decreasing = TRUE)]

# Subset to those in the range 33-66% methylation
intermediate.5 <- promoters.5[
  promoters.5$Meth.promoter > 1/3 & promoters.5$Meth.promoter < 2/3
  ]

length(intermediate.5)
length(intermediate.5)/length(promoters.5)

head(intermediate.5)

tmp <- mcols(intermediate.5)

write.csv(tmp, "intermediate_promoters/topCGcontent.csv", row.names = FALSE)

# Preprocess methylation calls ----

# Load sites ready for smoothiin
BS <- readRDS(file.path(outdir, "BS.rmZero.rds"))

# Collapse replicates
BS.infection <- collapseBSseq(
  BSseq = BS,
  columns = ifelse(grepl("^C", colnames(BS)), "Control", "M. bovis")
)
rm(BS)

BS.infection.smooth <- BSmooth(BS.infection, mc.cores = 2, verbose = TRUE)
BS.infection.smooth
saveRDS(BS.infection.smooth, file.path(outdir, "BS.infection.smooth.rds"))

rm(BS.infection)

# Save line colour in phenoData slot
BS.infection.smooth$Infection <- sampleNames(BS.infection.smooth)
BS.infection.smooth$col <- ifelse(
  BS.infection.smooth$Infection == "Control", "blue", "red"
)
pData(BS.infection.smooth)

# Identify the genes associated with "defense response" (top enriched) ----

defRepGenes <- select(org.Bt.eg.db, "GO:0006952", "ENSEMBL", "GOALL")[,"ENSEMBL"]

promTopDensity <- intermediate.5[which(names(intermediate.5) %in% defRepGenes),]
length(promTopDensity)

head(names(promTopDensity), n=10)

# Plot regions ----

ensGenes2plot <- head(names(promTopDensity), n=10)

genes.gr <- readRDS(file.path(genomeDir, "biomart_ensembl_genes.rds"))

immune.gr <- genes.gr[genes.gr$ensembl_gene_id %in% ensGenes2plot]

exons.df <- readRDS(file.path(genomeDir, "biomart_ensembl_exons.rds"))
exons.immune.df <- subset(exons.df, ensembl_gene_id %in% ensGenes2plot)
exons.immune.gr <- GRanges(
  seqnames = exons.immune.df$chromosome_name,
  ranges = IRanges(
    start = exons.immune.df$exon_chrom_start,
    end = exons.immune.df$exon_chrom_end,
    names = exons.immune.df$ensembl_gene_id
  ),
  strand = exons.immune.df$strand,
  external_gene_name = exons.immune.df$external_gene_name
)
exons.immune.grl <- split(exons.immune.gr, names(exons.immune.gr))

for (ensGene in ensGenes2plot){
  plotGR <- immune.gr[ensGene]
  geneName <- plotGR$external_gene_name
  annotGR <- list(
    Promoter = promoters(plotGR, 1500, 500),
    Gene = plotGR,
    Exons = exons.immune.grl[[ensGene]]
  )
  names(annotGR)[2] <- ifelse(geneName == "", ensGene, geneName)
  pdf(
    file = file.path(immuneDir, sprintf(
      "%s-%s_gene_exons_promoter.pdf", geneName, ensGene)),
    width = 9, height = 7
  )
  plotRegion(
    BS.infection.smooth, plotGR, extend = 10E3,
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
}

for (ensGene in ensGenes2plot){
  geneGR <- immune.gr[ensGene]
  promoterGR <- promoters(geneGR, 1500, 500)
  geneName <- geneGR$external_gene_name
  annotGR <- list(
    Promoter = promoterGR,
    Gene = geneGR,
    Exons = exons.immune.grl[[ensGene]]
  )
  names(annotGR)[2] <- ifelse(geneName == "", ensGene, geneName)
  pdf(
    file = file.path(immuneDir, sprintf(
      "%s-%s_promoter.pdf", geneName, ensGene)),
    width = 9, height = 7
  )
  plotRegion(
    BS.infection.smooth, promoterGR, extend = 5E3,
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
}

# Repeat for imprinted genes ----

imprinted.symbol <-
  c("PLAGL1", "SNRPN", "MEST", "PEG10", "GNAS", "NNAT", "NAPIL5")

imprinted.symbol <-
  imprinted.symbol[imprinted.symbol %in% genes.gr$external_gene_name]

imprinted.gr <- genes.gr[genes.gr$external_gene_name %in% imprinted.symbol]
names(imprinted.gr) <- imprinted.gr$external_gene_name

exons.imprinted.df <- subset(exons.df, external_gene_name %in% imprinted.symbol)
exons.imprinted.gr <- GRanges(
  seqnames = exons.df$chromosome_name,
  ranges = IRanges(
    start = exons.df$exon_chrom_start,
    end = exons.df$exon_chrom_end,
    names = exons.df$external_gene_name
  ),
  strand = exons.df$strand
)

exons.imprinted.grl <- split(exons.imprinted.gr, names(exons.imprinted.gr))

for (plotSymbol in names(imprinted.gr)){
  plotGR <- imprinted.gr[plotSymbol]
  annotGR <- list(
    Promoter = promoters(plotGR, 1500, 500),
    Gene = plotGR,
    Exons = exons.imprinted.grl[[plotSymbol]]
  )
  names(annotGR) <- c("Promoter", plotSymbol, "Exons")
  pdf(
    file = file.path(imprintedDir, sprintf(
      "%s_gene_exons_promoter.pdf", plotSymbol)),
    width = 9, height = 7
  )
  plotRegion(
    BS.infection.smooth, plotGR, extend = 10E3,
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
}

for (plotSymbol in names(imprinted.gr)){
  geneGR <- imprinted.gr[plotSymbol]
  promoterGR <- promoters(geneGR, 1500, 500)
  annotGR <- list(
    Promoter = promoterGR,
    Gene = geneGR,
    Exons = exons.imprinted.grl[[plotSymbol]]
  )
  names(annotGR) <- c("Promoter", plotSymbol, "Exons")
  pdf(
    file = file.path(imprintedDir, sprintf(
      "%s_promoter.pdf", plotSymbol)),
    width = 9, height = 7
  )
  plotRegion(
    BS.infection.smooth, promoterGR, extend = 5E3,
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
}
