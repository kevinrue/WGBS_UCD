# Libraries ---------------------------------------------------------------

library(SummarizedExperiment)
library(org.Bt.eg.db)

# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'
genomeDir <- "bostaurus"
plotDir <- "collapsedGroups"
outDir.2kb <- file.path(plotDir, "intermediate_shortlisted_CGI/2kb")
outDir.4kb <- file.path(plotDir, "intermediate_shortlisted_CGI/4kb")

# Create folders ----

if (!dir.exists(outDir.2kb)){
  dir.create(outDir.2kb)
}
if (!dir.exists(outDir.4kb)){
  dir.create(outDir.4kb)
}

# Import promoter information ----

promoters.5 <- readRDS(file.path(outdir, "promoters.5.rds"))

# Import CGI coordinates ----

CGI.gr <- readRDS(file = file.path(outdir, 'CpG.gr.rds'))
seqlevels(CGI.gr) <- gsub("^chr(Un_)?", "", seqlevels(CGI.gr))

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

# Import identifiers of genes to plot ----

library(readxl)
shortlist <- read_excel("~/Desktop/WGBS/WGBS_UCD/exp_data/shortlisted promoters for KR 27-6-17.xlsx")


# Produce plots ----

ensGenes2plot <- shortlist$ensembl_gene_id

all(ensGenes2plot %in% names(genes.gr))
setdiff(ensGenes2plot, names(genes.gr))
which(!ensGenes2plot %in% names(genes.gr))

for (ensGene in ensGenes2plot){
  message(match(ensGene, ensGenes2plot))
  geneGR <- genes.gr[ensGene]
  geneName <- geneGR$external_gene_name
  promoterGR <- promoters(geneGR, upstream = 1500, downstream = 500)
  cgiGR <- subsetByOverlaps(CGI.gr, resize(geneGR, width(geneGR) + 4E3, fix="center"))
  annotGR <- list(
    CGIs = cgiGR,
    Promoter = promoterGR,
    Gene = geneGR,
    Exons = exons.grl[[ensGene]]
  )
  CGcount <- mcols(promoters.5[ensGene])[,"CG.promoter"]
  MethPromoter <- format(
    mcols(promoters.5[ensGene])[,"Meth.promoter"], digits = 2)
  names(annotGR)[3] <- ifelse(geneName == "", ensGene, geneName)
  fileOut.2kb <- file.path(outDir.2kb, sprintf(
    "%i-%s-%s-%s_gene_exons_promoter.pdf",
    CGcount, MethPromoter, geneName, ensGene))
  fileOut.4kb <- file.path(outDir.4kb, sprintf(
    "%i-%s-%s-%s_gene_exons_promoter.pdf",
    CGcount, MethPromoter, geneName, ensGene))
  message(fileOut.2kb)
  pdf(file = fileOut.2kb, width = 9, height = 7)
  plotRegion(
    BS.infection.smooth, promoterGR, extend = 2E3, # 2kb extended
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
  pdf(file = fileOut.4kb, width = 9, height = 7)
  plotRegion(
    BS.infection.smooth, promoterGR, extend = 4E3, # 4kb extended
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
}
