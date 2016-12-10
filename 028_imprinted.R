
# Packages ----

library(bsseq)

# Settings ----

outDir <- "bsseq"
genomeDir <- "bostaurus"
plotDir <- "imprinted_promoters"

# Create folders ----

if (!dir.exists(plotDir)){
  dir.create(plotDir)
}

# Load methylation calls ----

BS.keepLoci <- readRDS(file = file.path(outDir, "BS.keepLoci.rds"))

BS.keepLoci$col <- c("blue", "red")[as.numeric(BS.keepLoci$Infection)]

# Imprinted genes ----

imprinted.symbol <-
  c("PLAGL1", "SNRPN", "MEST", "PEG10", "GNAS", "NNAT", "NAPIL5")

genes.gr <- readRDS(file.path(genomeDir, "biomart_ensembl_genes.rds"))

imprinted.symbol <-
  imprinted.symbol[imprinted.symbol %in% genes.gr$external_gene_name]

imprinted.gr <- genes.gr[genes.gr$external_gene_name %in% imprinted.symbol]
names(imprinted.gr) <- imprinted.gr$external_gene_name

exons.df <- readRDS(file.path(genomeDir, "biomart_ensembl_exons.rds"))
exons.df <- subset(exons.df, external_gene_name %in% imprinted.symbol)
exons.gr <- GRanges(
  seqnames = exons.df$chromosome_name,
  ranges = IRanges(
    start = exons.df$exon_chrom_start,
    end = exons.df$exon_chrom_end,
    names = exons.df$external_gene_name
  ),
  strand = exons.df$strand
)

exons.grl <- split(exons.gr, names(exons.gr))

for (plotSymbol in names(imprinted.gr)){
  plotGR <- imprinted.gr[plotSymbol]
  annotGR <- list(
    Promoter = promoters(plotGR, 1500, 500),
    Gene = plotGR,
    Exons = exons.grl[[plotSymbol]]
  )
  names(annotGR) <- c("Promoter", plotSymbol, "Exons")
  pdf(
    file = file.path(plotDir, sprintf("%s_gene_exons_promoter.pdf", plotSymbol)),
    width = 9, height = 7
  )
  plotRegion(
    BS.keepLoci, plotGR, extend = 10E3,
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
    Exons = exons.grl[[plotSymbol]]
  )
  names(annotGR) <- c("Promoter", plotSymbol, "Exons")
  pdf(
    file = file.path(plotDir, sprintf("%s_promoter.pdf", plotSymbol)),
    width = 9, height = 7
  )
  plotRegion(
    BS.keepLoci, promoterGR, extend = 5E3,
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
}
