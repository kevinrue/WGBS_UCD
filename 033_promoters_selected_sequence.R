# Libraries ---------------------------------------------------------------

library(SummarizedExperiment)
library(org.Bt.eg.db)

# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'
genomeDir <- "bostaurus"
plotDir <- "collapsedGroups"
outDir <- file.path(plotDir, "sequences")

# Create folders ----

if (!dir.exists(outDir)){
  dir.create(outDir)
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

# Import sequences ----

sequences <- DNAStringSet(c(
  NOS2 = "CTGCAGGGCCAGGCAGATGGTATTCCGATGACCAGCGGGTTCCTCGCCACCTTGCTTCAGGTCAACCGTGAAGCTCGCCAGGTCCCATTCAAGAGGCTCACAGGGTTTGTGCTTTGTCCTTATTGAGTTATGTGTTCTTCCTTTAGACAGTTACGTTGTTTTGTAGAATTAAATCCATAGGATCCCCGCATGCAACTCCAAAAGAGGAACAGAGGTCGCAAAACCCCGATGAAGTGGCAGTTCATGACTTACGAGCGTTATCCCCTGCACTGATGCCGCGGTGTGATAAACCAACCACACTGAGAAG",
  C1QB = "AGAATCTGAATCAGGGTTTCTGACTGCAAACCCAGCCCTCGGGCAGCAGGCGGGCCACATCTCCACTTCGTGGAGGAACACCGCAGTTCCTGCAGAAGTTCATGAGGTCAGGGCCCTTGGCAGGCGCTCTCTGCCACCACTGTGGACATTCCCCTGGAACACCCTACTCCCCACAACAGCTGCTGAGCTGGCGAAACCAAACTAGCCCCCCACCCATCCTATTTACAGTAAACCCCTCCCGATTGCAGAAATGGGACCTGAAAGTGCCT",
  IL2RA = "CCAGGGCATCATGGTGAGAGAATTAAGGCTAGTACTGATTGTAAGTCTTGGGAGAATATTCAGGAGCTCTGAGCGTGGGGCAGTCTGTGCTTCCTCTTCAATTATCTGAATATTGGAGACCCCACCCGAAGACCCAGTGGACTTCTTCAAGGGCTCCGTGGTCAGACAGCTTTAGGAAATTCCCATCTAGGGCCAAGACTGCTTGGCATGACCTGGCGCCCCTCCATCTTGAGTGGCACCTGCCGACTCCCCGGGGAGCTGCGAGTCCGGGACCTCAGTTTCACGGTGACATCAGATTATGAAACTCTAGTTGAGACCACTGCCAAGAAGGGCTTGCTCACCCTGCTTTAACGGCAGTGGGAATCTCCCTGTCTTTTT"
))

# Test: find coordinates of sequence in genome ----

library(BSgenome.Btaurus.UCSC.bosTau6)

# NOS2
gr <- matchPattern(reverseComplement(sequences[["NOS2"]]), BSgenome.Btaurus.UCSC.bosTau6[["chr19"]], fixed = "subject")
NOS2.pcr <- GRanges("19", IRanges(start(gr), end(gr), names = "NOS2"), strand = "-")

# C1QB
gr <- matchPattern(reverseComplement(sequences[["C1QB"]]), BSgenome.Btaurus.UCSC.bosTau6[["chr2"]], fixed = "subject")
C1QB.pcr <- GRanges("2", IRanges(start(gr), end(gr), names = "C1QB"), strand = "-")

# IL2RA
gr <- matchPattern(sequences[["IL2RA"]], BSgenome.Btaurus.UCSC.bosTau6[["chr13"]], fixed = "subject")
IL2RA.pcr <- GRanges("13", IRanges(start(gr), end(gr), names = "IL2RA"), strand = "+")

PCR.gr  <- c(NOS2.pcr, C1QB.pcr, IL2RA.pcr)

# Produce plots ----

ensGenes2plot <- subset(
  shortlist,
  external_gene_name %in% c("NOS2", "C1QB", "IL2RA")
  )$ensembl_gene_id

all(ensGenes2plot %in% names(genes.gr))
setdiff(ensGenes2plot, names(genes.gr))
which(!ensGenes2plot %in% names(genes.gr))

for (ensGene in ensGenes2plot){
  message(match(ensGene, ensGenes2plot))
  geneGR <- genes.gr[ensGene]
  geneName <- geneGR$external_gene_name
  pcrGR <- PCR.gr[geneName]
  promoterGR <- promoters(geneGR, upstream = 1500, downstream = 500)
  cgiGR <- subsetByOverlaps(CGI.gr, resize(geneGR, width(geneGR) + 4E3, fix="center"))
  annotGR <- list(
    PCR = pcrGR,
    CGIs = cgiGR,
    Promoter = promoterGR,
    Gene = geneGR,
    Exons = exons.grl[[ensGene]]
  )
  CGcount <- mcols(promoters.5[ensGene])[,"CG.promoter"]
  MethPromoter <- format(
    mcols(promoters.5[ensGene])[,"Meth.promoter"], digits = 2)
  names(annotGR)[4] <- ifelse(geneName == "", ensGene, geneName)
  fileOut <- file.path(outDir, sprintf(
    "%i-%s-%s-%s_gene_exons_promoter.pdf",
    CGcount, MethPromoter, geneName, ensGene))
  message(fileOut)
  pdf(file = fileOut, width = 9, height = 7)
  plotRegion(
    BS.infection.smooth, promoterGR, extend = 4E3, # 4kb extended
    annoTrack = annotGR,
    stat.ylim = 0:1
  )
  dev.off()
}
