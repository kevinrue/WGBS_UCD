# Libraries ---------------------------------------------------------------

library(SummarizedExperiment)
library(org.Bt.eg.db)

options(DelayedArray.block.size=45E7)

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

promoters.10 <- readRDS(file.path(outdir, "promoters.10.rds"))

# Import CGI coordinates ----

CGI.gr <- readRDS(file = file.path(outdir, 'CpG.gr.rds'))
seqlevels(CGI.gr) <- gsub("^chr(Un_)?", "", seqlevels(CGI.gr))

# Import smoothed methylation collapsed by group ----

# BS.infection.smooth <- readRDS(file.path(outdir, "BS.infection.smooth.rds"))
# BS.infection.smooth <- updateObject(BS.infection.smooth)
# saveRDS(BS.infection.smooth, file.path(outdir, "BS.infection.smooth_DelayedMatrix.rds"))
BS.infection.smooth <- readRDS(file.path(outdir, "BS.infection.smooth_DelayedMatrix.rds"))

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

# Import sequences ----

sequences <- DNAStringSet(c(
  TNF = "TAGAGAAGCCCACTCAGAATCCGAGCGGGCGGAGTGTAGGAAGTATCCTTGATGCCTGGGTGTCCCCAACTTTCCAAACCCCCGCCCCCGCGATGGAGAAGAAACCGAGACAGAAGGTGTAGGGCCCGCTACCGCTTCCTCCAGATGAGCTCATGGGTTTCTCCACCAAGGAAGTTTTCCGCTGG",
  IL12A = "CAACCAGAGCGCTAGGCTGGTTACTCACTGCGAAGCGGGCACATGCTGAGCGGAGCGGCGGGGACGCGGAACCGAGCCGGCAGTTGGACGCAGACCGGTGCACGCGGCAGGTGAGGGTGGTGGTTGGGAGGCCAAACCAGGGGTCACATTTTTAT",
  TLR2 = "GGGGATGCCAGCGGATCCTAATTCCTGACCGACGTACCTGGGACTTGCGCGGCCTTGCAGCGCCTTCCACAGCCTCCGGCCGGGAGCGGCCCGGGAAAAGCGCGGGAACGTGCGCACCCCCTCCTCGCGGGTGCGGGACCGCCGGTTCCGCGGAGTGCGCGTAACCCCTGTGGCCCAGCGCGCCGCCGCGCTTCCCCACGGTCTCCGGCGGGGACCGTGACCCGGGTGCTGCCCGGGTCGGAGGAGGGCGCTGGGGC",
  NFKB2 = "CCTGGTGGTGGGAGAGGTGTCGCGACCCGTCCGAGGTGGGTCCGGCCGGGAGAGAATCCTGAACCGGAGCCGCCGCCGCGGTGAGTGGCCGGGTTCAGACCCCTGGGTGGTGGGACACCGGCAAGGGTGGGAGGAGG"
))

# Test: find coordinates of sequence in genome ----

library(BSgenome.Btaurus.UCSC.bosTau6)

# TNF
gr <- matchPattern(reverseComplement(sequences[["TNF"]]), BSgenome.Btaurus.UCSC.bosTau6[["chr23"]], fixed = "subject")
length(gr)
TNF.pcr <- GRanges("23", IRanges(start(gr), end(gr), names = "TNF"), strand = "-")

# IL12A
# weird, Ensembl says rev strand, did Alan give me the rc?
gr <- matchPattern(sequences[["IL12A"]], BSgenome.Btaurus.UCSC.bosTau6[["chr1"]], fixed = "subject")
length(gr)
IL12A.pcr <- GRanges("1", IRanges(start(gr), end(gr), names = "IL12A"), strand = "+")

# TLR2
# weird, Ensembl says rev strand, did Alan give me the rc?
gr <- matchPattern(sequences[["TLR2"]], BSgenome.Btaurus.UCSC.bosTau6[["chr17"]], fixed = "subject")
length(gr)
TLR2.pcr <- GRanges("17", IRanges(start(gr), end(gr), names = "TLR2"), strand = "+")

# NFKB2
gr <- matchPattern(sequences[["NFKB2"]], BSgenome.Btaurus.UCSC.bosTau6[["chr26"]], fixed = "subject")
length(gr)
NFKB2.pcr <- GRanges("26", IRanges(start(gr), end(gr), names = "NFKB2"), strand = "+")

PCR.gr  <- c(TNF.pcr, IL12A.pcr, TLR2.pcr, NFKB2.pcr)

# Produce plots ----

ensGenes2plot <- with(
  genes.gr,
  data.frame(
    gene_name = c("TNF","IL12A","TLR2","NFKB2"),
    gene_id = ensembl_gene_id[match(c("TNF","IL12A","TLR2","NFKB2"), hgnc_symbol)],
    stringsAsFactors = FALSE
  )
)

all(ensGenes2plot$gene_id %in% names(genes.gr))
setdiff(ensGenes2plot$gene_id, names(genes.gr))
which(!ensGenes2plot$gene_id %in% names(genes.gr))

for (ensGene in ensGenes2plot$gene_id){
  message(match(ensGene, ensGenes2plot$gene_id))
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
  CGcount <- mcols(promoters.10[ensGene])[,"CG.promoter"]
  MethPromoter <- format(
    mcols(promoters.10[ensGene])[,"Meth.promoter"], digits = 2)
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
