
# Dependencies ----

library(BSgenome.Btaurus.UCSC.bosTau6)
# names(BSgenome.Btaurus.UCSC.bosTau6)
# Note the "chr" prefix (UCSC-style)
library(biomaRt)
library(bsseq)
library(xlsx)

# Set parameters ----

sequences <- c(
  TNFa = "TAGAGAAGCCCACTCAGAATCCGAGCGGGCGGAGTGTAGGAAGTATCCTTGATGCCTGGGTGTCCCCAACTTTCCAAACCCCCGCCCCCGCGATGGAGAAGAAACCGAGACAGAAGGTGTAGGGCCCGCTACCGCTTCCTCCAGATGAGCTCATGGGTTTCTCCACCAAGGAAGTTTTCCGCTGG",
  
  IL12A = "CAACCAGAGCGCTAGGCTGGTTACTCACTGCGAAGCGGGCACATGCTGAGCGGAGCGGCGGGGACGCGGAACCGAGCCGGCAGTTGGACGCAGACCGGTGCACGCGGCAGGTGAGGGTGGTGGTTGGGAGGCCAAACCAGGGGTCACATTTTTAT",
  
  TLR2 = "GGGGATGCCAGCGGATCCTAATTCCTGACCGACGTACCTGGGACTTGCGCGGCCTTGCAGCGCCTTCCACAGCCTCCGGCCGGGAGCGGCCCGGGAAAAGCGCGGGAACGTGCGCACCCCCTCCTCGCGGGTGCGGGACCGCCGGTTCCGCGGAGTGCGCGTAACCCCTGTGGCCCAGCGCGCCGCCGCGCTTCCCCACGGTCTCCGGCGGGGACCGTGACCCGGGTGCTGCCCGGGTCGGAGGAGGGCGCTGGGGC",
  
  NFKB2 = "CCTGGTGGTGGGAGAGGTGTCGCGACCCGTCCGAGGTGGGTCCGGCCGGGAGAGAATCCTGAACCGGAGCCGCCGCCGCGGTGAGTGGCCGGGTTCAGACCCCTGGGTGGTGGGACACCGGCAAGGGTGGGAGGAGG"
)

outdir <- 'bsseq'

# Prepare sequences in Bioconductor format --------------------------------

DNAset <- DNAStringSet(x = sequences)
DNAset <- append(x = DNAset, values = reverseComplement(DNAset))

# width(DNAset)

# UCSC/BLAT approach ------------------------------------------------------

# https://genome.ucsc.edu/cgi-bin/hgBlat?hgsid=500819951_dRPcuiRElemGq9YK1y9QhKv89cWv&command=start
# Paste the sequences and retrieve the coordinates of each longest match
# which should have a width equal to the full length of the corresponding
# sequence.

# Looks like UCSC/BLAT gives correct coordinates
# Note that the example below (TNFa) is on the complementary strand
TNFa.blat <- GRanges(
  seqnames = "chr23",
  ranges = IRanges(start = 27536746, end = 27536930, names = "TNFa"),
  strand = "-")
TNFa.seq <- getSeq(BSgenome.Btaurus.UCSC.bosTau6, TNFa.blat)
reverseComplement(TNFa.seq) == sequences[["TNFa"]]

IL12A.blat <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = 108291940, end = 108292094, names = "IL12A"),
  strand = "+")
IL12A.seq <- getSeq(BSgenome.Btaurus.UCSC.bosTau6, IL12A.blat)
IL12A.seq == sequences[["IL12A"]]

TLR2.blat <- GRanges(
  seqnames = "chr17",
  ranges = IRanges(start = 3962393, end = 3962649, names = "TLR2"),
  strand = "+")
TLR2.seq <- getSeq(BSgenome.Btaurus.UCSC.bosTau6, TLR2.blat)
TLR2.seq == sequences[["TLR2"]]

NFKB2.blat <- GRanges(
  seqnames = "chr26",
  ranges = IRanges(start = 22890256, end = 22890392, names = "NFKB2"),
  strand = "+")
NFKB2.seq <- getSeq(BSgenome.Btaurus.UCSC.bosTau6, NFKB2.blat)
NFKB2.seq == sequences[["NFKB2"]]


# Import methylation calls ------------------------------------------------

# Import all sites covered by at least one read in at least one sample
BS.rmZero.rds <- readRDS(file = file.path(outdir, "BS.rmZero.rds"))

# GRanges listing all regions of interest
targets.gr <- c(TNFa.blat, IL12A.blat, TLR2.blat, NFKB2.blat)

# Edit from UCSC to Ensembl chromosome naming
seqinfo(targets.gr, new2old=1:length(targets.gr), force=FALSE) <- Seqinfo(
  seqnames = as.character(gsub("chr", "", seqnames(targets.gr))))

# Get methylation % in each region in each sample 
methGenes <- getMeth(
  BSseq = BS.rmZero.rds,
  regions = targets.gr,
  type = "raw", what = "perRegion")
rownames(methGenes) <- names(targets.gr)

# Summarise
summary(t(methGenes))

## One NaN suggests no coverage in the region
## Another 1 suggests a single call, considering that all other values are low

covGenes <- getCoverage(
  BSseq = BS.rmZero.rds, regions = targets.gr, type = "Cov",
  what =  "perRegionTotal")
rownames(covGenes) <- names(targets.gr)

mGenes <- getCoverage(
  BSseq = BS.rmZero.rds, regions = targets.gr, type = "M",
  what =  "perRegionTotal")
rownames(mGenes) <- names(targets.gr)

write.xlsx(
  x = methGenes,
  file = file.path(outdir, "pyroseq.xlsx"),
  sheetName = "Percentage")
write.xlsx(
  x = covGenes,
  file = file.path(outdir, "pyroseq.xlsx"),
  sheetName = "Coverage", append = TRUE)
write.xlsx(
  x = mGenes,
  file = file.path(outdir, "pyroseq.xlsx"),
  sheetName = "Methylated", append = TRUE)
