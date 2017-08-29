# Version 1.9.2
# devtools::install_github("Bioconductor-mirror/bsseq", ref = "bd60d13")

library(bsseq)
library(GenomicRanges)
library(BSgenome.Btaurus.UCSC.bosTau6)
library(DelayedArray)

options(DelayedArray.block.size=getOption("DelayedArray.block.size") * 10L)


# Parameters of analysis ----

outdir <- 'bsseq'
genomeDir <- 'bostaurus'

# Load source data ----

# Load CGI coordinates
CGI.gr <- readRDS(file = file.path(outdir, 'CpG.gr.rds'))
seqlevels(CGI.gr) <- gsub("^chr(Un_)?", "", seqlevels(CGI.gr))

# Load transcript coordinates
transcripts.gr <- readRDS(
  file.path(genomeDir, "biomart_ensembl_transcripts.rds"))
seqlevels(transcripts.gr) <- gsub("\\.[[:digit:]]$", "", seqlevels(transcripts.gr))
# head(genes.gr)

# Load gene coordinates
genes.gr <- readRDS(file.path(genomeDir, "biomart_ensembl_genes.rds"))
seqlevels(genes.gr) <- gsub("\\.[[:digit:]]$", "", seqlevels(genes.gr))
# head(genes.gr)

# Load BSgenome for chromosome lengths
genome <- BSgenome.Btaurus.UCSC.bosTau6
seqlevels(genome) <- gsub("chr", "", seqlevels(genome))
seqlevels(genome) <- gsub("^M$", "MT", seqlevels(genome))
#seqlengths(genome)


# Generate intergenic probes ----

# Generate 5 kb tiles of the genome
# Truncate the last tile rather than letting tiles span across chromosomes
tiles <- unlist(tileGenome(
  seqlengths = seqinfo(genome),
  tilewidth = 5E3,
  cut.last.tile.in.chrom = TRUE))
# summary(width(tiles))
# sum(width(tiles) < 5E3)

# Keep tiles over 2 kb from nearest gene (i.e. intergenic)
tilesDistanceToGene <- distanceToNearest(x = tiles, subject = genes.gr)
# head(tilesDistanceToGene)
intergenic.gr <- tiles[which(mcols(tilesDistanceToGene)[,"distance"] > 2E3)]

rm(tilesDistanceToGene, tiles)

# Generate gene body probes ----

tilesGenes <- unlist(tile(x = genes.gr, width = 5E3))

# Generate promoters ranges (wrt. CGI) ----

# Generate promoter coordinates
proms.gr <- promoters(x = transcripts.gr, upstream = 1500, downstream = 500)
# Remove duplicated promoters (shared by multiple transcripts)
proms.gr <- proms.gr[!duplicated(proms.gr)]
# which(duplicated(proms.gr))[1]
# head(proms.gr)

# Promoters overlapping CGI, or not
promOverCGI.bool <- overlapsAny(query = proms.gr, subject = CGI.gr)
promsOverCGI.gr <- proms.gr[promOverCGI.bool]
promsNotOverCGI.gr <- proms.gr[!promOverCGI.bool]
rm(promOverCGI.bool)

# Generate CGI ranges (wrt. promoters) ----

# CGI overlapping promoters, or not
CGIOverProm.bool <- overlapsAny(query = CGI.gr, subject = proms.gr)
CGIOverProm.gr <- CGI.gr[CGIOverProm.bool]
CGINotOverProm.gr <- CGI.gr[!CGIOverProm.bool]
rm(CGIOverProm.bool)

# Load methylation calls ----

BS.collapsed <- readRDS(file.path(outdir, "BS.collapse.all_DelayedMatrix.rds"))

# saveRDS(BS.collapsed, file.path(outdir, "BS.collapsed_DelayedMatrix.rds"))
# BS.infection <- readRDS(file.path(outdir, "BS.infection_DelayedMatrix.rds"))

# Keep only CG >= 5 coverage
BS.2 <- BS.collapsed[
  as.logical(getCoverage(BS.collapsed, what = "perBase") >= 2),
  ]

length(BS.2)
length(BS.2) / length(BS.collapsed) # 99.6% (>=2), 98.3% (>=5)

# Process CGI *not* overlapping promoters ----

# Count loci >= 2 calls in Control
mcols(CGINotOverProm.gr)[,"loci"] <- countOverlaps(
  query = CGINotOverProm.gr, subject = BS.2)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methCGINotOverProm <- getMeth(
  BSseq = BS.collapsed,
  regions = CGINotOverProm.gr,
  type = "raw",
  what = "perRegion")
colnames(methCGINotOverProm) <- "Meth"
values(CGINotOverProm.gr) <- cbind(
  values(CGINotOverProm.gr), DataFrame(as.matrix(methCGINotOverProm)))

rm(methCGINotOverProm)

cgi.noProm <- as.data.frame(CGINotOverProm.gr, row.names = NULL)


# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Non-Promoter CGIs",
    Infection = "Control",
    Value = subset(
      x = values(CGINotOverProm.gr),
      subset = loci >= 10,
      select = "Meth",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Save values for plot (M. bovis)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Non-Promoter CGIs",
    Infection = "M. bovis",
    Value = subset(
      x = values(CGINotOverProm.gr),
      subset = lociMbovis >= 4,
      select = "meth_Mbovis",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

ggData <- readRDS(file = file.path(outdir, "Peat_et_al-ggData.rds"))
