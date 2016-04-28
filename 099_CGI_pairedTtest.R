# library(broom)
# library(ggbio)

# Count methylations/calls per CpG island ---------------------------------

# seqlevels(CpG.gr)
# # Rename BS.rmZero contigs to match names in CpG.gr
# seqlevels(BS.rmZero) <- gsub(
#   pattern = '^GJ',
#   replacement = 'Un_GJ',
#   x = seqlevels(BS.rmZero))
# seqlevels(BS.rmZero) <- gsub(
#   pattern = '^',
#   replacement = 'chr',
#   x = seqlevels(BS.rmZero))

# Subset both objects to overlapping contigs (save memory)
# seqlevels.intersect <- intersect(
#   seqlevels(BS.rmZero),
#   seqlevels(CpG.gr)
# )
# seqlevels(CpG.gr, force = TRUE) <- seqlevels.intersect
# seqlevels(BS.rmZero, force = TRUE) <- seqlevels.intersect

# For each CpG island, count the number of CG with at least a call
# CpG.gr$CG.covered <- countOverlaps(query = CpG.gr, subject = BS.rmZero)

# CpG.gr[1,]
# BS.rmZero[1,]
# Find/count methylation calls overlapping CpG islands
# subsetByOverlaps(query = BS.rmZero[1,], subject = CpG.gr[1,])
# getCoverage(subsetByOverlaps(query = BS.rmZero[1,], subject = CpG.gr[1,]))
# colSums(getCoverage(subsetByOverlaps(query = BS.rmZero[1,], subject = CpG.gr[1,])))
# getBSseq(subsetByOverlaps(query = BS.rmZero[1,], subject = CpG.gr[1,]), type = "M")
# getMeth(subsetByOverlaps(query = BS.rmZero[1,], subject = CpG.gr[1,]), type = "raw")
# colSums(getMeth(subsetByOverlaps(query = BS.rmZero[1,], subject = CpG.gr[1,])))
# tmp.base <- getMeth(BSseq = BS.rmZero, regions = CpG.gr, type = "raw", what = "perBase")
# tmp.region <- getMeth(BSseq = BS.rmZero, regions = CpG.gr, type = "raw", what = "perRegion")
# 
# Function to test differential methylation level in a paired design
# Pairs of samples must be ordered: group1rep1, group1rep2, ..., group2rep1, group2rep2, ...
# region.paired.t.test <- function(x, factor = rep(1:(length(x)/2), 2)){
#   values1 <- x[1:(length(x)/2)]
#   values2 <- x[(length(x)/2+1):length(x)]
#   if (sum(!is.na(values1+values2)) < 2){
#     return(rep(NA,6))
#   }
#   tidy(t.test(
#     x = values1,
#     y = values2,
#     alternative = "two.sided",
#     paired = TRUE))
# }
# 
# Run the paired t-test on all CpG islands with at least two pairs of samples
# covered by at least one read
# region.paired.t.test(x = tmp.region)
# 
# paired.t.test <- do.call(
#   rbind,
#   apply(X = tmp.region, MARGIN = 1, FUN = region.paired.t.test))
# QQ-plot
# plot(
#   x = -log10(sort(paired.t.test$p.value)),
#   y = -log10(sort(runif(n = sum(!is.na(paired.t.test$p.value)), min = 0, max = 1))))
# rm(paired.t.test)
