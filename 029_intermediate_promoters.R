
# Libraries ---------------------------------------------------------------

library(SummarizedExperiment)

# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'

# Import previous data ----------------------------------------------------

promoters.5 <- readRDS(file.path(outdir, "promoters.5.rds"))


# Analysis -----------------------------------------------------------------

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

write.csv(tmp, "intermediate_promoters/topCGcontent.csv")
