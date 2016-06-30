
# Libraries ---------------------------------------------------------------

library(ggplot2)
library(broom)

# Variables ---------------------------------------------------------------

outdir <- "bsseq"

# Import objects ----------------------------------------------------------

CGI.gr <- readRDS(file.path(outdir, "CpG.gr.rds"))
BS.unstranded <- readRDS(file.path(outdir, "BS.unstranded.rds"))

# Match seqlevels ---------------------------------------------------------

seqlevels(CGI.gr) <- gsub("^chr(Un_)?", "", seqlevels(CGI.gr))
seqlevels(BS.unstranded) <- gsub("(\\.1)?", "", seqlevels(BS.unstranded))

# Basic stats per CGI -----------------------------------------------------

# Count calls in each CGI
mcols(CGI.gr)[,"TotalCov"] <- rowSums(getCoverage(
  BSseq = BS.unstranded, regions = CGI.gr, what = "perRegionTotal"))

# Add up coverage in each CGI
mcols(CGI.gr)[,"AvgCov"] <- rowSums(getCoverage(
  BSseq = BS.unstranded, regions = CGI.gr, what = "perRegionAverage"))
summary(mcols(CGI.gr)[,"AvgCov"])

saveRDS(object = CGI.gr, file = file.path(outdir, "CGI.gr.rds"))

mcols(CGI.gr)[,"AvgMeth"] <- rowMeans(getMeth(
  BSseq = BS.unstranded, type = "raw", regions = CGI.gr, what = "perRegion"))
summary(mcols(CGI.gr)[,"AvgMeth"])

# Plot the average methylation level of CGI -------------------------------

# (across all samples)
ggplot(as.data.frame(mcols(CGI.gr))) +
  stat_bin(mapping = aes(x = AvgMeth), bins = 20) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 10E3))
  
ggsave(filename = "CGI_avgMeth.pdf", path = outdir, width = 8, height = 6)

# Per sample CGI stats ----------------------------------------------------

meth.CGI <- getMeth(
  BSseq = BS.unstranded, type = "raw", regions = CGI.gr, what = "perRegion")

shoreLeft <- flank(x = CGI.gr, width = width(CGI.gr), start = TRUE)
shoreRight <- flank(x = CGI.gr, width = width(CGI.gr), start = FALSE)

covShoreLeft <- getCoverage(BSseq = BS.unstranded, regions = shoreLeft, type = "Cov", what = "perRegionTotal")
covShoreRight <- getCoverage(BSseq = BS.unstranded, regions = shoreRight, type = "Cov", what = "perRegionTotal")
covShores <- covShoreLeft + covShoreRight
rm(covShoreLeft, covShoreRight)

mShoreLeft <- getCoverage(BSseq = BS.unstranded, regions = shoreLeft, type = "M", what = "perRegionTotal")
mShoreRight <- getCoverage(BSseq = BS.unstranded, regions = shoreRight, type = "M", what = "perRegionTotal")
mShores <- mShoreLeft + mShoreRight
rm(mShoreLeft, mShoreRight)

meth.shores <- mShores / covShores
rm(mShores, shoreLeft, shoreRight)
# Keep covShore, to filter shores with too little coverage

# Compare total coverage of CGI and shores
summary(CGI.gr$TotalCov)
summary(rowSums(covShores))

# Paired t-test comparing the Meth% of CG and shores
meth.diff <- meth.CGI - meth.shores
# meth.diff

tTestNaN <- function(x){
  if(sum(!is.nan(x) & !is.na(x)) < 2){
    return(NA)
  }
  tidy(t.test(
    x = x))
}
ttests <- do.call(rbind, apply(X = meth.diff, MARGIN = 1, FUN = tTestNaN))
ttests$BH <- p.adjust(p = ttests$p.value, method = "BH")

mcols(CGI.gr) <- cbind(mcols(CGI.gr), ttests)

# Volcano plot
ggplot(as.data.frame(mcols(CGI.gr))[which(!is.na(mcols(CGI.gr)[,"p.value"])),]) +
  geom_point(
    mapping = aes(
      x = estimate,
      y = -log10(p.value),
      colour = (BH <= 0.05 & abs(estimate) > 0.1))) +
  scale_x_continuous(limits = c(-1,1))

# QQ plot
# ggplot() +
#   geom_point(
#     data = data.frame(
#       Expected = sort(-log10(ppoints(ttests$p.value))),
#       Observed = sort(-log10(ttests$p.value))
#     ),
#     mapping = aes(x = Expected, y = Observed))

