
# Libraries ---------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(goseq)
library(org.Bt.eg.db)
library(BS)

# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'
genomeDir <- 'bostaurus'

# Import previous data ----------------------------------------------------

BS <- readRDS(file.path(outdir, "BS.rmZero.rds"))

sampleNames(BS)
colnames(BS)
colnames(getCoverage(BS))
rownames(colData(BS))
colData(BS)

genes <- readRDS(file.path(genomeDir, "biomart_ensembl_genes.rds"))

# Compare gene body to promoter -------------------------------------------

# Combine all samples
BS.collapsed <- collapseBSseq(BSseq = BS, columns = rep("Merged", 16))

# Keep only CG >= 2 coverage
BS.2 <- BS.collapsed[getCoverage(BS.collapsed, what = "perBase") >= 2,]

length(BS.2)
length(BS.2) / length(BS.collapsed) # 99.5%

# Count CG per promoter
genes$CG.promoter <- countOverlaps(query = promoters(genes), subject = BS.2)

# Meth % per promoter
genes$Meth.promoter <- as.vector(getMeth(
  BSseq = BS.2, regions = promoters(genes), type = "raw", what = "perRegion"))

# WARNING: NAs when 0 CG in the region considered
# Subset to genes where both promoter and genes have >= 10 CG

promoters.5 <- subset(genes, CG.promoter >= 5)

saveRDS(promoters.5, file = file.path(outdir, "promoters.5.rds"))
write.csv(x = promoters.5, file = file.path(outdir, "promoters.5.csv"))

gg.data <- melt(
  data = as.data.frame(mcols(promoters.5)),
  id.vars = "ensembl_gene_id",
  measure.vars = c("Meth.promoter"),
  variable.name = "Region",
  value.name = "Methylation")

gg.data$Region = factor(
  x = gg.data$Region,
  levels = c("Meth.promoter"),
  labels = c("Promoter"))


# Plot 
ggplot(
  data = gg.data,
  mapping = aes(x = Region, y = Methylation)) +
  geom_violin(
    draw_quantiles = c(0.25, 0.5, 0.75),
    mapping = aes(fill = Region)
  ) +
  theme(
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5))
  )

ggsave(
  filename = file.path(outdir, "Promoters.pdf"),
  width = 9, height = 7)

# Highest and lowest promoter methylation ---------------------------------

length(promoters.5)

assayed.genes <- names(promoters.5)

quickGO <- function(target){
  gene.vector = as.integer(assayed.genes %in% target)
  names(gene.vector) = assayed.genes
  
  pwf = data.frame(
    DEgenes = gene.vector,
    bias.data = rep(1, length(gene.vector)),
    pwf = rep(1, length(gene.vector)))
  
  GO.nobias = goseq(pwf, "bosTau6", "ensGene", method="Hypergeometric")
  
  GO.nobias$over_rep_FDR <- p.adjust(
    p = GO.nobias$over_represented_pvalue, method = "BH")
  
  GO.nobias$under_rep_FDR <- p.adjust(
    p = GO.nobias$under_represented_pvalue, method = "BH")
  
  GO.nobias
  }

countGenes = data.frame(Methylation = character(), Count = integer())

# Lowest methylated promoters

subGr <- subset(promoters.5, Meth.promoter < 0.01)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "0 - 0.01", Count = length(subGr)))
prom.0_1 <- quickGO(names(subGr))
write.csv(x = prom.0_1, file = file.path(outdir, "prom.0_1.csv"))

subGr <- subset(promoters.5, Meth.promoter > 0.01 & Meth.promoter < 0.1)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "0.01 - 0.1", Count = length(subGr)))
prom.1_10 <- quickGO(names(subGr))
write.csv(x = prom.1_10, file = file.path(outdir, "prom.1_10.csv"))

subGr <- subset(promoters.5, Meth.promoter > 0.1 & Meth.promoter < 0.2)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "0.1 - 0.2", Count = length(subGr)))
prom.10_20 <- quickGO(names(subGr))
write.csv(x = prom.10_20, file = file.path(outdir, "prom.10_20.csv"))

subGr <- subset(promoters.5, Meth.promoter > 0.2 & Meth.promoter < 1/3)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "0.2 - 1/3", Count = length(subGr)))
prom.20_33 <- quickGO(names(subGr))
write.csv(x = prom.20_33, file = file.path(outdir, "prom.20_33.csv"))

# Mid-meth promoters

subGr <- subset(promoters.5, Meth.promoter > 1/3 & Meth.promoter < 2/3)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "1/3 - 2/3", Count = length(subGr)))
prom.33_66 <- quickGO(names(subGr))
write.csv(x = prom.33_66, file = file.path(outdir, "prom.33_66.csv"))

# Highest methylated promoters

subGr <- subset(promoters.5, Meth.promoter > 2/3 & Meth.promoter < 0.8)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "2/3 - 0.8", Count = length(subGr)))
prom.66_80 <- quickGO(names(subGr))
write.csv(x = prom.66_80, file = file.path(outdir, "prom.66_80.csv"))

subGr <- subset(promoters.5, Meth.promoter > 0.8 & Meth.promoter < 0.9)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "0.8 - 0.9", Count = length(subGr)))
prom.80_90 <- quickGO(names(subGr))
write.csv(x = prom.80_90, file = file.path(outdir, "prom.80_90.csv"))

subGr <- subset(promoters.5, Meth.promoter > 0.9 & Meth.promoter < 0.99)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "0.9 - 0.99", Count = length(subGr)))
prom.90_99 <- quickGO(names(subGr))
write.csv(x = prom.90_99, file = file.path(outdir, "prom.90_99.csv"))

subGr <- subset(promoters.5, Meth.promoter > 0.99)
countGenes <- rbind(
  countGenes, data.frame(Methylation = "0.99 - 1", Count = length(subGr)))
prom.99_100 <- quickGO(names(subset(promoters.5, Meth.promoter > 0.9)))
write.csv(x = prom.99_100, file = file.path(outdir, "prom.99_100.csv"))

write.csv(
  x = countGenes,
  file = file.path(outdir, "Promoters_countGenes.csv"),
  row.names = FALSE)