
# Libraries ---------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(goseq)
library(org.Bt.eg.db)

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

# Keep only CG >= 10 coverage
BS.10 <- BS.collapsed[getCoverage(BS.collapsed) >= 10,]

length(BS.10)
length(BS.10) / length(BS.collapsed) # 95%

# Count CG per gene
genes$CG.gene <- countOverlaps(query = genes, subject = BS.10)
# Count CG per promoter
genes$CG.promoter <- countOverlaps(query = promoters(genes), subject = BS.10)

# Meth % per gene
genes$Meth.gene <- as.vector(getMeth(
  BSseq = BS.10, regions = genes, type = "raw", what = "perRegion"))
# Meth % per promoter
genes$Meth.promoter <- as.vector(getMeth(
  BSseq = BS.10, regions = promoters(genes), type = "raw", what = "perRegion"))

# WARNING: NAs when 0 CG in the region considered
# Subset to genes where both promoter and genes have >= 10 CG

genes.10 <- subset(genes, CG.gene >= 10 & CG.promoter >= 10)

write.csv(x = genes.10, file = file.path(outdir, "genes.prom.10.csv"))

gg.data <- melt(
  data = as.data.frame(mcols(genes.10)),
  id.vars = "ensembl_gene_id",
  measure.vars = c("Meth.gene", "Meth.promoter"),
  variable.name = "State",
  value.name = "Methylation")

gg.data$State = factor(x = gg.data$State, levels = c("Meth.gene", "Meth.promoter"), labels = c("Gene body", "Promoter"))
factor

# Plot 
ggplot(
  data = gg.data,
  mapping = aes(x = State, y = Methylation)) +
  geom_violin(
    draw_quantiles = c(0.25, 0.5, 0.75),
    mapping = aes(fill = State)
  ) +
  theme(
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5))
  )

ggsave(
  filename = file.path(outdir, "Gene-promoter.pdf"),
  width = 9, height = 7)

# Highest and lowest promoter methylation ---------------------------------

length(genes.10)

assayed.genes <- names(genes.10)

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

# Lowest methylated promoters
prom.0_1 <- quickGO(names(subset(genes.10, Meth.promoter < 0.01)))
write.csv(x = prom.0_1, file = file.path(outdir, "prom.0_1.csv"))

prom.0_10 <- quickGO(names(subset(genes.10, Meth.promoter < 0.1)))
write.csv(x = prom.0_10, file = file.path(outdir, "prom.0_10.csv"))

prom.10_20 <- quickGO(names(subset(genes.10, Meth.promoter > 0.1 & Meth.promoter < 0.2)))
write.csv(x = prom.10_20, file = file.path(outdir, "prom.10_20.csv"))

prom.20_30 <- quickGO(names(
  subset(genes.10, Meth.promoter > 0.2 & Meth.promoter < 0.3)))
write.csv(x = prom.10_20, file = file.path(outdir, "prom.20_30.csv"))

# Highest methylated promoters
prom.100_90 <- quickGO(names(subset(genes.10, Meth.promoter > 0.9)))
write.csv(x = prom.100_90, file = file.path(outdir, "prom.100_90.csv"))

prom.33_66 <- quickGO(names(
  subset(genes.10, Meth.promoter > 1/3 & Meth.promoter < 2/3)))
write.csv(x = prom.33_66, file = file.path(outdir, "prom.33_66.csv"))
