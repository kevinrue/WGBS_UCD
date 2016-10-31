
library(gdata)
library(ggplot2)

RNAseq.file <- "RNAseq/srep13629-s4.xls"
promoters.file <- "bsseq/promoters.5.csv"
outFolder <- "bsseq"

RNAseq.data <- read.xls(RNAseq.file, sheet="Worksheet 1")
promoters.data <- read.csv(promoters.file, row.names = 1)

meta.data <- merge(
  RNAseq.data, promoters.data,
  by = "ensembl_gene_id"
)
names(meta.data)

p.meth.FDR <- ggplot(meta.data, aes(x=Meth.promoter, y=-log10(FDR_24H))) +
  geom_point()
ggsave(
  file.path(outFolder, "meth.prom-vs-log10FDR(point).pdf"),
  p.meth.FDR, width = 7, height = 7
)

p.meth.FDR <- ggplot(meta.data, aes(x=Meth.promoter, y=-log10(FDR_24H))) +
  stat_binhex(bins = 100)
ggsave(
  file.path(outFolder, "meth.prom-vs-log10FDR(binhex).pdf"),
  p.meth.FDR, width = 7, height = 7
)

p.meth.FDR <- ggplot(meta.data, aes(x=Meth.promoter, y=logFC_24H)) +
  geom_point(aes(colour=-log10(FDR_24H)))
ggsave(
  file.path(outFolder, "meth.prom-vs-logFC-log10FDR(point).pdf"),
  p.meth.FDR, width = 7, height = 7
)

p.meth.FDR <- ggplot(meta.data, aes(x=Meth.promoter, y=logFC_24H)) +
  stat_binhex(bins = 100)
ggsave(
  file.path(outFolder, "meth.prom-vs-logFC-log10FDR(binhex).pdf"),
  p.meth.FDR, width = 7, height = 7
)