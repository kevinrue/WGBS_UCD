
library(biomaRt)
library(GenomicRanges)

listMarts(host = "Mar2016.archive.ensembl.org")
mart <- useMart(host = "Mar2016.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")

listDatasets(mart = mart)[,1:2]
mart <- useMart(host = "Mar2016.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "btaurus_gene_ensembl")

listAttributes(mart = mart, page = "feature_page")[,1:2]

umd31_genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "hgnc_symbol", "external_gene_name"), mart = mart)

genes.gr <- GRanges(
  seqnames = umd31_genes$chromosome_name,
  ranges = IRanges(
    start = umd31_genes$start_position,
    end = umd31_genes$end_position,
    names = umd31_genes$ensembl_gene_id),
  strand = umd31_genes$strand)
mcols(genes.gr) <- umd31_genes[,c("hgnc_symbol", "external_gene_name")]

write.csv(x = umd31_genes, file = file.path(outdir, "biomart_ensembl_genes.csv"))
saveRDS(object = genes.gr, file = file.path(outdir, "biomart_ensembl_genes.rds"))
