
# Libraries ---------------------------------------------------------------

library(biomaRt)
library(GenomicRanges)

# Parameters --------------------------------------------------------------

outdir <- "bostaurus"

if (file.access(outdir) != 0)
  dir.create(outdir)

# Find the right mart -----------------------------------------------------

listMarts(host = "mar2016.archive.ensembl.org")
mart <- useMart(
  host = "mar2016.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")

listDatasets(mart = mart)[,1:2]
mart <- useMart(
  host = "mar2016.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "btaurus_gene_ensembl")

attributePages(mart = mart)
listAttributes(mart = mart, page = "feature_page")[,1:2]

# Import gene information -------------------------------------------------

umd31_genes <- getBM(
  attributes = c(
    "chromosome_name", "start_position", "end_position", "strand",
    "gene_biotype",
    "ensembl_gene_id", "hgnc_symbol", "external_gene_name"),
  mart = mart)

# Good: each ensembl_gene_id is unique with respect to this information
any(table(umd31_genes$ensembl_gene_id) > 1)

# Format as Granges -------------------------------------------------------

genes.gr <- GRanges(
  seqnames = umd31_genes$chromosome_name,
  ranges = IRanges(
    start = umd31_genes$start_position,
    end = umd31_genes$end_position,
    names = umd31_genes$ensembl_gene_id),
  strand = umd31_genes$strand)
mcols(genes.gr) <- umd31_genes[
  ,c("ensembl_gene_id", "hgnc_symbol", "external_gene_name")]

# Export ------------------------------------------------------------------

write.csv(x = umd31_genes, file = file.path(outdir, "biomart_ensembl_genes.csv"))
saveRDS(object = genes.gr, file = file.path(outdir, "biomart_ensembl_genes.rds"))

# Import/export transcript info -------------------------------------------

listAttributes(mart = mart, page = "structure")[,1:2]

umd31_transcripts <- getBM(
  attributes = c(
    "chromosome_name", "transcript_start", "transcript_end", "strand",
    "gene_biotype",
    "transcription_start_site",
    "5_utr_start", "5_utr_end",
    "3_utr_start", "3_utr_end",
    "ensembl_transcript_id", "external_gene_name"),
  mart = mart)

# TODO: multiple transcripts per gene will count multiple times the same CG
any(table(umd31_transcripts$ensembl_transcript_id) > 1)
which(table(umd31_transcripts$ensembl_transcript_id) > 1)
subset(umd31_transcripts, ensembl_transcript_id == "ENSBTAT00000028017")
# For each transcriptID, keep the first one, ordered by increased count of NAs

write.csv(
  x = umd31_transcripts,
  file = file.path(outdir, "biomart_ensembl_transcripts.csv"))

# Format as Granges -------------------------------------------------------

transcripts.gr <- GRanges(
  seqnames = umd31_transcripts$chromosome_name,
  ranges = IRanges(
    start = umd31_transcripts$transcript_start,
    end = umd31_transcripts$transcript_end,
    names = umd31_transcripts$ensembl_transcript_id),
  strand = umd31_transcripts$strand)
mcols(transcripts.gr) <- umd31_transcripts[
  ,c("gene_biotype", "transcription_start_site",
     "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end",
     "external_gene_name")]

saveRDS(
  object = transcripts.gr,
  file = file.path(outdir, "biomart_ensembl_transcripts.rds"))

# Import/export exonic info -------------------------------------------

umd31_exons <- getBM(
  attributes = c(
    "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand",
    "gene_biotype", "is_constitutive",
    "genomic_coding_start", "genomic_coding_end",
    "ensembl_exon_id", "external_gene_name", "ensembl_gene_id"),
  mart = mart)

# For each exonID, keep the first one, ordered by increased count of NAs
# transcripts.leastNA <- by(data = umd31_transcripts, INDICES = umd31_transcripts$ensembl_transcript_id, FUN = function(x){x[which.min(rowSums(is.na(x))),]})

# NA.rows <- apply(X = umd31_transcripts, MARGIN = 1, FUN = function(x){sum(is.na(x))})
# keep <- tapply(NA.rows, umd31_transcripts$ensembl_transcript_id, which.min)

write.csv(
  x = umd31_exons,
  file = file.path(outdir, "biomart_ensembl_exons.csv"))
saveRDS(
  object = umd31_exons,
  file = file.path(outdir, "biomart_ensembl_exons.rds"))



