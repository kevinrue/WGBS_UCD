#!/usr/local/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(stringr))

option_list <- list(
  make_option(
    c('-o', '--outfile'), type='character', default=NULL,
    help='Output file [input.txt -> input.csv]'
  )
  )

parser.description = '
Formats the output of 003_fastqc_summarisation.sh for convenient
analysis of trim_galore trimmed reads in Excel Pivot tables.

Mandatory arguments:
\tinput_file
\t\tFile produced by the script 003_fastqc_summarisation.sh
'

# Define the parser
parser <- OptionParser(
  usage = "%prog [options] input_file",
  option_list = option_list,
  description = parser.description
)

# Parse the CL arguments
args = parse_args(parser, positional_arguments = 1)

file.summary <- args$args[1]

if (file.access(file.summary) != 0){
  stop('File not found: ', file.summary)
}

if (is.null(args$options$outfile)){
  outfile <- gsub('txt$', 'csv', file.summary)
} else {
  if (file.access(dirname(args$options$outfile)) != 0){
    stop('Folder not found: ', dirname(args$options$outfile))
  }
  outfile <- args$options$outfile
}

fastqc <- read.table(
  file = file.summary,
  sep = '\t',
  col.names = c('Value', 'QC', 'File', 'Batch'))

files <- as.character(fastqc$File)

filesplit <- strsplit2(files, '_')

fastqc$Lane <- apply(filesplit, 1, function(x){x[grep('^L00', x)]})

fastqc$Sample <- strsplit2(files, '_')[,1]

fastqc$Read <- factor(
  str_detect(files, 'R1'),
  levels = c(TRUE, FALSE),
  labels = c('R1', 'R2'))

fastqc$Treatment <- factor(
  str_detect(files, 'NOT_BS'),
  levels = c(TRUE, FALSE),
  labels = c('NOT BS', 'BS'))

fastqc$Infection <- factor(
  str_detect(files, '^C'),
  levels = c(TRUE, FALSE),
  labels = c('Control', 'M. bovis'))

fastqc$Unpaired <- factor(
  str_detect(files, 'unpaired'),
  levels = c(TRUE, FALSE),
  labels = c('Unpaired', 'Paired'))

write.csv(x = fastqc, file = outfile, row.names = FALSE)
