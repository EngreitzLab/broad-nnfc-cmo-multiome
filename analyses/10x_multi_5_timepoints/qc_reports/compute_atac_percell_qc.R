## Per-cell ATAC QC metrics for one pseudobulk (timepoint): fragment count, TSS enrichment and
## nucleosome signal, from that pseudobulk's fragments file. Every barcode in the file is kept
## (no cell-calling filter): the fragments are already restricted to the pseudobulk's cells.
##
## Usage:
##   Rscript compute_atac_percell_qc.R <day> <fragments.tsv.gz> <gene_bounds.bed> <out.tsv.gz> [work_dir]
##
##   day              sample name used by ArchR (e.g. d0)
##   fragments        pseudobulk fragments (chr, start, end, barcode, count), bgzip-compressed
##   gene_bounds.bed  gene bounds BED6 (chr, start, end, symbol, score, strand); the analysis used
##                    scE2G resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.bed
##   out              output TSV.gz: barcode, atac_fragment_count, atac_tss_enrichment,
##                    atac_nucleosome_signal
##   work_dir         ArchR working directory (default: <dirname(out)>/archr_work_<day>)
##
## Requires R with ArchR (1.0.3 was used), BSgenome.Hsapiens.UCSC.hg38 and data.table.
suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: Rscript compute_atac_percell_qc.R <day> <fragments.tsv.gz> <gene_bounds.bed> <out.tsv.gz> [work_dir]")
}
day <- args[1]
frag_file <- normalizePath(args[2])
gene_bed <- normalizePath(args[3])
out_file <- file.path(normalizePath(dirname(args[4])), basename(args[4]))
work_dir <- if (length(args) >= 5) args[5] else file.path(dirname(out_file), paste0("archr_work_", day))

message("[", day, "] building gene/genome annotation from ", gene_bed)
genes_dt <- fread(gene_bed)
setnames(genes_dt, 1:6, c("chr", "start", "end", "symbol", "score", "strand"))
genes_dt <- genes_dt[!duplicated(symbol)]
genes_gr <- GRanges(genes_dt$chr, IRanges(genes_dt$start + 1, genes_dt$end), strand = genes_dt$strand, symbol = genes_dt$symbol)

tss_gr <- resize(genes_gr, width = 1, fix = "start")  # TSS point per gene, strand-aware

genome_annotation <- createGenomeAnnotation(genome = BSgenome.Hsapiens.UCSC.hg38)
gene_annotation <- createGeneAnnotation(TSS = tss_gr, exons = genes_gr, genes = genes_gr)

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

message("[", day, "] createArrowFiles (minTSS=0, minFrags=0: no cell-calling filter, keep every barcode)")
arrow <- createArrowFiles(
  inputFiles = frag_file,
  sampleNames = day,
  geneAnnotation = gene_annotation,
  genomeAnnotation = genome_annotation,
  minTSS = 0,
  minFrags = 0,
  addTileMat = FALSE,
  addGeneScoreMat = FALSE,
  force = TRUE,
  threads = 1
)

message("[", day, "] extracting per-cell QC metadata (nFrags, TSSEnrichment) from ArchR")
proj <- ArchRProject(
  ArrowFiles = arrow,
  outputDirectory = file.path(work_dir, "proj"),
  geneAnnotation = gene_annotation,
  genomeAnnotation = genome_annotation,
  copyArrows = FALSE
)
cd <- as.data.frame(getCellColData(proj))
cd <- data.table(barcode = sub("^[^#]+#", "", rownames(cd)), atac_fragment_count = cd$nFrags, atac_tss_enrichment = cd$TSSEnrichment)

message("[", day, "] computing nucleosome signal from fragment lengths (ArchR has no per-cell metric for this)")
# mono-nucleosomal (147-294 bp) vs nucleosome-free (<147 bp) fragment count ratio, the same
# convention as Signac's NucleosomeSignal / 10x's fragment-size QC.
frag <- fread(frag_file, header = FALSE, col.names = c("chr", "start", "end", "barcode", "count"))
frag[, flen := end - start]
nuc <- frag[, .(
  mono = sum(flen >= 147 & flen <= 294),
  nfr = sum(flen < 147)
), by = barcode]
nuc[, atac_nucleosome_signal := ifelse(nfr > 0, mono / nfr, NA_real_)]

out <- merge(cd, nuc[, .(barcode, atac_nucleosome_signal)], by = "barcode", all.x = TRUE)
fwrite(out, out_file, sep = "\t", compress = "gzip")
message("[", day, "] wrote ", nrow(out), " barcodes -> ", out_file)
