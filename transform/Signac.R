
library(Signac)
library(Seurat)
library(SeuratObject)
library(sceasy)
library(EnsDb.Hsapiens.v86)
library(zellkonverter)

sce <- readH5AD("rna.h5ad")
rna = sce@assays@data@listData[["X"]]
sce2 <- readH5AD("atac.h5ad")
data.atac = sce2@assays@data@listData[["X"]]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(data.atac), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
data.atac <- data.atac[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = data.atac,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC"
)


DefaultAssay(pbmc) = "ATAC"
gene.activities <- GeneActivity(object=pbmc, features = NULL)

scgeneactivity = CreateSeuratObject(counts = gene.activities)
sceasy::convertFormat(scgeneactivity, from="seurat", to="anndata",
                      outFile='scGAM_Signac.h5ad')
