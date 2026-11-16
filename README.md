# Single-cell-ECCITE-Seq-IEL
# ECCITE-seq data processing (RNA + ADT + HTO + TCR)
R workflows for single-cell ECCITE-Seq (RNA+VDJ+ADT) on intraepithelial lymphocytes of mice ileum under SD and WD. This includes data processing, integration, and downstream analyses.
Per-sample processing of single-cell ECCITE-seq data from mouse gut lymphocytes: gene expression (RNA), surface protein (ADT), sample hashtags (HTO) and TCR α/β V(D)J. This repository currently contains two samples, ielsw1 and ielsw2.

Pipeline

ECCITE_Seq_data_processing.R runs the same steps on each sample:

Load Cell Ranger RNA/ADT matrices and V(D)J contig/clonotype tables.
Build the Seurat object with RNA, ADT (CLR-normalized) and HTO assays. No ambient RNA correction is applied.
QC and filtering: cells with < 200 genes, > 15% mitochondrial reads, or > 75,000 UMIs are removed. Ribosomal and hemoglobin percentages are computed for QC only.
Hashtag demultiplexing: HTODemux (positive quantile 0.85 and 0.99) and MULTIseqDemux are compared. MULTIseqDemux calls are used for the final assignment; negative and doublet cells are removed.
RNA doublet removal with scDblFinder (expected rate 7.5%).
Normalization, PCA and UMAP on RNA only (30 PCs).
TCR processing: productive contigs are summarized per cell (chains detected) and joined with clonotype information.
Input layout
data/inputs/
├── <sample>_rna_adt/outs/filtered_feature_bc_matrix/
└── <sample>_vdj/outs/
    ├── filtered_contig_annotations.csv
    └── clonotypes.csv
Output

Written to data/outputs/ for each sample:

File	Content
<sample>_processed.rds	Filtered Seurat object (RNA, ADT, HTO) with demultiplexing and doublet calls in the metadata
<sample>_vdj_summary.csv	One row per cell: number of chains and TRA/TRB/TRG/TRD detection
<sample>_vdj_full.csv	One row per productive contig, with clonotype information
<sample>_QC.pdf	QC, HTO demultiplexing and UMAP plots

barcode_sample in the VDJ tables is <sample>_<barcode>, matching cell names after merging samples with add.cell.ids.

Requirements

R with Seurat (v5), dplyr, ggplot2, patchwork, SingleCellExperiment and scDblFinder.

Usage

Set data_dir, out_dir and samples at the top of the script, then run:

r
source("ECCITE_Seq_data_processing.R")
Design notes

Ambient RNA and mitochondrial filter

SoupX is not used. In a clean whole-cell dataset like this one, it can over-correct and remove real signal.
The mitochondrial cutoff is 15% instead of the more common 20%.

Hashtag demultiplexing

A high fraction of negative cells suggests insufficient hashtag labeling, poor antibody conjugation, or dead cells.
HTODemux applies one global cutoff at a chosen quantile of each hashtag's count distribution. Cells above it are called positive.
MULTIseqDemux sweeps thresholds and can rescue weakly tagged cells, so it is used for the final calls.

Embedding

UMAP is RNA-only. A WNN (RNA + ADT) embedding is not used because limited ADT granularity and sequencing-depth differences, which CLR normalization corrects poorly, add noise.

TCR

Only TCR α/β was amplified and sequenced, so γδ T cells are not expected to carry VDJ calls. TCR-negative cells are expected to include γδ T cells and CD8αα IELs.
γδ T cells can still be identified from RNA and from ADT (the Universal TotalSeq-C cocktail includes anti-mouse TCR γ/δ; feature names TCR-gd / TCRgd).
The junction is the CDR3 region: the recombined sequence between the V and J segments (plus D, if present). CDJ is CDR3 plus one flanking amino acid on each side (last FWR3 residue + CDR3 + first FWR4 residue).

Downstream conventions (applied in later analysis, not in this script)

CDR3 calls are kept only if supported by > 2 UMIs and > 5 reads.
Clonotypes may be defined by CDR3 or CDJ.
A clonotype is called expanded when ≥ 2 cells share the same TCR. scRNA-seq samples only a small part of the repertoire, so a TCR seen in two cells is already likely clonal, and a clone making up 1% of cells is considered highly expanded.
Further reading
https://pubmed.ncbi.nlm.nih.gov/35320704/
https://pmc.ncbi.nlm.nih.gov/articles/PMC11608099/
https://pmc.ncbi.nlm.nih.gov/articles/PMC7677182/
