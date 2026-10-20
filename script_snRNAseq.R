
####################################################################
# Single-nucleus RNA-seq (snRNA-seq) Analysis Pipeline
# Project: Perilesional tissue of mouse GG/Ctrl at p5 and p30
# Samples: 2 samples per group. Each sample containes nuclei from 3 different mice. 

# 
# This R script performs preprocessing, integration, and analysis of 
# snRNA-seq data using Seurat and related packages. 
#
# Author: Silvia Cases Cunillera


####################################################################
# ------------------------- #
# 0. Install & Load Packages
# ------------------------- #

# CRAN packages
install.packages("devtools")
# Load core libraries
library(Seurat)
library(SeuratObject)
library(hdf5r)
library(patchwork)
library(stringr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(EnhancedVolcano)
library(tidyr)





# ================================
# 1. Setup sample metadata
# ================================
base_path <- "path_to_cell_Ranger_outs"

samples <- data.frame(
  group = c(
    "Mouse_GG_P5",
    "Mouse_GG_P5",
    "Mouse_Control_P5",
    "Mouse_Control_P5",
    "Mouse_GG_P30",
    "Mouse_GG_P30",
    "Mouse_Control_P30",
    "Mouse_Control_P30"
  ),
  sample_id = c(
    "SC142-1",
    "SC141-4",
    "SC142-5",
    "SC141-10",
    "SC145-3",
    "SC134-1-4",
    "SC134-1-1",
    "SC134-2-2"
  ),
  timepoint = c(
    "P5","P5","P5","P5",
    "P30","P30","P30","P30"
  ),
  folder = c(
    "DE76NGSUKBR131210_1_counts",
    "DE49NGSUKBR131211_2_counts",
    "DE22NGSUKBR131212_3_counts",
    "DE92NGSUKBR131213_4_counts",
    "DE65NGSUKBR131214_5_counts",
    "DE38NGSUKBR131215_6_counts",
    "DE11NGSUKBR131216_7_counts",
    "DE81NGSUKBR131217_8_counts"
  ),
  fq1 = c(
    "DE76NGSUKBR131210_1_1.fq.gz",
    "DE49NGSUKBR131211_2_1.fq.gz",
    "DE22NGSUKBR131212_3_1.fq.gz",
    "DE92NGSUKBR131213_4_1.fq.gz",
    "DE65NGSUKBR131214_5_1.fq.gz",
    "DE38NGSUKBR131215_6_1.fq.gz",
    "DE11NGSUKBR131216_7_1.fq.gz",
    "DE81NGSUKBR131217_8_1.fq.gz"
  ),
  fq2 = c(
    "DE76NGSUKBR131210_1_2.fq.gz",
    "DE49NGSUKBR131211_2_2.fq.gz",
    "DE22NGSUKBR131212_3_2.fq.gz",
    "DE92NGSUKBR131213_4_2.fq.gz",
    "DE65NGSUKBR131214_5_2.fq.gz",
    "DE38NGSUKBR131215_6_2.fq.gz",
    "DE11NGSUKBR131216_7_2.fq.gz",
    "DE81NGSUKBR131217_8_2.fq.gz"
  ),
  stringsAsFactors = FALSE
)

# ================================
# 2. Load individual samples
# ================================
seurat_list <- list()

for (i in seq_len(nrow(samples))) {
  sample_name <- paste0(samples$group[i], "_", samples$sample_id[i])
  data_dir <- file.path(base_path, samples$folder[i], "filtered_feature_bc_matrix")
  
  message("📦 Loading: ", sample_name, "  <-- ", data_dir)
  stopifnot(dir.exists(data_dir))
  
  counts <- Read10X(data.dir = data_dir)
  obj <- CreateSeuratObject(counts = counts, min.features = 100, project = sample_name)
  
  obj$group <- samples$group[i]
  obj$sample_id <- samples$sample_id[i]
  obj$timepoint <- samples$timepoint[i]
  
  seurat_list[[sample_name]] <- obj
}


# Merge all samples

merged_all <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = "Mouse_GG_Control"
)

print(table(merged_all$group))

# ================================
# 3. Join Assay5 layers (Seurat v5)
# ================================
merged_all <- JoinLayers(merged_all, assay = "RNA")
DefaultAssay(merged_all) <- "RNA"

# ================================
# 4. Gene filtering (>=10 cells)
# ================================
counts_mat <- LayerData(merged_all[["RNA"]], layer = "counts")
keep_genes <- Matrix::rowSums(counts_mat > 0) >= 10
merged_all <- subset(merged_all, features = rownames(counts_mat)[keep_genes])
cat("Genes kept: ", nrow(merged_all), "\n")

library(Seurat)

# Calculate mitochondrial percentage (mouse uses "mt-" prefix)
merged_all[["percent.mt"]] <- PercentageFeatureSet(merged_all, pattern = "^mt-")

# Ribosomal genes
merged_all[["percent.ribo"]] <- PercentageFeatureSet(merged_all, pattern = "^Rpl|^Rps")

# Hemoglobin genes (optional for brain tissue but harmless)
merged_all[["percent.hb"]] <- PercentageFeatureSet(merged_all, pattern = "^Hb[ab]")

VlnPlot(
  merged_all,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  group.by = "group",
  pt.size = 0,
  ncol = 4
) 


# ================================
# 5. Split by group for integration
# ================================
merged_all$cells <- colnames(merged_all)
sample_list <- SplitObject(merged_all, split.by = "group")
sample_list <- Filter(function(x) ncol(x) > 0, sample_list)

# ================================
# 6. Normalize & find HVGs per sample
# ================================
sample_list <- lapply(sample_list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  x <- ScaleData(x, features = VariableFeatures(x))
  x <- RunPCA(x, features = VariableFeatures(x))
  return(x)
})

# ================================
# 7. Integration (RPCA)
# ================================
features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = sample_list, anchor.features = features, reduction = "rpca")
integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(integrated) <- "integrated"

# ================================
# 8. Dimensionality reduction & clustering
# ================================

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 50)
ElbowPlot(integrated, ndims = 50)

integrated <- RunUMAP(integrated, dims = 1:15, n.neighbors = 20, min.dist = 0.2)
integrated <- FindNeighbors(integrated, dims = 1:15)
integrated <- FindClusters(integrated, resolution = 0.5, algorithm = 4)

DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)


# ================================
# 9. Sanity plots
# ================================

Idents(integrated) <- "seurat_clusters"

# ------------------------------
# Clusters
# ------------------------------
DimPlot(
  integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 2
) +
  ggtitle("Clusters") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# ------------------------------
# Groups (all samples)
# ------------------------------
DimPlot(
  integrated,
  reduction = "umap",
  group.by = "group",
  label = FALSE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 2
) +
  ggtitle("UMAP: All Groups") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ==============================================================
#  Tumor vs Control by Timepoint (P5 / P30)
# ==============================================================

# Subset integrated object for P5 and P30
p5_subset  <- subset(integrated, subset = timepoint == "P5")
p30_subset <- subset(integrated, subset = timepoint == "P30")

# ------------------------------
# P5: Tumor (Mouse_GG_P5) vs Control (Mouse_Control_P5)
# ------------------------------
DimPlot(
  p5_subset,
  reduction = "umap",
  group.by = "group",
  label = FALSE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 2,
  cols = c("#F94144", "#277DA1")  # red-blue palette
) +
  ggtitle("UMAP: Tumor vs Control (P5)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ------------------------------
# P30: Tumor (Mouse_GG_P30) vs Control (Mouse_Control_P30)
# ------------------------------
DimPlot(
  p30_subset,
  reduction = "umap",
  group.by = "group",
  label = FALSE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 2,
  cols = c("#F9844A", "#43AA8B")  # orange-green palette
) +
  ggtitle("UMAP: Tumor vs Control (P30)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# ------------------------- #
# 10. Cell Type Marker Visualization (fixed for Seurat v5)
# ------------------------- #

DefaultAssay(integrated) <- "RNA"

marker_list <- list(
  Neurons = c("Rbfox3", "Stmn2", "Snap25", "Elavl4", "Tubb3"),                                   # pan-neuronal
  Excitatory_Neurons = c("Slc17a7", "Satb2", "Tbr1", "Camk2a", "Grin1", "Neurod6"),              # glutamatergic
  Inhibitory_Neurons = c("Gad1", "Gad2", "Sst", "Pvalb", "Vip", "Npy"),                          # GABAergic
  Astrocytes = c("Aqp4", "Gfap", "Slc1a3", "Aldh1l1", "S100b", "Fgfr3"),                         # astrocytic
  Oligodendrocytes = c("Mbp", "Plp1", "Mog", "Mag", "Mobp", "Sox10"),                            # myelinating oligodendrocytes
  OPCs = c("Pdgfra", "Cspg4", "Olig2", "Olig1", "Gpr17"),                                 # oligodendrocyte progenitors
  Microglia = c("P2ry12", "Cx3cr1", "Tmem119", "Csf1r", "Aif1", "Trem2"),                        # microglial core
  Endothelial = c("Cldn5", "Pecam1", "Flt1", "Kdr", "Vwf", "Esam")                             # vascular endothelial
)


out_dir <- "path_to_folder_analysis"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)



# -------------------------
#  11. Violin & Feature Plots by Cluster (Regions Combined)
# -------------------------

# Subfolder for violin plots
violin_dir <- file.path(out_dir, "plots", "violin_plots")
dir.create(violin_dir, recursive = TRUE, showWarnings = FALSE)

for (cell_type in names(marker_list)) {
  for (gene in marker_list[[cell_type]]) {
    
    if (!(gene %in% rownames(integrated))) {
      message("Skipping missing gene: ", gene)
      next
    }
    
    cat("Plotting (by cluster, regions combined):", cell_type, "-", gene, "\n")
    
    # ---- Violin Plot (by Cluster, with dots) ----
    p1 <- tryCatch({
      VlnPlot(
        integrated,
        features = gene,
        group.by = "seurat_clusters",   # group by clusters
        pt.size = 1,                    # dot size for single-cell expression
        combine = TRUE
      ) +
        scale_y_continuous(trans = "log1p") +
        ggtitle(paste0(cell_type, ": ", gene, " (by Cluster — Regions Combined)")) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
        )
    }, error = function(e) NULL)
    if (!is.null(p1)) print(p1)
    
    # ---- Feature Plot (by Cluster on UMAP) ----
    p2 <- tryCatch({
      fplot <- FeaturePlot(
        integrated,
        features = gene,
        reduction = "umap",
        pt.size = 2,           # increased point size
        raster = TRUE
      ) +
        ggtitle(paste0(cell_type, ": ", gene, " (Feature Plot by Cluster — Regions Combined)")) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      fplot
    }, error = function(e) NULL)
    if (!is.null(p2)) print(p2)
    
    # ---- Save outputs ----
    if (!is.null(p1))
      ggsave(
        file.path(violin_dir, paste0(gene, "_VlnPlot_byCluster_combined.pdf")),
        plot = p1, width = 7, height = 5
      )
    
    
    if (!is.null(p2))
      ggsave(
        file.path(violin_dir, paste0(gene, "_FeaturePlot_byCluster_combined.pdf")),
        plot = p2, width = 7, height = 4
      )
  }
}



# Replace "GeneName" with the gene you want to check (e.g. "Slc17a7")
gene_of_interest <- "GeneName"

# Check if the gene exists in your dataset
if (gene_of_interest %in% rownames(integrated)) {
  VlnPlot(
    integrated,
    features = gene_of_interest,
    group.by = "seurat_clusters",
    pt.size = 0.5,              # dot size for individual nuclei
    cols = "lightblue"          # optional fill color
  ) +
    scale_y_continuous(trans = "log1p") +
    ggtitle(paste0("Expression of ", gene_of_interest, " across clusters")) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
} else {
  message("⚠️ Gene not found in this dataset: ", gene_of_interest)
}



# ---- Quick canonical check ----
#Show DotPlot of cannonical markers
DotPlot(integrated, features = marker_list, group.by = "seurat_clusters") +
  RotatedAxis() +
  ggtitle("Canonical Marker Expression per Seurat Cluster") +
  theme_minimal(base_size = 9) +  # smaller overall font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )



# ================================
# 12. Marker gene identification
# ================================
DefaultAssay(integrated) <- "RNA"
if (!is.null(Layers(integrated[["RNA"]]))) {
  cat("🔧 Joining RNA layers to ensure full data access...\n")
  integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])
}

cat("🔍 Finding markers for all clusters...\n")
all_markers <- FindAllMarkers(
  integrated,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

write.csv(
  all_markers,
  file = "path_to_folder_analysis/all_markers/all_markers.csv",
  row.names = FALSE
)



markers <- read.csv(
  "path_to_folder_analysis/all_markers/all_markers.csv",
  stringsAsFactors = FALSE
)


## Another way to annotate cell types based on the cluster_markers_pre_annotation

brain_markers <- list(
  Excitatory_Neuron  = c("SLC17A7", "SATB2", "CAMK2A", "TBR1"),
  Inhibitory_Neuron  = c("GAD1", "GAD2", "SLC6A1", "PVALB", "SST", "VIP"),
  Astrocyte          = c("AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GFAP"),
  Oligodendrocyte    = c("PLP1", "MBP", "MOG", "MOBP"),
  OPC                = c("PDGFRA", "CSPG4", "OLIG1", "OLIG2", "VCAN"),
  Microglia          = c("P2RY12", "CX3CR1", "CSF1R", "AIF1", "C1QA"),
  Endothelial        = c("CLDN5", "VWF", "PECAM1", "FLT1"),
  Pericyte_VLMC      = c("PDGFRB", "RGS5", "KCNJ8", "ABCC9", "DCN", "COL1A1"),
  T_Cell             = c("CD3D", "CD3E", "TRAC", "CD4", "CD8A")
)

subset(markers, cluster == 3)[1:20, c("gene", "avg_log2FC")]

for (cl in sort(unique(markers$cluster))) {
  cat("\n==============================\n")
  cat("CLUSTER", cl, "\n")
  cat("==============================\n")
  
  print(
    subset(markers, cluster == cl)[1:20, c("gene", "avg_log2FC")]
  )
}



# ================================
# 13. Annotation example
# ================================
new.cluster.ids <- c(
  "1"  = "Ex. Neuron",
  "2"  = "Ex. Neuron",
  "3"  = "Ex. Neuron",
  "4"  = "Astrocyte",
  "5"  = "Ex. Neuron",
  "6"  = "Ex. Neuron",
  "7"  = "In. Neuron",
  "8"  = "In. Neuron",
  "9"  = "Ex. Neuron",
  "10" = "In. Neuron", 
  "11"  = "OPC",
  "12"  = "Endothelial",
  "13"  = "Astrocyte",
  "14"  = "Oligodendrocyte",
  "15"  = "Microglia"
)

integrated <- RenameIdents(integrated, new.cluster.ids)
integrated$celltype <- Idents(integrated)

umap_global <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 2,
  label.size = 6
) +
  scale_color_manual(values = celltype_colors, drop = FALSE) +
  ggtitle("UMAP: Annotated Cell Types") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),   # removes inner grid
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20
    ),
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

umap_global



library(patchwork)

p5_plot <- DimPlot(
  subset(integrated, subset = timepoint == "P5"),
  reduction = "umap",
  group.by = "celltype",
  split.by = "group",
  label = FALSE,
  repel = FALSE,
  raster = TRUE,
  pt.size = 2
) +
  scale_color_manual(values = celltype_colors, drop = FALSE) +
  ggtitle("P5 — Cell Types (Tumor vs Control)") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),   # removes inner grid
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p5_plot


p30_plot <- DimPlot(
  subset(integrated, subset = timepoint == "P30"),
  reduction = "umap",
  group.by = "celltype",
  split.by = "group",
  label = FALSE,
  repel = FALSE,
  raster = TRUE,
  pt.size = 2
) +
  scale_color_manual(values = celltype_colors, drop = FALSE) +
  ggtitle("P30 — Cell Types (Tumor vs Control)") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),   # removes inner grid
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  )


p30_plot



DotPlot(integrated, features = marker_list, group.by = "celltype") +
  RotatedAxis() +
  ggtitle("Canonical Marker Expression per Annotated Cell Type") +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )




##Save RDS data
saveRDS(integrated,
        file = "path_to_folder_analysis/GG_integrated_annotated_v5.RDS")



# ============================================
# 14. Composition Analysis: % Cell Types per Group
# ============================================

# Make sure metadata columns exist
meta <- integrated@meta.data %>%
  mutate(
    group_simple = case_when(
      grepl("GG", group, ignore.case = TRUE) ~ "Tumor",
      grepl("Control", group, ignore.case = TRUE) ~ "Control",
      TRUE ~ "Other"
    ),
    combo = paste0(timepoint, "_", group_simple)
  )

# Count cells per celltype and group
cell_composition <- meta %>%
  group_by(combo, celltype) %>%
  summarise(n = n()) %>%
  mutate(percentage = n / sum(n) * 100)

# Ensure consistent ordering
cell_composition$combo <- factor(
  cell_composition$combo,
  levels = c("P5_Tumor", "P5_Control", "P30_Tumor", "P30_Control")
)


# ------------------------------------------------------------
# 🧬 Custom, High-Contrast Color Palette
# ------------------------------------------------------------
celltype_colors <- c(
  "Ex. Neuron"       = "#E41A1C",  # vivid red
  "In. Neuron"       = "#FB9A06",  # amber-orange
  "Astrocyte"        = "#4DAF4A",  # bright green
  "Oligodendrocyte"  = "#377EB8",  # blue
  "OPC"              = "#984EA3",  # purple
  "Microglia"        = "#E7298A",  # magenta-pink
  "Endothelial"      = "#BDBDBD"   # neutral grey
)

# ------------------------------------------------------------
# 🧠 Stacked Barplot
# ------------------------------------------------------------
ggplot(cell_composition, aes(x = combo, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = celltype_colors) +
  ggtitle("Cell Type Composition (%) by Condition and Timepoint") +
  xlab("") + ylab("Percentage of Cells") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank()
  )



# ============================================
# 15. Differential Expression: Tumor vs Control per Cell Type & Timepoint
# ============================================

# ====================================================
# 🧬 Fix Seurat v5: Join RNA layers for DE analysis
# ====================================================

DefaultAssay(integrated) <- "RNA"

if (!is.null(Layers(integrated[["RNA"]]))) {
  cat("🔧 Joining RNA layers into a single counts/data matrix...\n")
  integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])
}

# confirm
Layers(integrated[["RNA"]])


library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "celltype"

de_outdir <- file.path(out_dir, "DEG")
if (!dir.exists(de_outdir)) dir.create(de_outdir, recursive = TRUE)

logfc_threshold <- 0.5
fdr_threshold <- 0.05

# --------------------------------------------
# Function for per-celltype comparisons
# --------------------------------------------
run_de_analysis <- function(timepoint_sel, celltype_sel) {
  subset_obj <- subset(integrated, subset = timepoint == timepoint_sel & celltype == celltype_sel)
  
  if (length(unique(subset_obj$group_simple)) < 2) {
    message("Skipping ", celltype_sel, " at ", timepoint_sel, " (one condition missing)")
    return(NULL)
  }
  
  Idents(subset_obj) <- "group_simple"
  
  # Run DE without pre-filtering by logFC
  degs <- FindMarkers(
    subset_obj,
    ident.1 = "Tumor",
    ident.2 = "Control",
    assay = "RNA",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1
  )
  
  if (nrow(degs) == 0) return(NULL)
  degs$gene <- rownames(degs)
  
  # Save *unfiltered* results for volcano plots
  out_file_all <- file.path(de_outdir, paste0("DEG_allGenes_", celltype_sel, "_", timepoint_sel, ".csv"))
  write.csv(degs, out_file_all, row.names = FALSE)
  
  # Filter for significant DEGs (for violin plots)
  sig_degs <- degs %>% filter(p_val_adj < fdr_threshold & abs(avg_log2FC) > logfc_threshold)
  
  if (nrow(sig_degs) == 0) {
    message("No significant DEGs for ", celltype_sel, " at ", timepoint_sel)
    return(degs)
  }
  
  # Save filtered table
  out_file_sig <- file.path(de_outdir, paste0("DEG_significant_", celltype_sel, "_", timepoint_sel, ".csv"))
  write.csv(sig_degs, out_file_sig, row.names = FALSE)
  
  # Plot top genes (by abs logFC)
  top_genes <- head(sig_degs$gene[order(-abs(sig_degs$avg_log2FC))], 6)
  p_list <- list()
  
  for (gene in top_genes) {
    p <- VlnPlot(
      subset_obj,
      features = gene,
      group.by = "group_simple",
      pt.size = 0.5
    ) +
      scale_fill_manual(values = c("Tumor" = "#E41A1C", "Control" = "#377EB8")) +
      ggtitle(paste0(celltype_sel, " — ", timepoint_sel, " — ", gene)) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        legend.position = "none"
      )
    p_list[[gene]] <- p
  }
  
  combined <- wrap_plots(p_list, ncol = 3)
  ggsave(
    file.path(de_outdir, paste0("Violin_", celltype_sel, "_", timepoint_sel, ".pdf")),
    plot = combined, width = 10, height = 5
  )
  
  return(degs)
}

# ====================================================
# Add Tumor vs Control label
# ====================================================
integrated$group_simple <- case_when(
  grepl("GG", integrated$group, ignore.case = TRUE) ~ "Tumor",
  grepl("Control", integrated$group, ignore.case = TRUE) ~ "Control",
  TRUE ~ "Other"
)

table(integrated$timepoint, integrated$group_simple)

# Run DE for all cell types and timepoints
celltypes <- levels(integrated$celltype)
timepoints <- c("P5", "P30")

all_degs <- list()
for (tp in timepoints) {
  for (ct in celltypes) {
    message("🔍 Running DE for ", ct, " — ", tp)
    res <- run_de_analysis(tp, ct)
    if (!is.null(res)) {
      all_degs[[paste(ct, tp, sep = "_")]] <- res
    }
  }
}

# Save combined summary
if (length(all_degs) > 0) {
  all_degs_df <- bind_rows(all_degs, .id = "comparison")
  write.csv(all_degs_df, file.path(de_outdir, "All_DEGs_summary.csv"), row.names = FALSE)
}

# ============================================
# Volcano plots per cell type × timepoint (show all genes)
# ============================================
library(ggplot2)

volcano_outdir <- file.path(out_dir, "DEG_VolcanoPlots")
if (!dir.exists(volcano_outdir)) dir.create(volcano_outdir, recursive = TRUE)

for (name in names(all_degs)) {
  df <- all_degs[[name]]
  if (is.null(df) || nrow(df) == 0) next
  
  if (!all(c("avg_log2FC", "p_val_adj") %in% colnames(df))) next
  
  df <- df[complete.cases(df[, c("avg_log2FC", "p_val_adj")]), ]
  if (nrow(df) == 0) next
  
  df$gene <- rownames(df)
  
  # Assign color groups
  df$Regulation <- "NotSig"
  df$Regulation[df$avg_log2FC > logfc_threshold & df$p_val_adj < fdr_threshold] <- "Up"
  df$Regulation[df$avg_log2FC < -logfc_threshold & df$p_val_adj < fdr_threshold] <- "Down"
  
  # Parse cell type and timepoint
  parts <- strsplit(name, "_")[[1]]
  cell_type <- paste(parts[-length(parts)], collapse = "_")
  time_point <- parts[length(parts)]
  
  message("📊 Volcano (all genes): ", cell_type, " — ", time_point)
  
  # Volcano plot showing all genes
  # Volcano plot showing all genes
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Regulation)) +
    geom_point(alpha = 0.8, size = 2.2) +   # slightly larger points
    scale_color_manual(
      values = c("Down" = "#377EB8", "NotSig" = "grey80", "Up" = "#E41A1C")
    ) +
    geom_vline(
      xintercept = c(-logfc_threshold, logfc_threshold),
      linetype = "dashed",
      color = "black",
      linewidth = 0.6
    ) +
    geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype = "dashed",
      color = "black",
      linewidth = 0.6
    ) +
    labs(
      title = paste0(cell_type, " — ", time_point),
      x = expression("log"[2] * " Fold Change (Tumor vs Control)"),
      y = expression("-log"[10] * " Adjusted p-value")
    ) +
    
    # >>> Bigger overall font size <<<
    theme_minimal(base_size = 16) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 20                      # larger title
      ),
      axis.title = element_text(
        face = "bold",
        size = 18                      # axis titles
      ),
      axis.text = element_text(
        size = 14                      # axis number labels
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),  # legend labels larger
      legend.position = "right"
    )
  
  # Save PDF
  ggsave(
    filename = file.path(
      volcano_outdir,
      paste0("Volcano_", cell_type, "_", time_point, ".pdf")
    ),
    plot = p,
    width = 7,
    height = 6
  )
}


deg_counts <- lapply(names(all_degs), function(name) {
  df <- all_degs[[name]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  # Extract cell type + timepoint
  parts <- strsplit(name, "_")[[1]]
  timepoint <- parts[length(parts)]
  cell_type <- paste(parts[-length(parts)], collapse = "_")
  
  # Count significant DEGs
  sig_count <- df %>%
    filter(p_val_adj < fdr_threshold & abs(avg_log2FC) > logfc_threshold) %>%
    nrow()
  
  data.frame(
    celltype = cell_type,
    timepoint = timepoint,
    sig_deg = sig_count
  )
})

deg_counts <- bind_rows(deg_counts)
deg_counts


plot_deg_pie <- function(deg_counts, timepoint_sel) {
  
  df <- deg_counts %>%
    filter(timepoint == timepoint_sel)
  
  # If some cell types have 0 DEGs, include them anyway
  # so they appear in the legend.
  df <- df %>%
    complete(celltype = names(celltype_colors), fill = list(sig_deg = 0))
  
  if (nrow(df) == 0) {
    stop(paste("No DEG data available for", timepoint_sel))
  }
  
  p <- ggplot(df, aes(x = "", y = sig_deg, fill = celltype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    scale_fill_manual(
      values = celltype_colors,
      drop = FALSE           # <<< CRITICAL: prevents ggplot from dropping legend levels
    ) +
    ggtitle(paste("Significant DEGs per Cell Type —", timepoint_sel)) +
    theme_void(base_size = 13) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_blank(),
      legend.text  = element_text(size = 10)
    )
  
  return(p)
}



# -------------------------------
# Generate both pie charts of DEG
# -------------------------------
pie_P5  <- plot_deg_pie(deg_counts, "P5")
pie_P30 <- plot_deg_pie(deg_counts, "P30")

pie_P5
pie_P30

# Optional: Save to PDF
ggsave(file.path(out_dir, "DEG_PieChart_P5.pdf"),  pie_P5,  width = 6, height = 6)
ggsave(file.path(out_dir, "DEG_PieChart_P30.pdf"), pie_P30, width = 6, height = 6)



library(tidyverse)

# ================================
# 16. Load + Plot GO Terms
# ================================

library(tidyverse)
library(janitor)

go_folder <- "path_to_folder_analysis/DEG_analysis/GO_BP"
go_files <- list.files(go_folder, pattern = "\\.csv$", full.names = TRUE)

go_outdir <- file.path(go_folder, "GO_plots")
if (!dir.exists(go_outdir)) dir.create(go_outdir)


# =======================================================
# Function to load + plot GO terms (row 2 = real header)
# =======================================================
plot_go_terms <- function(file, outfile, title_text) {
  
  # Read entire CSV
  raw <- read.csv(file, header = FALSE, check.names = FALSE)
  
  # Row 2 = real column names
  colnames(raw) <- raw[2, ]
  
  # Data starts from row 3
  df <- raw[-c(1, 2), ] %>% clean_names()
  
  # Rename columns explicitly based on Panther fixed structure
  df <- df %>% rename(
    term              = 1,
    list_count        = 2,
    go_count          = 3,
    expected          = 4,
    fold_enrichment   = 5,
    direction         = 6,
    raw_p             = 7,
    fdr               = 8
  )
  
  # Convert FDR to numeric
  df$fdr <- as.numeric(df$fdr)
  
  # Clean Fold Enrichment values (handles >100, > 250, etc.)
  df$fold_enrichment <- df$fold_enrichment %>%
    str_replace_all("[^0-9\\.]", "") %>%   # keep only digits + decimal
    as.numeric()
  
  # Replace NA fold enrichment with max FE
  max_fe <- max(df$fold_enrichment, na.rm = TRUE)
  df$fold_enrichment[is.na(df$fold_enrichment)] <- max_fe
  
  # ======================================
  # Sort by Fold Enrichment (descending)
  # ======================================
  df <- df %>%
    arrange(desc(fold_enrichment)) %>%   # highest FE first
    slice_head(n = 25)                   # top 25 GO terms
  
  # ======================================
  # Horizontal bar plot
  # ======================================
  p <- ggplot(df, aes(x = fold_enrichment,
                      y = reorder(term, fold_enrichment))) +   # MAX FE ON TOP
    
    geom_col(fill = "grey60") +
    
    labs(
      title = title_text,
      x = "Fold Enrichment",
      y = ""
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold")
    )
  
  # Save the PDF
  ggsave(outfile, p, width = 10, height = 8)
  
  return(p)
}



# =======================================================
# LOOP THROUGH ALL GO FILES
# =======================================================
for (file in go_files) {
  
  fname <- basename(file) %>% gsub(".csv", "", .)
  title_text <- paste0(gsub("_", " ", fname), " — GO Biological Processes")
  outfile <- file.path(go_outdir, paste0(fname, "_GOplot.pdf"))
  
  message("📊 Plotting: ", title_text)
  
  plot_go_terms(file, outfile, title_text)
}

