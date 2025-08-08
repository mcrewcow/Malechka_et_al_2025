
Olga_rgc <- readRDS('D://all_rgc_Olga.rds')
table(Olga_rgc$Timepoint)
Olga_rgc$Timepoint <- Olga_rgc$Condition #transfer ONC label
if ("Timepoint" %in% colnames(Olga_rgc@meta.data) & "condition" %in% colnames(Olga_rgc@meta.data)) {
  # Replace NA or empty values in Timepoint with corresponding values from condition
  Olga_rgc@meta.data$Timepoint[is.na(Olga_rgc@meta.data$Timepoint) | Olga_rgc@meta.data$Timepoint == ""] <-
    Olga_rgc@meta.data$condition[is.na(Olga_rgc@meta.data$Timepoint) | Olga_rgc@meta.data$Timepoint == ""]
} else {
  print("Either Timepoint or condition column is missing in metadata")
}
table(Olga_rgc$condition)
Olga_rgc <- subset(Olga_rgc, subset = condition == '5 days', invert = T)
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T) #, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score")
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}
Olga_rgc <- ProcessInt(Olga_rgc)

Olga_rgc <- SetIdent(Olga_rgc, value= 'condition')
Olga_rgc <- RenameIdents(Olga_rgc, '0.5 day' = 'ONC <=4 days',
                         '1 day' = 'ONC <=4 days',
                         '2 days' = 'ONC <=4 days',
                         '4 days' = 'ONC <=4 days',
                         '7 days' = 'ONC >4 days',
                         '14 days' = 'ONC >4 days')
Olga_rgc$Timepoint <- Olga_rgc@active.ident
DimPlot(Olga_rgc, group.by = 'Timepoint', raster = F, split.by = 'Timepoint', ncol = 2) + NoLegend()
DimPlot(Olga_rgc, group.by = 'Timepoint', raster = F)

Olga_rgc$condition <- factor(
  Olga_rgc$condition,
  levels = c(
    "Development E10-P8",
    "Adult 4-52 wk",
    "0.5 day",
    "1 day",
    "2 days",
    "4 days",
    "7 days",
    "14 days"
  )
)

# 2) Set the default assay to RNA (if not already):
DefaultAssay(MG_Olga) <- "RNA"

# 1) Grab expression data for desired genes
expr_mat <- GetAssayData(Olga_rgc, slot = "data", assay = "RNA")
DotPlot(Olga_rgc, features = c("Lrp1", "Gas6", "Tulp1", "Mfge8", "Anxa1", "Calr", "C1qa", "C3", "Gla", "Pros1", "Ptdss1",

                               # Do-not-eat-me
                               "Cd47",

                               # Find-me
                               "Cx3cl1",

                               # Additional
                               "Xkr8", "Timd4", "Adgrb1", "Mertk", "Axl", "C1qb", "C1qc",
                               "Stab2", "Scarf1", "Cd300lf", "Plscr1", "Plscr3", "Ano6", "Rab7", "Rab5a", "Rab35", "Ptx3",
                               "Crp", "Cd14", "Lrp8", "Spp1", "Thbs1", "Cd93", "Sirpa"), group.by = 'Timepoint', assay = 'RNA')
my_genes <- c("Lrp1", "Gas6", "Tulp1", "Mfge8", "Anxa1", "Calr", "C1qa", "C3", "Gla", "Pros1", "Ptdss1",

              # Do-not-eat-me
              "Cd47",

              # Find-me
              "Cx3cl1",

              # Additional
              "Xkr8", "Timd4", "Adgrb1", "Mertk", "Axl", "C1qb", "C1qc",
              "Stab2", "Scarf1", "Cd300lf", "Plscr1", "Plscr3", "Ano6", "Rab7", "Rab5a", "Rab35", "Ptx3",
              "Crp", "Cd14", "Lrp8", "Spp1", "Thbs1", "Cd93", "Sirpa"
)


expr_subset <- expr_mat[my_genes, , drop = FALSE]

# 2) Identify the timepoint of each cell
cell_timepoints <- factor(Olga_rgc$condition)
tps_in_order <- levels(cell_timepoints)

# 3) Create an empty matrix to hold average expression for Genes x Timepoints
avg_mat <- matrix(
  nrow = length(my_genes),
  ncol = length(tps_in_order),
  dimnames = list(my_genes, tps_in_order)
)

# 4) Compute the mean expression for each gene within each timepoint
for (tp in tps_in_order) {
  cells_in_tp <- which(cell_timepoints == tp)
  if (length(cells_in_tp) == 0) {
    # no cells in that timepoint
    avg_mat[, tp] <- NA
  } else {
    avg_mat[, tp] <- rowMeans(expr_subset[, cells_in_tp, drop = FALSE])
  }
}

# 5) Min–max scale each gene's row across the timepoints
scaled_avg_mat <- t(apply(avg_mat, 1, function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  if (x_max == x_min) {
    # If all timepoint averages are the same, set that row's scaled values to 0
    return(rep(0, length(x)))
  } else {
    return( (x - x_min) / (x_max - x_min) )
  }
}))

# 6) Make a heatmap
my_palette <- colorRampPalette(c("blue", "white", "red"))(50)
heatmap(
  scaled_avg_mat,
  Colv = NA,         # do not cluster columns
  Rowv = TRUE,         # do not cluster rows
  scale = "none",
  col = my_palette,
  margins = c(8, 6)  # adjust margins as needed
)
