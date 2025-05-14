
MG_Olga$Timepoint <- factor(
  MG_Olga$Timepoint,
  levels = c(
    "Development E10-P8",
    "Adult 4-52 wk",
    "12h_afterONC",
    "24h_afterONC",
    "48h_afterONC",
    "4d_afterONC",
    "1w_afterONC",
    "2w_afterONC",
    "MB-induced glaucoma"
  )
)

# 2) Set the default assay to RNA (if not already):
DefaultAssay(MG_Olga) <- "RNA"

# 1) Grab expression data for desired genes
expr_mat <- GetAssayData(MG_Olga, slot = "data", assay = "RNA")

my_genes <- c('Aif1',"P2ry12","Tmem119","Ccl4","Ccl5","Spp1",
              "Igfbpl1", "Apoe",'Mki67','Cdc20','Igf1','Cd74','Fasl','Anxa1','Fpr2','Alx1',
              'Timd4','Stab1','Stab2','Olr1','Jmjd6','Ager')

my_genes <- c('Apoe','Trem2','Lpl','Bin1','Rin3','Cd2ap','Zyx','Anxa1','Anxa3','Anxa11')


expr_subset <- expr_mat[my_genes, , drop = FALSE]

# 2) Identify the timepoint of each cell
cell_timepoints <- factor(MG_Olga$Timepoint)
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
apoe_vals <- avg_mat["Apoe", ]
apoe_vals

df_apoe <- data.frame(
  Timepoint  = factor(names(apoe_vals), levels = colnames(avg_mat)),
  Expression = as.numeric(apoe_vals)
)
library(ggplot2)

ggplot(df_apoe, aes(x = Timepoint, y = Expression, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw() +ylim(2,5.5) # Plot the dots

MG_Olga$tech <- 'technical'
RidgePlot(MG_Olga, features = "Apoe", group.by = 'tech') +
  ggtitle("Ridge Plot of Apoe Expression") +
  theme_minimal()
RidgePlot(MG_Olga, features = "Apoe") +
  ggtitle("Ridge Plot of Apoe Expression") +
  theme_minimal()

MG_Olga_ApoeHigh <- subset(MG_Olga, subset = Apoe >= 4)

expr_mat <- GetAssayData(MG_Olga_ApoeHigh, slot = "data", assay = "RNA")

my_genes <- c("P2ry12","Tmem119","Ccl4","Ccl5","Spp1",
              "Igfbpl1", "Apoe",'Mki67','Cdc20','Igf1','Cd74','Fasl')
expr_subset <- expr_mat[my_genes, , drop = FALSE]

# 2) Identify the timepoint of each cell
cell_timepoints <- factor(MG_Olga_ApoeHigh$Timepoint)
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

apoe_vals <- avg_mat["Apoe", ]
apoe_vals

df_apoe <- data.frame(
  Timepoint  = factor(names(apoe_vals), levels = colnames(avg_mat)),
  Expression = as.numeric(apoe_vals)
)
library(ggplot2)

ggplot(df_apoe, aes(x = Timepoint, y = Expression, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw() +ylim(2,5.5) # Plot the dots


# 1) Store the counts in objects
tableApoeHigh <- table(MG_Olga_ApoeHigh$Timepoint)
tableAll      <- table(MG_Olga$Timepoint)

# 2) Calculate percentages
propApoeHigh <- 100 * (tableApoeHigh / tableAll)

# 3) Create a small summary table/data frame
percentage_df <- data.frame(
  Timepoint          = names(tableAll),
  n_ApoeHigh         = as.vector(tableApoeHigh),
  n_Total            = as.vector(tableAll),
  Percent_ApoeHigh   = round(as.vector(propApoeHigh), 2)
)

percentage_df
percentage_df$Timepoint <- factor(
  percentage_df$Timepoint,
  levels = c(
    "Development E10-P8",
    "Adult 4-52 wk",
    "12h_afterONC",
    "24h_afterONC",
    "48h_afterONC",
    "4d_afterONC",
    "1w_afterONC",
    "2w_afterONC",
    "MB-induced glaucoma"
  )
)
ggplot(percentage_df, aes(x = Timepoint, y = Percent_ApoeHigh, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()






