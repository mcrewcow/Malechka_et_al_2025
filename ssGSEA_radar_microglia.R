

library(escape)
DefaultAssay(MG_Olga) <- 'RNA'
gene.sets1 <- getGeneSets(library = "C5", gene.sets = c('GOBP_MICROGLIAL_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                                                        'GOBP_MICROGLIAL_CELL_MEDIATED_CYTOTOXICITY',
                                                        'GOBP_MICROGLIAL_CELL_MIGRATION',
                                                        'GOBP_MICROGLIAL_CELL_PROLIFERATION',
                                                        'GOBP_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION',
                                                        'GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE',
                                                        'GOBP_ACUTE_INFLAMMATORY_RESPONSE',
                                                        'GOBP_NEUROINFLAMMATORY_RESPONSE',
                                                        'GOBP_PHAGOCYTOSIS',
                                                        'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
                                                        'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II',
                                                        'GOBP_MICROGLIA_DIFFERENTIATION',
                                                        'GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION',
                                                        'GOBP_MYELOID_CELL_DEVELOPMENT',
                                                        'GOBP_MYELOID_PROGENITOR_CELL_DIFFERENTIATION'
),species = 'Mus musculus')
ES <- enrichIt(obj = MG_Olga,
               gene.sets = gene.sets1,
               groups = 1000)
MG_Olga<- AddMetaData(MG_Olga, ES)
ES2 <- data.frame(MG_Olga[[]], Idents(MG_Olga))
colnames(ES2)[ncol(ES2)] <- "cluster"


library(dplyr)
library(fmsb)

# Calculate the average of each pathway grouped by condition
avg_pathway <- ES2 %>%
  group_by(Timepoint) %>%
  summarise(across(starts_with("GOBP"), mean, .names = "avg_{.col}"))

# Min-max normalization for each pathway
normalized_pathways <- avg_pathway %>%
  mutate(across(starts_with("avg_GOBP"),
                ~ (. - min(.)) / (max(.) - min(.)),
                .names = "scaled_{.col}"))

# Prepare data for radar chart
radar_data <- normalized_pathways %>%
  pivot_longer(-Timepoint, names_to = "Pathway", values_to = "Average") %>%
  pivot_wider(names_from = Timepoint, values_from = Average)

# Add min and max rows (required for fmsb radar chart)
radar_data <- rbind(
  max = rep(1, ncol(radar_data) - 1),  # Max normalized value
  min = rep(0, ncol(radar_data) - 1),  # Min normalized value
  radar_data
)
radarchart(
  radar_data[-1],  # Remove Pathway column for plotting
  axistype = 1,  # Semi-transparent fill
  plwd = 2,  # Line width for polygons
  cglcol = "grey", cglty = 1, cglwd = 0.8,  # Grid style
  axislabcol = "black",  # Axis label color
  vlcex = 0.8  # Vertex label size
)
title("Polygon Chart of Pathways by Condition")
radar_data[3,]
# Define the pathways (replace with actual pathway names if available)
# Shorten or wrap pathway names
radar_data <- radar_data[-14,]#twice
radar_data <- radar_data[-10,]
radar_data <- radar_data[-25,]
radar_data <- radar_data[-21,]
pathways <- strwrap(radar_data$Pathway[3:13], width = 20)  # Wrap text to 20 characters

# Add the legend for pathways
legend(
  "bottomright",  # Adjust position as needed
  legend = pathways,
  col = "black",  # All axes are black in radar chart
  lty = 1,        # Line type
  lwd = 2,        # Line width
  bty = "n",      # No border for legend
  title = "Pathways",
  cex = 0.7       # Adjust text size for better fit
)

pathways <- radar_data$Pathway[3:13]
pathway_colors <- rainbow(length(pathways))  # Generate distinct colors for each pathway

# Plot the radar chart
radarchart(
  radar_data[-1],  # Remove Pathway column for plotting
  axistype = 1,
  pcol = pathway_colors,  # Assign distinct colors for each pathway
  plty = 1,               # Solid line style
  plwd = 2,               # Line width for polygons
  cglcol = "grey",        # Grid line color
  cglty = 1, cglwd = 0.8, # Grid style
  axislabcol = "black",   # Axis label color
  vlcex = 0.8             # Vertex label size
)



# Add the legend for pathways with distinct colors
legend(
  "bottomright",  # Adjust position as needed
  legend = pathways,
  col = pathway_colors,  # Assign distinct colors to each pathway
  pch = 16,              # Matching point style
  pt.cex = 1.2,          # Matching point size
  lty = 1,               # Line type
  lwd = 2,               # Line width
  bty = "n",             # No border for legend
  title = "Pathways",
  cex = 0.7              # Adjust text size for better fit
)


normalized_pathways$Timepoint <- factor(
  normalized_pathways$Timepoint,
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

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_ACUTE_INFLAMMATORY_RESPONSE, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()


ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_MICROGLIAL_CELL_MEDIATED_CYTOTOXICITY, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_MICROGLIAL_CELL_MIGRATION, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_MICROGLIAL_CELL_PROLIFERATION, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_MICROGLIA_DIFFERENTIATION, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_MYELOID_CELL_DEVELOPMENT, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_MYELOID_PROGENITOR_CELL_DIFFERENTIATION, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_NEUROINFLAMMATORY_RESPONSE, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_PHAGOCYTOSIS, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

ggplot(normalized_pathways, aes(x = Timepoint, y = normalized_pathways$avg_GOBP_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION, group = 1)) +
  geom_line(color = "blue") +         # Connect the dots with a line
  geom_point(size = 3, color = "red") + theme_bw()

df_for_heatmap <- normalized_pathways %>%
  select(Timepoint, starts_with("scaled_"))

# Convert the scaled columns to a numeric matrix
data_mat <- as.matrix(df_for_heatmap %>% select(-Timepoint))

# Assign row names to be the Timepoint
rownames(data_mat) <- df_for_heatmap$Timepoint

data_mat_t <- t(data_mat)

pheatmap(
  data_mat_t,
  cluster_rows = TRUE,  # set to TRUE if you want to cluster timepoints
  cluster_cols = FALSE,  # set to TRUE if you want to cluster pathways
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of Scaled Pathway Scores by Timepoint"
)
