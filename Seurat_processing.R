MG <- readRDS('D://ALL_MICROGLIA_OLGA_MILICA.rds')
head(MG)
table(MG$condition)
DimPlot(MG, group.by = 'condition')
FeaturePlot(MG, features = c('Apoe'))
MG_Olga <- subset(MG, subset = condition == 'Host upon transplantation', invert = T)
DefaultAssay(MG_Olga)

#found one batch is outlying - this is ONC 5d
table(MG_Olga$Timepoint)
MG_Olga <- subset(MG_Olga, subset = Timepoint == '5d_afterONC', invert = T)
ProcessInt <- function(data.integrated){
data.integrated <- ScaleData(data.integrated, verbose = T) #, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score")
data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
data.integrated <- FindClusters(data.integrated, resolution = 1)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}
MG_Olga <- ProcessInt(MG_Olga)

DefaultAssay(MG_Olga) <- 'RNA'
DimPlot(MG_Olga, group.by = 'condition')
DimPlot(MG_Olga, group.by = 'condition', split.by = 'condition', ncol = 2) + NoLegend()
FeaturePlot(MG_Olga, features = c('Tmem119','P2ry12','Apoe','Spp1'), ncol = 2)


MG_Olga <- SetIdent(MG_Olga, value = 'condition')
markersr <- FindAllMarkers(MG_Olga, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
test <- markersr %>%
group_by(cluster) %>%
  slice_max(n=25, order_by = avg_log2FC)
markersr %>%
group_by(cluster) %>%
top_n(n=20, wt = avg_log2FC) -> top20
DefaultAssay(MG_Olga) <- 'RNA'
DoHeatmap(MG_Olga, features = top20$gene) + NoLegend()
test

DotPlot(MG_Olga, features = rev(c('Iba1','Cartpt', 'Nnat', 'Mdk', 'Mest', 'mt-Atp6',
'mt-Nd2', 'P2ry12', 'Hspa1b', 'H2-Ab1', 'Igfbpl1','Apoe','Timd4','Stab1','Stab2','Olr1','Jmjd6','Ager')), cols = c('blue','red'))

DotPlot(MG_Olga, features = rev(c('Fasl')), cols = c('blue','red'))
