#!/usr/bin/Rscript

#load custom functions & packages
source("./customFunctions_Seuratv5.R")


#load in preprocessed data
seu.obj <- readRDS("../output/s3/bm_cd34_analysis_231211_v5_integrated_res0.6_dims45_dist0.1_neigh10_S3.rds")
outName <- "allCells"


### Fig extra - umap by minorIdent before removal of disconnected cell types
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated.harmony",
              group.by = "minorIdent",
              pt.size = 0.1,
              label = T,
              ncol = 2,
              label.box = T,
              repel = F,
              split.by = "cellSource"
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_umap_harmony_all_minorIdent.png"), width = 14, height = 7)


### Fig extra - umap by clusterID before removal of disconnected cell types
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated.harmony",
              group.by = "clusterID_integrated.harmony",
              pt.size = 0.1,
              label = T,
              ncol = 2,
              label.box = T,
              repel = F,
              split.by = "cellSource"
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_umap_harmony_all_clusID.png"), width = 14, height = 7)


### Fig extra - umap highlighting c21 before removal of disconnected cell types
Idents(seu.obj) <- "clusterID_integrated.harmony"
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated.harmony",
              group.by = "clusterID_integrated.harmony",
              pt.size = 0.1,
              label = T,
              ncol = 2,
              cells.highlight = WhichCells(seu.obj, ident = "21"),
              label.box = F,
              repel = F,
              split.by = "cellSource"
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_umap_harmony_all_clusID.png"), width = 14, height = 7)


### Supp data - generate violin plots of defining features
# vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_integrated.harmony", numOfFeats = 24, outName = "bm_cd34_clusterID_integrated.harmony",
#             outDir = "../output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
#             min.pct = 0.25, only.pos = T)


#remove disconnected populations that do not connect to lineages (T cells, Macrophage, and Plasma cells)
seu.obj <- subset(seu.obj, invert = T, subset = clusterID_integrated.harmony %in% c(16,19,20,21,22,23))



#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(dout = "../output/s2/", outName = "bm_cd34_subset_analysis_231211_v5", runAllMethods = TRUE,
                        indReClus = T, seu.obj = seu.obj)


# #use clustree to identify clustering parameters that appear most appropriate
# clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = experiment, 
#             test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- readRDS("../output/s2/bm_cd34_subset_analysis_231211_v5_integrated_S2.rds")
#complete data visualization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "integrated.cca", 
                        final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                        prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.cca",
                        saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)

#complete data visualization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "integrated.harmony", 
                        final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                        prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.harmony",
                        saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)

#complete data visualization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "integrated.joint", 
                        final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                        prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.joint",
                        saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)

#complete data visualization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "integrated.rcpa", 
                        final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                        prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.rcpa",
                        saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)

# saveRDS(seu.obj, "../output/s3/bm_cd34_subset_analysis_231211_v5_integrated_res0.6_dims45_dist0.1_neigh10_S3.rds")
seu.obj <- readRDS("../output/s3/bm_cd34_subset_analysis_231211_v5_integrated_res0.6_dims45_dist0.1_neigh10_S3.rds")

Idents(seu.obj) <- "clusterID_integrated.harmony"

clusId <- c("neut", "neut","neut","neut","neut","mono","mono","eryth","eryth","HPSC","mast","tcell","bcell","bcell","bcell","bcell","bcell")
names(clusId) <- c(0,1,4,5,6,9,8,13,15,2,10,16,14,12,7,3,11)

seu.obj <- RenameIdents(seu.obj, clusId) 
seu.obj$braches <- Idents(seu.obj)

# ### Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "braches", numOfFeats = 24, outName = "bm_cd34_subset_braches",
            outDir = "../output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
            min.pct = 0.25, only.pos = T)


pi <- autoDot(seu.integrated.obj = seu.obj, inFile = "../output/viln/allCells/bm_cd34_subset_braches_gene_list.csv", groupBy = "braches",
                     MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1,
                    filterTerm = "ENSCAFG"
                    )

ggsave(paste0("../output/", outName, "/", outName, "_autodot_branches.png"), width = 7, height = 14)



# ### Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_integrated.harmony", numOfFeats = 24, outName = "bm_cd34_subset_clusterID_integrated.harmony",
            outDir = "../output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
            min.pct = 0.25, only.pos = T)

lapply(c("umap.integrated.harmony", "umap.integrated.joint", "umap.integrated.rcpa"), function(x){
    pi <- DimPlot(seu.obj, 
                  reduction = x,
                  group.by = "minorIdent",
                  pt.size = 0.1,
                  label = T,
                  ncol = 2,
                  label.box = T,
                  repel = F,
                  split.by = "cellSource"
    )
    p <- formatUMAP(plot = pi) + NoLegend()
    ggsave(paste0("../output/", outName, "/", outName, "_",x, ".png"), width = 14, height = 7)
})


### Fig supp: run SlingShot
sce.obj <- as.SingleCellExperiment(seu.obj)
rd1 <- Embeddings(seu.obj, reduction = "pca")[,1:2]
rd2 <- Embeddings(seu.obj, reduction = "umap.integrated.harmony")
colnames(rd2) <- c('UMAP1', 'UMAP2')

#assign origin -- left most hpsc population
start.clus <- '2'
end.clus <- c('1','8','15','11','16','7','10')

reducedDims(sce.obj) <- SimpleList(PCA = rd1, UMAP = rd2)

sce.obj <- slingshot(sce.obj, clusterLabels = 'clusterID_integrated.harmony', reducedDim = 'UMAP', start.clus = start.clus, end.clus = end.clus)

#identify lineages
lin1 <- getLineages(Embeddings(seu.obj, reduction = "umap.integrated.harmony"), sce.obj$clusterID_integrated.harmony, start.clus = start.clus, end.clus = end.clus)
branchData <- SlingshotDataSet(lin1)@lineages

#plot the lineages
plot <- DimPlot(seu.obj, 
              reduction = "umap.integrated.harmony", 
              group.by = "clusterID_integrated.harmony",
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE
 )

p <- cleanSling(plot = plot, shape = 21, labCol = "black", size = 8, alpha = 1, rm.na = T, branchData = branchData) + NoLegend()
ggsave(paste("../output/", outName, "/", outName, "_bracnch.png", sep = ""), width = 7, height = 7)

#get the meta data over - to the seurat object
metaData <- as.data.frame(sce.obj@colData@rownames)
colnames(metaData) <- "barcode"
metaData$slingPseudotime_1 <- sce.obj$slingPseudotime_1
metaData$slingPseudotime_2 <- sce.obj$slingPseudotime_2
metaData$slingPseudotime_3 <- sce.obj$slingPseudotime_3
metaData$slingPseudotime_4 <- sce.obj$slingPseudotime_4
metaData$slingPseudotime_5 <- sce.obj$slingPseudotime_5
metaData$slingPseudotime_6 <- sce.obj$slingPseudotime_6
# metaData$slingPseudotime_7 <- sce.obj$slingPseudotime_7

metaData <- metaData %>% rowwise %>% mutate(pseudoTime = mean(c(slingPseudotime_1,slingPseudotime_2,slingPseudotime_3,slingPseudotime_4,slingPseudotime_5,slingPseudotime_6),na.rm=TRUE)) #,slingPseudotime_7

seuMeta <- seu.obj@meta.data %>% mutate(barcode = rownames(.))
seuMeta <- seuMeta[ ,!grepl("slingPseudotime|pseudoTime", colnames(seuMeta))]

newMeta <- seuMeta %>% left_join(metaData, by = 'barcode')
rownames(newMeta) <- newMeta$barcode
seu.obj@meta.data <- newMeta

#plot the lineages on original umap rep
features <- c("slingPseudotime_1","slingPseudotime_2","slingPseudotime_3","slingPseudotime_4","slingPseudotime_5","slingPseudotime_6")
titles <- c("Lineage 1","Lineage 2","Lineage 3","Lineage 4","Lineage 5","Lineage 6")

p <- prettyFeats(seu.obj = seu.obj, reduction = "umap.integrated.harmony",nrow = 2, ncol = 3, features = features, color = "black", order = F, titles = titles, noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

ggsave(paste("../output/", outName, "/", outName, "_pseudoTime.png", sep = ""), width = 9, height = 6)


features <- "pseudoTime"
titles <- "All time"

p <- prettyFeats(seu.obj = seu.obj, reduction = "umap.integrated.harmony",nrow = 1, ncol = 1, features = features, color = "black", order = F, titles = titles, noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

ggsave(paste("../output/", outName, "/", outName, "_pseudoTime_all.png", sep = ""), width = 3, height = 3)




### play around with the parameters
#complete data visualization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "subbed_integrated.harmony", 
                        final.dims = 45, final.res = 1.4, stashID = "clusterID", algorithm = 3, min.dist = 0.2, n.neighbors = 20,
                        prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.harmony",
                        saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)



















