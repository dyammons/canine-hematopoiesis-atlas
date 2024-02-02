#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin additional preprocessing   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in preprocessed data
seu.obj <- readRDS("../output/s3/allCellslabTransfer_S3.rds")
outName <- "allCells"
outName_full <- "240201_bm_cd34_clusterID_integrated.harmony"
clusID_main <- "clusterID_integrated.harmony"


### Umap by minorIdent before removal of disconnected cell types
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


### Umap by clusterID before removal of disconnected cell types
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated.harmony",
              group.by = clusID_main,
              pt.size = 0.1,
              label = T,
              ncol = 2,
              label.box = T,
              repel = F,
              split.by = "cellSource"
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_umap_harmony_all_clusID.png"), width = 14, height = 7)


### Supp data - generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = clusID_main, numOfFeats = 24, outName = outName,
            outDir = paste0("../output/viln/", outName, "/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"),
            min.pct = 0.25, only.pos = T)


### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "../output/cb_input/", 
                markers = "../output/viln/allCells/240201_bm_cd34_clusterID_integrated.harmony_clusterID_integrated.harmony_gene_list.csv",
                reduction = "umap.integrated.harmony",  
                colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                "majorID", clusID_main, "name", "cellSource"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                          "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                          "CD4", "MS4A1", "PPBP", "HBM")
)    

#remove disconnected populations that do not connect to lineages (T cells = 14 & 19, Macrophage = 18, endothelial = 20, and Plasma cells = 17)
seu.obj <- subset(seu.obj, invert = T, subset = clusterID_integrated.harmony %in% c(14,19,18,20,17))

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(seu.obj = seu.obj, dout = "../output/s2/", outName = outName_full, runAllMethods = TRUE,
                        indReClus = T)

# #use clustree to identify clustering parameters that appear most appropriate
# clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = experiment, 
#             test_dims = c(45,40,35), algorithm = 3, prefix = "RNA_snn_res.")

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}

saveRDS(seu.obj, "../output/s3/240201_bm_cd34_removed_disconnected.rds")

################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end additional preprocessing   ######## <<<<<<<<<<<<<<
################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Begin lineage analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# #load in the object
seu.obj <- readRDS("../output/s3/240201_bm_cd34_removed_disconnected.rds")
outName <- "allCells_clean"
outName_full <- "240201_bm_cd34_clusterID_integrated.harmony"
clusID_main <- "clusterID_integrated.harmony"

#load metadata
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "cellSource")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$clusterID_integrated.harmony <- factor(seu.obj$clusterID_integrated.harmony, levels = 0:max(as.numeric(names(table(seu.obj$clusterID_integrated.harmony)))))
seu.obj <- convertTOclusID(seu.obj = seu.obj, metaSlot = "majorID", newMetaName = "major_clusID")


### UMAP of the reduced data
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated.harmony", 
              group.by = clusID_main,
              pt.size = 0.1,
              label = T,
              label.box = T,
              repel = F,
)
p <- cusLabels(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_clusID_main_UMAP.png"), width = 7, height = 7)


### Supp data - generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = clusID_main, numOfFeats = 24, outName = outName,
            outDir = paste0("../output/viln/", outName, "/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"),
            min.pct = 0.25, only.pos = T)


### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "../output/cb_input/", 
                markers = paste0("../output/viln/", outName, "/", outName, "_", clusID_main, "_gene_list.csv"),
                reduction = "umap.integrated.harmony",  
                colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                "majorID", clusID_main, "name", "cellSource"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                          "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                          "CD4", "MS4A1", "PPBP", "HBM")
)    








# majorColz <- c("#DFA74D","#B7DBF0", "#C59979", "#ECD58E", "#D89183" )

# #rename idents
# Idents(seu.obj) <- "clusterID_integrated.harmony"
# clusId <- c("neut","neut","neut","neut","neut","mono","mono","eryth","eryth","HPSC","mast","bcell","bcell","bcell","bcell","bcell","bcell")
# names(clusId) <- c(0,1,4,5,6,9,8,13,15,2,10,16,14,12,7,3,11)
# seu.obj <- RenameIdents(seu.obj, clusId) 
# seu.obj$braches <- Idents(seu.obj)
# levels(seu.obj$braches) # "neut"  "mono"  "eryth" "HPSC"  "mast"  "bcell"

# plot <- DimPlot(seu.obj, 
#               reduction = "umap.integrated.harmony", 
#               group.by = "braches",
#               pt.size = 0.25,
#               label = TRUE,
#               label.box = TRUE
#  )

# p <- formatUMAP(plot = plot) + NoLegend()
# ggsave(paste("../output/", outName, "/", outName, "_bracnch.png", sep = ""), width = 7, height = 7)


# ### Generate violin plots of defining features
# vilnPlots(seu.obj = seu.obj, groupBy = "braches", numOfFeats = 24, outName = "bm_cd34_subset_braches",
#             outDir = "../output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
#             min.pct = 0.25, only.pos = T)

# #gsea on susspected mast cells to confirm ID -- does not clearly confirm, susspect they are HPSCs
# markers.df <- read.csv("../output/viln/allCells/bm_cd34_subset_braches_gene_list.csv")
# susMastCells <- markers.df[markers.df$cluster == "mast", ] %>% filter(p_val_adj < 0.01) %>% pull(gene)
# p <- plotGSEA(
#     geneList = susMastCells,
#     category = "C8", species = "dog", upCol = "red", dwnCol = "blue",
#     pvalueCutoff = 0.05, subcategory = NULL, termsTOplot = 16, upOnly = T, trimTerm = T
# )
# ggsave(paste("../output/", outName, "/", outName, "_mast_cell_gsea.png", sep = ""), width = 10, height = 7)

# #rename idents to correct mast cell call to be HPSCs
# Idents(braches) <- "clusterID_integrated.harmony"
# seu.obj <- RenameIdents(seu.obj, "mast" = "HPSC") 
# seu.obj$braches <- Idents(seu.obj)
# levels(seu.obj$braches) # "HPSC"  "neut"  "mono"  "eryth" "bcell"

# #convert to cluster ID numbers
# seu.obj <- convertTOclusID(seu.obj = seu.obj, metaSlot = "braches")


# ### Generate violin plots of defining features
# vilnPlots(seu.obj = seu.obj, groupBy = "braches_clusID", numOfFeats = 24, outName = "bm_cd34_subset_braches",
#             outDir = "../output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
#             min.pct = 0.25, only.pos = T)


# pi <- autoDot(seu.integrated.obj = seu.obj, inFile = "../output/viln/allCells/bm_cd34_subset_braches_gene_list.csv", groupBy = "braches_clusID",
#                      MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1,
#                     filterTerm = "ENSCAFG"
#                     ) + theme(legend.box="vertical") + scale_fill_manual(values = majorColz)
# ggsave(paste0("../output/", outName, "/", outName, "_autodot_branches.png"), width = 5, height = 10)


# lapply(c("umap.integrated.harmony", "umap.integrated.joint", "umap.integrated.rcpa"), function(x){
#     pi <- DimPlot(seu.obj, 
#                   reduction = x,
#                   group.by = "braches_clusID",
#                   pt.size = 0.1,
#                   label = T,
#                   label.box = T,
#                   repel = F,
#     )
#     p <- cusLabels(plot = pi, smallAxes = T) + NoLegend()
#     ggsave(plot = p, paste0("../output/", outName, "/", outName, "_",x, ".png"), width = 7, height = 7)
# })


# #modify the reduction parameters
# seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "subbed_integrated.harmony", 
#                         final.dims = 45, final.res = 1.4, stashID = "clusterID", algorithm = 3, min.dist = 0.2, n.neighbors = 20,
#                         prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.harmony",
#                         saveRDS = F, return_obj = T, returnFeats = T,
#                         features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
#                                         "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
#                                         "CD4", "MS4A1", "PPBP","HBM")
# )


# ### plot the data
# pi <- DimPlot(seu.obj, 
#               reduction = "umap.integrated.harmony",
#               group.by = "braches_clusID",
#               pt.size = 0.1,
#               cols = majorColz,
#               label = T,
#               label.box = T,
#               repel = F,
# )
# p <- cusLabels(plot = pi, smallAxes = T) + NoLegend()
# ggsave(plot = p, paste0("../output/", outName, "/", outName, "_braches_clusID.png"), width = 7, height = 7)



# ### Fig supp: run SlingShot
# sce.obj <- as.SingleCellExperiment(seu.obj)
# rd1 <- Embeddings(seu.obj, reduction = "pca")[,1:2]
# rd2 <- Embeddings(seu.obj, reduction = "umap.integrated.harmony")
# colnames(rd2) <- c('UMAP1', 'UMAP2')

# #assign origin -- left most hpsc population
# start.clus <- '2'
# end.clus <- c('1','8','15','11','16','7','10')

# reducedDims(sce.obj) <- SimpleList(PCA = rd1, UMAP = rd2)

# sce.obj <- slingshot(sce.obj, clusterLabels = 'clusterID_integrated.harmony', reducedDim = 'UMAP', start.clus = start.clus, end.clus = end.clus)

# #identify lineages
# lin1 <- getLineages(Embeddings(seu.obj, reduction = "umap.integrated.harmony"), sce.obj$clusterID_integrated.harmony, start.clus = start.clus, end.clus = end.clus)
# branchData <- SlingshotDataSet(lin1)@lineages

# #plot the lineages
# plot <- DimPlot(seu.obj, 
#               reduction = "umap.integrated.harmony", 
#               group.by = "clusterID_integrated.harmony",
#               pt.size = 0.25,
#               label = TRUE,
#               label.box = TRUE
#  )

# p <- cleanSling(plot = plot, shape = 21, labCol = "black", size = 8, alpha = 1, rm.na = T, branchData = branchData) + NoLegend()
# ggsave(paste("../output/", outName, "/", outName, "_bracnch.png", sep = ""), width = 7, height = 7)

# #get the meta data over - to the seurat object
# metaData <- as.data.frame(sce.obj@colData@rownames)
# colnames(metaData) <- "barcode"
# metaData$slingPseudotime_1 <- sce.obj$slingPseudotime_1
# metaData$slingPseudotime_2 <- sce.obj$slingPseudotime_2
# metaData$slingPseudotime_3 <- sce.obj$slingPseudotime_3
# metaData$slingPseudotime_4 <- sce.obj$slingPseudotime_4
# metaData$slingPseudotime_5 <- sce.obj$slingPseudotime_5
# metaData$slingPseudotime_6 <- sce.obj$slingPseudotime_6
# # metaData$slingPseudotime_7 <- sce.obj$slingPseudotime_7

# metaData <- metaData %>% rowwise %>% mutate(pseudoTime = mean(c(slingPseudotime_1,slingPseudotime_2,slingPseudotime_3,slingPseudotime_4,slingPseudotime_5,slingPseudotime_6),na.rm=TRUE)) #,slingPseudotime_7

# seuMeta <- seu.obj@meta.data %>% mutate(barcode = rownames(.))
# seuMeta <- seuMeta[ ,!grepl("slingPseudotime|pseudoTime", colnames(seuMeta))]

# newMeta <- seuMeta %>% left_join(metaData, by = 'barcode')
# rownames(newMeta) <- newMeta$barcode
# seu.obj@meta.data <- newMeta

# #plot the lineages on original umap rep
# features <- c("slingPseudotime_1","slingPseudotime_2","slingPseudotime_3","slingPseudotime_4","slingPseudotime_5","slingPseudotime_6")
# titles <- c("Lineage 1","Lineage 2","Lineage 3","Lineage 4","Lineage 5","Lineage 6")

# p <- prettyFeats(seu.obj = seu.obj, reduction = "umap.integrated.harmony",nrow = 2, ncol = 3, features = features, color = "black", order = F, titles = titles, noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

# ggsave(paste("../output/", outName, "/", outName, "_pseudoTime.png", sep = ""), width = 9, height = 6)


# features <- "pseudoTime"
# titles <- "All time"

# p <- prettyFeats(seu.obj = seu.obj, reduction = "umap.integrated.harmony",nrow = 1, ncol = 1, features = features, color = "black", order = F, titles = titles, noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

# ggsave(paste("../output/", outName, "/", outName, "_pseudoTime_all.png", sep = ""), width = 3, height = 3)























