#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

##### prepare data set #####

######### MODIFY #########

#set output name -- recommend including data and sample size
experiment <- "bm_cd34_analysis_240131"
outName <- "allCells_SCT"

nFeature_RNA_high <- 5500
nFeature_RNA_low <- 100
percent.mt_high <- 12.5
nCount_RNA_high <- 75000
nCount_RNA_low <- 200

########## END MODIFY #########

# load10x(din = "../input/", dout = "../output/s1/", outName = experiment, testQC = FALSE, removeRBC_pal = FALSE,
#         nFeature_RNA_high = nFeature_RNA_high, nFeature_RNA_low = nFeature_RNA_low, percent.mt_high = percent.mt_high, 
#         nCount_RNA_high = nCount_RNA_high, nCount_RNA_low = nCount_RNA_low)

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(din = "../output/s1/", dout = "../output/s2/", outName = experiment, normalization.method = "SCT",
                         runAllMethods = TRUE, indReClus = F)

# #use clustree to identify clustering parameters that appear most appropriate
# clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = experiment, 
#             test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")


#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 30, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                           prefix = "SCT_snn_res.", assay = "SCT", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}


saveRDS(seu.obj, paste0("../output/s3/", outName,"_S3.rds"))


#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, reduction = "umap.integrated.harmony",
                    color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste0("../output/allCells/", experiment, "_QC_feats.png"), width = 9, height = 3)



seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "cellSource")


#load in independently analyzed and annotated data
seu.obj.bm <- readRDS("../output/s3/outputcombined_BM_res0.8_dims45_S3.rds")
seu.obj.bm <- loadMeta(seu.obj = seu.obj.bm, metaFile = "./metaData/BM_idents.csv", groupBy = "clusterID", metaAdd = "minorIdent")

seu.obj.cd34 <- readRDS("../output/s3/outputcombined_CD34_res0.3_dims35_S3.rds")
seu.obj.cd34 <- loadMeta(seu.obj = seu.obj.cd34, metaFile = "./metaData/CD34_idents.csv", groupBy = "clusterID", metaAdd = "minorIdent")


#modify the metadata rownames to ensure they match with the seu.obj cell barcodes
rownames(seu.obj.cd34@meta.data) <- paste0(seu.obj.cd34$orig.ident, "_", substr(rownames(seu.obj.cd34@meta.data), 1, nchar(rownames(seu.obj.cd34@meta.data))-4), "-1") 
rownames(seu.obj.bm@meta.data) <- paste0(seu.obj.bm$orig.ident, "_", substr(rownames(seu.obj.bm@meta.data), 1, nchar(rownames(seu.obj.bm@meta.data))-4), "-1") 

#extract the required cell barcodes and the corresponding idents -- store in named list
cd34_minorIdents <- seu.obj.cd34$minorIdent
names(cd34_minorIdents) <- rownames(seu.obj.cd34@meta.data)

#add the metadata to the seu.obj object
seu.obj <- AddMetaData(seu.obj, metadata = as.factor(c(cd34_minorIdents,seu.obj.bm$minorIdent)), col.name = "minorIdent")

#add a new metadata slot that contains info regarding where the cells came from
seu.obj$cellSource <- ifelse(grepl("BM",seu.obj$orig.ident), "BM","CD34")
seu.obj$cellSource <- droplevels(as.factor(seu.obj$cellSource))


seu.obj$minorIdent <- as.factor(ifelse(is.na(seu.obj$minorIdent), "remove", as.character(seu.obj$minorIdent)))
seu.obj <- subset(seu.obj, invert = T, minorIdent == "remove")
seu.obj$minorIdent <- droplevels(seu.obj$minorIdent)
outName <- "allCells"
#generate the plot to check label transfer was successful
### Fig 1a: plot inital cluster umap
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
ggsave(paste0("../output/", outName, "/", outName, "_minorIdent_UMAP.png"), width = 14, height = 7)

#save the processed object
saveRDS(seu.obj, "../output/s3/bm_cd34_analysis_231211_v5_integrated_res0.6_dims45_dist0.1_neigh10_S3.rds")



















