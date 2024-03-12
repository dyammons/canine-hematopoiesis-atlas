#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
library(monocle3)

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
                markers = paste0("../output/viln/", outName, "/", outName, "_", clusID_main, "_gene_list.csv"),
                reduction = "umap.integrated.harmony",  
                colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                "majorID", clusID_main, "name", "cellSource", "minorIdent"), 
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


#load in the preprocessed object to remove an additional disconnected cell type
seu.obj <- readRDS("../output/s3/240201_bm_cd34_removed_disconnected.rds")
outName <- "allCells"
outName_full <- "240204_bm_cd34_clusterID_integrated.harmony"
clusID_main <- "clusterID_integrated.harmony"

#ensure proper cluster is being removed
Idents(seu.obj) <- "clusterID_integrated.harmony"
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated.harmony",
              group.by = "clusterID_integrated.harmony",
              pt.size = 0.1,
              label = T,
              cells.highlight = WhichCells(seu.obj, ident = "18"),
              label.box = F,
              repel = F
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_umap_highlight_eos.png"), width = 7, height = 7)

#remove eosinophils as they are disconnected and now easily descernable
seu.obj <- subset(seu.obj, invert = T, subset = clusterID_integrated.harmony == 18)

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

#override preivous object
saveRDS(seu.obj, "../output/s3/240201_bm_cd34_removed_disconnected.rds")

seu.obj <- readRDS("../output/s3/240201_bm_cd34_removed_disconnected.rds")
for (x in list("integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.15, n.neighbors = 20,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}

#override preivous object
saveRDS(seu.obj, "../output/s3/240201_bm_cd34_removed_disconnected.rds")

#evaluate clustering resolution and adjust the res accordingly
clusTree(seu.obj = seu.obj, 
                     dout = "../output/clustree/", outName = paste0(outName, "_highRes_integrated.harmony"), 
                     test_dims = 45, 
                     resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1)), 
                     algorithm = 3, 
                     prefix = "RNA_snn_res.",
         reduction = "integrated.harmony"
                    )

#complete clustering at a higher resolution -- DO NOT RUN; can skip to load in the object
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_highRes_integrated.harmony"), 
                       final.dims = 45, final.res = 1.4, stashID = "clusterID2", algorithm = 3, min.dist = 0.15, n.neighbors = 20,
                       prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.harmony",
                       saveRDS = T, return_obj = T, returnFeats = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                    "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                    "CD4", "MS4A1", "PPBP","HBM")
                      )

################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end additional preprocessing   ######## <<<<<<<<<<<<<<
################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Begin preliminary analysis   ######## <<<<<<<<<<<<<<
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#load in the processed data for further work up
seu.obj <- readRDS("../output/s3/allCells_clean_highRes_integrated.harmony_res1.4_dims45_dist0.15_neigh20_S3.rds")
table(seu.obj$clusterID2_integrated.harmony, seu.obj$minorIdent)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "celltype")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "majCol")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "majCol", metaAdd = "labCol")
seu.obj$clusterID2_integrated.harmony <- factor(seu.obj$clusterID2_integrated.harmony, levels = 0:max(as.numeric(names(table(seu.obj$clusterID2_integrated.harmony)))))

outName <- "allCells_clean"
outName_full <- "240201_bm_cd34_clusterID_integrated.harmony"
clusID_main <- "clusterID2_integrated.harmony"
reduction <- "umap.integrated.harmony"


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
                                "majorID", clusID_main, "name", "cellSource", "minorIdent"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                          "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                          "CD4", "MS4A1", "PPBP", "HBM")
)    


#create umap
colz.df <- read.csv("./metaData/allCells_ID_disconected_highRes.csv")
pi <- DimPlot(seu.obj, 
              reduction = reduction, 
              group.by = "clusterID2_integrated.harmony",
              pt.size = 0.1,
              cols = colz.df$majCol,
              label = T,
              label.box = T,
              repel = F,
) + NoLegend()
p <- cusLabels(plot = pi, labCol = colz.df$labCol, smallAxes = T, textSize = 5, size = 9)
ggsave(paste0("../output/", outName, "/", outName, "_Idents_UMAP.png"), width = 7, height = 7)

pi <- DimPlot(seu.obj, 
              reduction = reduction, 
              group.by = "minorIdent",
              split.by = "minorIdent",
              ncol = 5,
              pt.size = 0.1,
              label = F,
              label.box = F,
              repel = F,
) + NoLegend()
p <- formatUMAP(plot = pi)
ggsave(paste0("../output/", outName, "/", outName, "_minorIdents_UMAP.png"), width = 14, height = 14)


#set levels for ordering
seu.obj$clusterID2_integrated.harmony <- factor(seu.obj$clusterID2_integrated.harmony, levels = rev(c(11,5,4,
                                                                                                  24,19,16,
                                                                                                  12,27,15,25,
                                                                                                  9,23,10,21,6,8,0,3,2,26,
                                                                                                  28,20,14,13,7,22,1,18,17
                                                                                                  ))
                                               )

### Supp data - generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID2_integrated.harmony", numOfFeats = 24, outName = outName,
            outDir = paste0("../output/viln/", outName, "/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"),
            min.pct = 0.25, only.pos = T, returnViln = F)


colz <- colz.df$majCol[(1+rev(c(11,5,4,
                            24,19,16,
                            12,27,15,25,
                            9,23,10,21,6,8,0,3,2,26,
                            28,20,14,13,7,22,1,18,17
                            )))]


#make auto dot plot
pi <- autoDot(seu.integrated.obj = seu.obj, inFile = "../output/viln/allCells_clean/allCells_clean_clusterID2_integrated.harmony_gene_list.csv", groupBy = "clusterID2_integrated.harmony",
              MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1, n_feat = 3, filterTerm = "ENSCAFG", 
              lvls = rev(c(11,5,4,
                                                                                                  24,19,16,
                                                                                                  12,27,15,25,
                                                                                                  9,23,10,21,6,8,0,3,2,26,
                                                                                                  28,20,14,13,7,22,1,18,17
                                                                                                  ))
                    ) + #theme(legend.box="vertical", legend.position = "right") + 
scale_fill_manual(values = colz) 
ggsave(paste0("../output/", outName, "/", outName, "_autodot_branches.png"), width = 16, height = 10)


### Fig 3e - heatmap by cluster

seu.obj$type <- paste0(seu.obj$clusterID2_integrated.harmony, "-", seu.obj$majorID)



#extract metadata and data
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA$data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#get cell type expression averages - do clus avg expression by sample
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL     

df <- rownames(clusAvg_expression) %>% as.data.frame()
colnames(df) <- "type"
df$clus <- unlist(lapply(df$type, function(x){strsplit(x,"-")[[1]][1]}))
df$branch <- unlist(lapply(df$type, function(x){strsplit(x,"-")[[1]][2]}))
df$branch <- factor(df$branch, levels = c("HSC", "Erythroid", "Monocytic", "Lymphocytic", "Neutrophilic"))
df$clus <- factor(df$clus, levels = levels(seu.obj$clusterID2_integrated.harmony))
df$odor <- as.numeric(df$clus)
df <- df %>% arrange(odor)

#load in feats that define
sig.df <- read.csv("../output/viln/allCells_clean/allCells_clean_clusterID2_integrated.harmony_gene_list.csv")
sig.df <- sig.df %>% filter(p_val_adj < 0.01)
sig.df$cluster <- factor(sig.df$cluster, levels = levels(seu.obj$clusterID2_integrated.harmony))
sig.df <- sig.df %>% left_join(df[ ,2:3], by = c("cluster" = "clus"))
sig.df$odor <- as.numeric(sig.df$cluster)

#extract labels to plot - tweak to be by branch??
lab.df <- sig.df[!grepl("ENSCAF", sig.df$gene), ]
text_list <- rev(split(lab.df$gene, lab.df$branch))
text_list <- lapply(text_list, function(x){c(paste0(paste(x[1:3], collapse = ", "),","), 
                                             paste0(paste(x[4:6], collapse = ", "),","),
                                             paste0(paste(x[7:9], collapse = ", "),","),
                                             paste(x[10:12], collapse = ", ")
                                            )})

#finish ordering
sig.df <- sig.df[!duplicated(sig.df[,"gene"]), ]
sig.df <- sig.df %>% arrange(odor, desc(p_val_adj))
geneOrder <- rev(sig.df$gene)

#filter matrix for feats that define and scale by row
clusAvg_expression <- clusAvg_expression[ ,colnames(clusAvg_expression) %in% sig.df$gene]
mat_scaled <- t(apply(t(clusAvg_expression), 1, scale))
colnames(mat_scaled) <- rownames(clusAvg_expression)
length(geneOrder) == nrow(mat_scaled)
mat_scaled <- mat_scaled[ ,match(df$type, colnames(mat_scaled))]
mat_scaled <- mat_scaled[match(geneOrder, rownames(mat_scaled)), ]        

#set annotations
branch <- c("#AFDBF6", "#C8E1CC", "#C1C1C1", "#F97C7C", "#B18EEA")
names(branch) <- unique(seu.obj$majorID)
clus <- colz
names(clus) <- levels(seu.obj$clusterID2_integrated.harmony)
heat_col <- viridis(option = "magma",100)

ha <- HeatmapAnnotation(
    Branch = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"-")[[1]][2]})),
    Cluster = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"-")[[1]][1]})),
    border = TRUE,
    col = list(Branch = branch, Cluster = clus),
    show_legend = c(FALSE, TRUE),
    annotation_legend_param = list(
        Branch = list(
            nrow = 1
        ),
        Cluster = list(
            nrow = 3
        )
    )
)
ra <- rowAnnotation(foo = anno_empty(border = FALSE, width = max_text_width(unlist(text_list)) + unit(4, "mm")))


#plot the data
ht <- Heatmap(
    mat_scaled,
    name = "mat",
    cluster_rows = F,
    row_title_gp = gpar(fontsize = 24),
    show_row_names=F,
    col=heat_col,
    cluster_columns = F,
    top_annotation = ha,
    show_column_names = F,
    right_annotation = ra,
    cluster_row_slices=F,
    column_split = df$branch,
    row_split = factor(rev(sig.df$branch), levels = rev(sig.df[!duplicated(sig.df[,"branch"]), ]$branch)),
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
            title = "Scaled expression",
            direction = "horizontal"
        )
)

#save the plot
png(file = paste0("../output/", outName, "/", outName, "_heatinUp.png"), width=5000, height=4000, res=400)
par(mfcol=c(1,1))    
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), merge_legend = TRUE, heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom")

for(i in 1:length(levels(df$branch))) {
    decorate_annotation("Branch", slice = i, {
        grid.text(levels(df$branch)[i], just = "center")
    })
}

for(i in 1:length(rev(sig.df[!duplicated(sig.df[,"branch"]), ]$branch))) {
    decorate_annotation("foo", slice = i, {
        grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left",
                  gp = gpar(fontsize = 10))
    })
}

# #problem that the heatmap is sliced by branch not cluster
# for(i in 1:length(levels(df$clus))) {
#     decorate_annotation("Cluster", slice = i, {
#         grid.text(levels(df$clus)[i], just = "center")
#     })
# }

dev.off()















