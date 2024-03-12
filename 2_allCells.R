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

################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end additional preprocessing   ######## <<<<<<<<<<<<<<
################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Begin preliminary analysis   ######## <<<<<<<<<<<<<<
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# #load in the object
seu.obj <- readRDS("../output/s3/240201_bm_cd34_removed_disconnected.rds")
outName <- "allCells_clean"
outName_full <- "240201_bm_cd34_clusterID_integrated.harmony"
clusID_main <- "clusterID_integrated.harmony"
reduction <- "umap.integrated.harmony"

#load metadata
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")
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
                                "majorID", clusID_main, "name", "cellSource", "minorIdent"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                          "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                          "CD4", "MS4A1", "PPBP", "HBM")
)    


##################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Begin lineage analysis (monocle3)  ######## <<<<<<<<<<<<<<
##################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### Do not run -- use code block below
### The noted code completes monocle precrocessing then run trajectory analysis.
### This yeilded similar results to the later apporach, so for sake of consitency
### the second approach was used as it runs monocle using the Seurat UMAP embeddings.

# #create cell data set from raw data in Seurat object
# cnt_mat <- seu.obj@assays$RNA@layers$counts
# rownames(cnt_mat) <- rownames(seu.obj)
# colnames(cnt_mat) <- colnames(seu.obj)
# gene_annotation <- as.data.frame(rownames(seu.obj))
# gene_annotation$gene_short_name <- rownames(seu.obj)
# rownames(gene_annotation) <- rownames(seu.obj)
# cds <- new_cell_data_set(cnt_mat,
#                          cell_metadata = seu.obj@meta.data,
#                          gene_metadata = gene_annotation)


# #preprocess using monocle
# cds <- preprocess_cds(cds, num_dim = 45)
# cds <- align_cds(cds, alignment_group = "orig.ident")
# cds <- reduce_dimension(cds, umap.n_neighbors = 10)
# p <- plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "minorIdent")
# p <- formatUMAP(p)
# ggsave(paste0("../output/", outName, "/", outName, "_cds_orig.ident_UMAP.png"), width = 7, height = 7)

# #get the pseudotime
# cds <- cluster_cells(cds, k = 30)
# cds <- learn_graph(cds, use_partition = FALSE)
# p <- plot_cells(cds,
#            color_cells_by = "minorIdent",
#                 label_principal_points = T, 
#            label_groups_by_cluster=FALSE,
#            label_leaves=F,
#            label_branch_points=F)
# ggsave(paste0("../output/", outName, "/", outName, "_cds_branch_UMAP.png"), width = 7, height = 7)

# cds <- order_cells(cds, root_pr_nodes = "Y_201")
# p <- plot_cells(cds,
#            color_cells_by = "pseudotime",
#            label_cell_groups=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,
#            graph_label_size=1.5)
# ggsave(paste0("../output/", outName, "/", outName, "_cds_pseudo_UMAP.png"), width = 7, height = 7)


#create cds object that retains the Seurat reductions
seu.obj@reductions$umap <- seu.obj@reductions$umap.integrated.harmony
cds <- SeuratWrappers::as.cell_data_set(seu.obj,
                                        assay = "RNA",
                                        reductions = "umap",
                                        default.reduction = "umap",
                                        graph = "RNA_snn",
                                        group.by = "clusterID_integrated.harmony")
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)
p <- plot_cells(cds,
           color_cells_by = "minorIdent",
                label_principal_points = T, 
           label_groups_by_cluster=FALSE,
           label_leaves=F,
           label_branch_points=F)
ggsave(paste0("../output/", outName, "/", outName, "_cds_branch_UMAP.png"), width = 7, height = 7)

cds <- order_cells(cds, root_pr_nodes = "Y_305")
p <- plot_cells(cds,
                color_cells_by = "pseudotime",
               label_cell_groups=FALSE,
               label_leaves=FALSE,
                show_trajectory_graph = T,
              label_branch_points = T,
              label_roots = FALSE,
           graph_label_size=1.5)
formatUMAP(p) + theme(axis.title = element_blank(),
                             panel.border = element_blank(),
                             plot.margin = unit(c(-7, -7, -7, -7), "pt")
                            ) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_cds_pseudo_UMAP.png"), width = 7, height = 7)

#look for degs overtime
deg_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
pr_deg_ids <- row.names(subset(deg_res, q_value < 0.01) %>% arrange(q_value))

p <- plot_cells(cds, genes=pr_deg_ids[1:4],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
ggsave(paste0("../output/", outName, "/", outName, "_cds_degs_UMAP.png"), width = 7, height = 7)

features <- pr_deg_ids[5:8]
p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 2, features = features, reduction = "umap.integrated.harmony",
                    color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste0("../output/", outName, "/", outName, "_cds_degs.png"), width = 7, height = 7)


cds <- preprocess_cds(cds, num_dim = 45) #find alternative to running pca again... although it seems to work ok
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble(cell = row.names(colData(cds)), 
                                cell_group = colData(cds)$majorID)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- str_c("Module ", row.names(agg_mat))

png(file = paste0("../output/", outName, "/", outName, "_modules_degs.png"), width = 1500, height = 3000, res = 400)
par(mfcol=c(1,1))         
pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
dev.off()

features <- gene_module_df %>% filter(module == 3) %>% pull(id) %>% head(4)
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 4, features = features, reduction = "umap.integrated.harmony",
                    color = "black", order = F, pt.size = 0.0000001, title.size = 10, noLegend = T)
ggsave(paste0("../output/", outName, "/", outName, "_cds_degs.png"), width = 8, height = 2)



AFD_genes <- c("gcy-8", "dac-1", "oig-8")
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$cell.type %in% c("AFD")]
AFD_lineage_cds <- order_cells(AFD_lineage_cds)

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5)

#save the .cds object
saveRDS(cds, file = "../output/s3/240220_CDS_bm_cd34_removed_disconnected.rds")


###################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Investigate critical branch points  ######## <<<<<<<<<<<<<<
###################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### This analysis is completed in sub-scripts 4_dge*.R
### The below code block creates a high-resolution clustered
### Seurat object.

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

#load in the processed data for further work up
# seu.obj <- readRDS("../output/s3/allCells_clean_highRes_integrated.harmony_res1.6_dims45_dist0.15_neigh20_S3.rds")
seu.obj <- readRDS("../output/s3/allCells_clean_highRes_integrated.harmony_res1.4_dims45_dist0.15_neigh20_S3.rds")
table(seu.obj$clusterID2_integrated.harmony, seu.obj$minorIdent)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "celltype")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "majCol")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "majCol", metaAdd = "labCol")
outName <- "allCells_clean"
outName_full <- "240201_bm_cd34_clusterID_integrated.harmony"
clusID_main <- "clusterID_integrated.harmony"
reduction <- "umap.integrated.harmony"


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


# ###################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# #######   Begin lineage analysis (slingshot)  ######## <<<<<<<<<<<<<<
# ###################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


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























