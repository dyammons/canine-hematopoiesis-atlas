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
                                "majorID", clusID_main, "name", "cellSource", "minorIdent"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                          "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                          "CD4", "MS4A1", "PPBP", "HBM")
)    



##################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Begin lineage analysis (monocle3)  ######## <<<<<<<<<<<<<<
##################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#create cell data set from raw data in Seurat object
cnt_mat <- seu.obj@assays$RNA@layers$counts
rownames(cnt_mat) <- rownames(seu.obj)
colnames(cnt_mat) <- colnames(seu.obj)
gene_annotation <- as.data.frame(rownames(seu.obj))
gene_annotation$gene_short_name <- rownames(seu.obj)
rownames(gene_annotation) <- rownames(seu.obj)
cds <- new_cell_data_set(cnt_mat,
                         cell_metadata = seu.obj@meta.data,
                         gene_metadata = gene_annotation)


#preprocess using monocle
cds <- preprocess_cds(cds, num_dim = 45)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds, umap.n_neighbors = 10)
p <- plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "minorIdent")
p <- formatUMAP(p)
ggsave(paste0("../output/", outName, "/", outName, "_cds_orig.ident_UMAP.png"), width = 7, height = 7)

#get the pseudotime
cds <- cluster_cells(cds, k = 30)
cds <- learn_graph(cds, use_partition = FALSE)
p <- plot_cells(cds,
           color_cells_by = "minorIdent",
                label_principal_points = T, 
           label_groups_by_cluster=FALSE,
           label_leaves=F,
           label_branch_points=F)
ggsave(paste0("../output/", outName, "/", outName, "_cds_branch_UMAP.png"), width = 7, height = 7)

cds <- order_cells(cds, root_pr_nodes = "Y_201")
p <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
ggsave(paste0("../output/", outName, "/", outName, "_cds_pseudo_UMAP.png"), width = 7, height = 7)



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

cds <- order_cells(cds, root_pr_nodes = "Y_567")
p <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
ggsave(paste0("../output/", outName, "/", outName, "_cds_pseudo_UMAP.png"), width = 7, height = 7)

#exctract pseudotime for use in Seurat
ptime <- pseudotime(cds, reduction_method = "UMAP")
seu.obj$ptime <- ptime


###################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Investigate critical branch points  ######## <<<<<<<<<<<<<<
###################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#create cell data set from raw data in Seurat object
#complete data visualization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_highRes_integrated.harmony"), 
                       final.dims = 45, final.res = 1.6, stashID = "clusterID2", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                       prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.harmony",
                       saveRDS = T, return_obj = T, returnFeats = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                    "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                    "CD4", "MS4A1", "PPBP","HBM")
                      )

seu.obj <- readRDS("../output/s3/allCells_clean_highRes_integrated.harmony_res1.6_dims45_dist0.1_neigh10_S3.rds")
#dge analysis between DC branch and precusor
btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_integrated.harmony", idents.1 = "17", idents.2 = "3", bioRep = "orig.ident", padj_cutoff = 0.05, lfcCut = 0.58, topn=c(20,20),
            minCells = 5, outDir = paste0("../output/", outName, "/"), title = "17v3", idents.1_NAME = "17", idents.2_NAME = "3", strict_lfc = T,
            returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = F, dwnSam = T, setSeed = 12, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3
                    )

#dge analysis between mono branch and precusor
btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_integrated.harmony", idents.1 = "30", idents.2 = "3", bioRep = "orig.ident", padj_cutoff = 0.05, lfcCut = 0.58, topn=c(20,20),
            minCells = 5, outDir = paste0("../output/", outName, "/"), title = "30v3", idents.1_NAME = "30", idents.2_NAME = "3", strict_lfc = T,
            returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = F, dwnSam = T, setSeed = 12, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3
                    )

#dge analysis between mono branch and precusor
btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_integrated.harmony", idents.1 = "30", idents.2 = "17", bioRep = "orig.ident", padj_cutoff = 0.05, lfcCut = 0.58, topn=c(20,20),
            minCells = 5, outDir = paste0("../output/", outName, "/"), title = "30v17", idents.1_NAME = "30", idents.2_NAME = "17", strict_lfc = T, 
            returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = F, dwnSam = T, setSeed = 12, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3
                    )

#collect output
files <- list.files(path = paste0("../output/", outName, "/"), pattern=".csv", all.files=FALSE,
                        full.names=T)
df.list <- lapply(files[c(1,3)], read.csv, header = T)
deg_res1 <- do.call(rbind, df.list) %>% mutate(grp = ifelse(log2FoldChange < 0, str_split(gs_base, "_", simplify = T)[ ,3], str_split(gs_base, "_", simplify = T)[ ,1])) %>% filter(grp == "3") %>% top_n(n = 20, wt = -padj)
deg_res2 <- read.csv(files[2]) %>% mutate(grp = ifelse(log2FoldChange < 0, str_split(gs_base, "_", simplify = T)[ ,3], str_split(gs_base, "_", simplify = T)[ ,1])) %>% group_by(grp) %>% top_n(n = 20, wt = -padj) %>% ungroup() %>% arrange(grp)


#plot heatmap
feats_forHeat <- c(deg_res2$gene[1:21], deg_res1$gene, deg_res2$gene[22:41])

#subset the data
seu.obj.sub <- subset(seu.obj, clusterID_integrated.harmony %in% c("3", "17", "30"))
seu.obj.sub$type <- paste0(seu.obj.sub$clusterID_integrated.harmony,"--",seu.obj.sub$name)

#extract metadata and data
metadata <- seu.obj.sub@meta.data
expression <- as.data.frame(t(seu.obj.sub@assays$RNA$data)) #use log noralized count
expression$anno_merge <- seu.obj.sub@meta.data[rownames(expression),]$type

#get cell type expression averages - do clus avg expression by sample
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL     

#scale data
mat <- scale(as.matrix(clusAvg_expression))
mat <- mat[,feats_forHeat]

#prep annotations
samp <- unique(seu.obj$colz)
names(samp) <- unique(seu.obj$name)
ha = HeatmapAnnotation(
    Sample = unlist(lapply(rownames(mat), function(x){strsplit(x,"--")[[1]][2]})),
    Cluster = unlist(lapply(rownames(mat), function(x){strsplit(x,"--")[[1]][1]})),
    border = TRUE,
    col = list(Sample = samp)
)

#plot the data
png(file = paste0("../output/", outName, "/", outName, "_dc_mono_branch_Heat_degs.png"), width=3000, height=4000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(t(mat), #name = "mat", #col = col_fun,
              name = "Scaled expression",
              cluster_rows = F,
              show_row_names=T,
              show_column_names=F,
              top_annotation = ha,
              col=viridis(option = "magma",100),
              cluster_columns = F
             )
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"))
dev.off()


#plot pseudotime in small region where deg analysis was completed
seu.obj$pseudotimee <- ifelse(seu.obj$clusterID_integrated.harmony %in% c("3", "17", "30"), seu.obj$pseudotime, NA)
p <- FeaturePlot(seu.obj, features = "pseudotimee", reduction = reduction) 
p <- formatUMAP(p) + scale_color_continuous(na.value = "grey80")
ggsave(paste0("../output/", outName, "/", outName, "_pseudo_highlight_UMAP.png"), width = 7, height = 7)

#highlight region by cluster
colArray <- as.data.frame(matrix(levels(seu.obj$clusterID_integrated.harmony), dimnames = list(NULL, "clusID")))
colArray$color <- gg_color_hue(nrow(colArray))
colArray <- colArray %>% mutate(high_col = ifelse(clusID %in% c("3", "17", "30"), color,"grey"))

pi <- DimPlot(seu.obj, 
        reduction = reduction, 
        group.by = "clusterID_integrated.harmony",
        cols = colArray$high_col,
        pt.size = 0.5,
        label = F,
        label.box = F
 )
umapHighLight <- formatUMAP(pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_pseudo_highlight_UMAP.png"), width = 7, height = 7)


#order cells by time then plot the genes in heatmap
seu.obj.sub <- subset(seu.obj.sub, downsample = min(table(droplevels(seu.obj.sub$clusterID_integrated.harmony))))
seu.obj.sub <- ScaleData(seu.obj.sub)
mat_scaled <- as.matrix(seu.obj.sub@assays$RNA$scale.data) #use scale data

order.df <- seu.obj.sub@meta.data[ ,c("clusterID_integrated.harmony", "pseudotime")] 
left.cells <- order.df %>% filter(clusterID_integrated.harmony == "17") %>% arrange(desc(pseudotime)) %>% rownames()
center.cells <- order.df %>% filter(clusterID_integrated.harmony == "3")  %>% rownames()
right.cells <- order.df %>% filter(clusterID_integrated.harmony == "30") %>% arrange(desc(pseudotime)) %>% rownames()

#scale data
mat_scaled <- mat_scaled[rownames(mat_scaled) %in% feats_forHeat,]
mat_scaled <- mat_scaled[ ,match(c(left.cells,center.cells,right.cells), colnames(mat_scaled))]

#prep annotations
samp <- unique(seu.obj$colz)
names(samp) <- unique(seu.obj$name)
ha = HeatmapAnnotation(
    Sample = unlist(lapply(rownames(mat), function(x){strsplit(x,"--")[[1]][2]})),
    Cluster = unlist(lapply(rownames(mat), function(x){strsplit(x,"--")[[1]][1]})),
    border = TRUE,
    col = list(Sample = samp)
)

#plot the data
png(file = paste0("../output/", outName, "/", outName, "_dc_mono_branch_Heat_degs.png"), width=3000, height=4000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(mat_scaled, #name = "mat", #col = col_fun,
              name = "Scaled expression",
              cluster_rows = F,
              show_row_names=T,
              show_column_names=F,
#               top_annotation = ha,
              col=viridis(option = "magma",100),
              cluster_columns = F
             )
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"))
dev.off()

Idents(seu.obj.sub) <- factor(droplevels(seu.obj.sub$clusterID_integrated.harmony), levels = c("17","3","30"))
p <- DoHeatmap(seu.obj.sub, cells = c(left.cells,center.cells,right.cells), features = feats_forHeat)
ggsave(paste0("../output/", outName, "/", outName, "_seurat_heat.png"), width = 7, height = 7)

#pick out TFs to highlight
files <- list.files(path = paste0("../output/", outName, "/"), pattern=".csv", all.files=FALSE,
                        full.names=T)
df.list <- lapply(files, read.csv, header = T)
df <- do.call(rbind, df.list)

tfs <- read.csv("./metaData/allTFs_hg38.txt", col.names = "tfs")

df <- df[df$gene %in% tfs$tfs, ]

#gene expression by pseudotime
features <- c("PIR","IRF4","TFEC", "NFE2") #only works with 1 feat rn
colorBy <- "clusterID2_integrated.harmony" #seu metadata slot to color data by
plotBy <- "ptime" #pseudotime data to use on the x axis -- match with the name of a metadata slot
clus <- c("17","30")
rmZeros <- F #option to remove zero values - probablaly can keep it at false
bin.size <- 5

seu.obj.sub <- subset(seu.obj, clusterID2_integrated.harmony %in% clus)
mulit.plot <- lapply(features, function(x){
    df <- FetchData(object = seu.obj.sub, vars = c('cellSource', 'name', colorBy, plotBy, x), layer = "data")
    colnames(df)[3] <- "cluster"
    colnames(df)[5] <- "feature"

    if(rmZeros){
        df <- df[df[ ,"feature"] != 0, ]
    }

    df <- df %>% na.omit() %>% group_by(cluster) %>% mutate(feat = x,
                                                            ptime = (ptime - sort(ptime)[3]) / (rev(sort(ptime))[3] - sort(ptime)[3]) 
                                                           ) %>% arrange(ptime) %>% mutate(avg = zoo::rollmean(feature, k = bin.size, fill = NA)) %>% ungroup()
    return(df)
})

df <- do.call(rbind, mulit.plot) %>% na.omit() 

p <- ggplot(data = df, aes(x = ptime, y = avg, color = feat, linetype = cluster)) + geom_smooth(se = FALSE)  + labs(title = "TFs by time", x = "Pseudotime", y = "Log2 normalized count") + guides(colour = guide_legend(ncol = 3)) #+ scale_colour_manual(values = colArray$colour)
ggsave(paste0("../output/", outName, "/", outName, "_feats_byTime.png"), width = 7, height = 7)

#  geom_point(data = df, aes(x = pseudotime, y = TCF4, colour = clusterID_integrated.harmony), alpha = 0.5) +

 + geom_point()





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























