#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
library(monocle3)

#load in preprocessed data
seu.obj <- readRDS("../output/s3/allCells_clean_highRes_integrated.harmony_res1.4_dims45_dist0.15_neigh20_S3.rds")
table(seu.obj$clusterID2_integrated.harmony, seu.obj$minorIdent)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "celltype")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "majCol")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "majCol", metaAdd = "labCol")
seu.obj$clusterID2_integrated.harmony <- factor(seu.obj$clusterID2_integrated.harmony, levels = 0:max(as.numeric(names(table(seu.obj$clusterID2_integrated.harmony)))))

#set output and plotting parameters
outName <- "allCells_clean"
outName_full <- "240201_bm_cd34_clusterID_integrated.harmony"
clusID_main <- "clusterID2_integrated.harmony"
reduction <- "umap.integrated.harmony"


######################################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Begin pseudotime analysis (monocle3)  ######## <<<<<<<<<<<<<<
######################################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


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
seu.obj@reductions$umap <- seu.obj@reductions[[reduction]]
cds <- SeuratWrappers::as.cell_data_set(seu.obj,
                                        assay = "RNA",
                                        reductions = "umap",
                                        default.reduction = "umap",
                                        graph = "RNA_snn",
                                        group.by = clusID_main)

#cluster cells and learn graph to prep data
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)

#plot data to ID name of root node 
p <- plot_cells(cds,
           color_cells_by = "celltype",
                label_principal_points = T, 
           label_groups_by_cluster=FALSE,
           label_leaves=F,
           label_branch_points=F)
ggsave(paste0("../output/", outName, "/", outName, "_cds_branch_UMAP.png"), width = 7, height = 7)

#set node and order cells based on time
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

#use glm to find degs overtime -- takes ~5 min to run, can skip and load in .csv with results
# deg_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
# write.csv(deg_res, file = paste0("../output/", outName, "/", outName, "_graph_test_res.csv"))
deg_res <- read.csv(file = paste0("../output/", outName, "/", outName, "_graph_test_res.csv"), row.names = 1)
pr_deg_ids <- row.names(subset(deg_res, q_value < 0.01) %>% arrange(q_value))

#group genes in modules and investigate DE
cds <- preprocess_cds(cds, num_dim = 45) #find alternative to running pca again? although it seems to work ok
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble(cell = row.names(colData(cds)), 
                        cell_group = colData(cds)$majorID
                       )
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- str_c("Module_", row.names(agg_mat))

png(file = paste0("../output/", outName, "/", outName, "_modules_degs.png"), width = 1500, height = 4000, res = 400)
par(mfcol = c(1,1))         
pheatmap(agg_mat, scale="column")
dev.off()

#save the .cds object
saveRDS(cds, file = "../output/s3/240314_CDS_bm_cd34_highRES.rds")


###################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   End pseudotime analysis (monocle3)  ######## <<<<<<<<<<<<<<
###################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# NOTE: can skip the monocle code and load in from checkpoint!
# Load from "../output/s3/240314_CDS_bm_cd34_highRES.rds"

############################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Beging plotting (monocle3)  ######## <<<<<<<<<<<<<<
############################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Resume
cds <- readRDS("../output/s3/240314_CDS_bm_cd34_highRES.rds")
ptime <- pseudotime(cds, reduction_method = "UMAP")
seu.obj$ptime <- ptime

### Run analysis on neutrophilic branch

features <- c("MPO","CAMP","CD4","SELL") #only works with 1 feat rn
colorBy <- "clusterID2_integrated.harmony" #seu metadata slot to color data by
plotBy <- "ptime" #pseudotime data to use on the x axis -- match with the name of a metadata slot
clus <- c("Neutrophilic")
rmZeros <- F #option to remove zero values - probablaly can keep it at false
bin.size <- 5


seu.obj.sub <- subset(seu.obj, majorID %in% clus)
single.plot <- lapply(features, function(x){
    df <- FetchData(object = seu.obj.sub, vars = c('cellSource', 'name', colorBy, plotBy, x), layer = "data")
    colnames(df)[3] <- "cluster"
    colnames(df)[5] <- "feature"

    if(rmZeros){
        df <- df[df[ ,"feature"] > 0, ]
    }

    df <- df %>% na.omit() %>% mutate(feat = x,
                                      ptime = (ptime - sort(ptime)[3]) / (rev(sort(ptime))[3] - sort(ptime)[3]) 
                                     ) %>% arrange(ptime) %>% mutate(avg = zoo::rollmean(feature, k = bin.size, fill = NA))
    return(df)
})

df <- do.call(rbind, single.plot) %>% na.omit() 

colz.df <- read.csv("./metaData/allCells_ID_disconected_highRes.csv")
colz.df <- colz.df %>% filter(majorID == clus) %>% arrange(desc(order))
df$feat <- factor(df$feat, levels = features)
p <- ggplot(data = df, aes(x = ptime, y = avg, color = cluster)) + theme_classic() + facet_wrap(vars(feat), ncol = 1) + geom_point() +
geom_smooth(se = FALSE, colour = "black")  + labs(x = "Pseudotime", y = "Log2 normalized count") + NoLegend() + scale_colour_manual(values = colz.df$majCol) +
scale_x_continuous(breaks = c(0,1))
ggsave(paste0("../output/", outName, "/", outName, "_feats_byTime.png"), width = 4, height = 6)


#create legend
cusLeg(legend = colz.df, colz = 2, rowz = 5, clusLabel = "clusterID2_integrated.harmony", legLabel = "celltype", groupLabel = "majorID", colorz = "majCol", dotSize = 8,
                   groupBy = "", sortBy = "order", labCol = "labCol", headerSize = 6, #add if len labCol == 1 then add col to the df
                   cellHeader = T, bump = 2, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 5, compress_x = 0.75, titleOrder = NULL, spaceBtwnCols = 0.4, breakGroups = T, horizontalWrap = F, returnData = F, overrideData = NULL
                     )
ggsave(paste0("../output/", outName, "/", outName, "_legend.png"), width = 4, height = 4)

#create umap with branch pseudotime
ptime.sub <- df %>% tibble::rownames_to_column(., "rownames") %>% filter(feat == features[1]) %>% pull(ptime, name = rownames)
seu.obj <- AddMetaData(seu.obj, metadata = ptime.sub, col.name = "ptime.sub")
p <- prettyFeats(seu.obj = seu.obj, 
                 reduction = "umap.integrated.harmony",
                 nrow = 1, ncol = 1, 
                 features = "ptime.sub", color = "black", 
                 order = F, titles = NA,
                 noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

p <- FeaturePlot(seu.obj, features = "ptime.sub", reduction = reduction) 
pi <- formatUMAP(p, smallAxes = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 5)) + scale_colour_viridis(na.value = "grey")
ggsave(paste0("../output/", outName, "/", outName, "_pTime.png"), width = 7, height = 7)

#make ridge plot
seu.obj.sub <- AddMetaData(seu.obj.sub, metadata = ptime.sub, col.name = "ptime.sub")
seu.obj.sub$celltype <- factor(seu.obj.sub$celltype, levels = rev(colz.df$celltype))
p <- RidgePlot(seu.obj.sub, features = "ptime.sub",
               group.by = "celltype",
               cols = colz.df$majCol) + NoLegend() + theme(axis.title = element_blank(),
                                                           plot.title = element_blank())
ggsave(paste0("../output/", outName, "/", outName, "_ridge.png"), width = 7, height = 4)


### Run analysis on Erythroid branch

features <- c("SELP","ADGRG1","HBM","NOSTRIN") #only works with 1 feat rn
colorBy <- "clusterID2_integrated.harmony" #seu metadata slot to color data by
plotBy <- "ptime" #pseudotime data to use on the x axis -- match with the name of a metadata slot
clus <- "Erythroid"
rmZeros <- F #option to remove zero values - probablaly can keep it at false
bin.size <- 5

#use glm to find degs overtime -- takes ~2 min to run, can skip and load in .csv with results
cds.sub <- cds[ ,colData(cds)$majorID %in% clus]
cds.sub <- cluster_cells(cds.sub)
cds.sub <- learn_graph(cds.sub, use_partition = FALSE)
deg_res <- graph_test(cds.sub, neighbor_graph = "principal_graph", cores = 4)
write.csv(deg_res, file = paste0("../output/", outName, "/", clus, "_graph_test_res.csv"))


seu.obj.sub <- subset(seu.obj, majorID %in% clus)
single.plot <- lapply(features, function(x){
    df <- FetchData(object = seu.obj.sub, vars = c('cellSource', 'name', colorBy, plotBy, x), layer = "data")
    colnames(df)[3] <- "cluster"
    colnames(df)[5] <- "feature"

    if(rmZeros){
        df <- df[df[ ,"feature"] > 0, ]
    }

    df <- df %>% na.omit() %>% mutate(feat = x,
                                      ptime = (ptime - sort(ptime)[3]) / (rev(sort(ptime))[3] - sort(ptime)[3]) 
                                     ) %>% arrange(ptime) %>% mutate(avg = zoo::rollmean(feature, k = bin.size, fill = NA))
    return(df)
})

df <- do.call(rbind, single.plot) %>% na.omit() 

colz.df <- read.csv("./metaData/allCells_ID_disconected_highRes.csv")
colz.df <- colz.df %>% filter(majorID == clus) %>% arrange(desc(order))
df$feat <- factor(df$feat, levels = features)
p <- ggplot(data = df, aes(x = ptime, y = avg, color = cluster)) + theme_classic() + facet_wrap(vars(feat), ncol = 1, scales = "free_y") + geom_point() +
geom_smooth(se = FALSE, colour = "black")  + labs(x = "Pseudotime", y = "Log2 normalized count") + NoLegend() + scale_colour_manual(values = colz.df$majCol) +
scale_x_continuous(breaks = c(0,1))
ggsave(paste0("../output/", outName, "/", outName, "_feats_byTime.png"), width = 4, height = 6)


#create legend
cusLeg(legend = colz.df, colz = 2, rowz = 2, clusLabel = "clusterID2_integrated.harmony", legLabel = "celltype", groupLabel = "majorID", colorz = "majCol", dotSize = 8,
                   groupBy = "", sortBy = "order", labCol = "labCol", headerSize = 6, #add if len labCol == 1 then add col to the df
                   cellHeader = T, bump = 2, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 5, compress_x = 0.75, titleOrder = NULL, spaceBtwnCols = 0.4, breakGroups = T, horizontalWrap = T, returnData = F, overrideData = NULL
                     )
ggsave(paste0("../output/", outName, "/", outName, "_legend.png"), width = 4, height = 4)

#create umap with branch pseudotime
seu.obj$ptime.sub <- NULL
ptime.sub <- df %>% tibble::rownames_to_column(., "rownames") %>% filter(feat == features[1]) %>% pull(ptime, name = rownames)
seu.obj <- AddMetaData(seu.obj, metadata = ptime.sub, col.name = "ptime.sub")
p <- prettyFeats(seu.obj = seu.obj, 
                 reduction = "umap.integrated.harmony",
                 nrow = 1, ncol = 1, 
                 features = "ptime.sub", color = "black", 
                 order = F, titles = NA,
                 noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

p <- FeaturePlot(seu.obj, features = "ptime.sub", reduction = reduction) 
pi <- formatUMAP(p, smallAxes = T) & scale_colour_viridis(na.value = "grey")
ggsave(paste0("../output/", outName, "/", outName, "_pTime.png"), width = 7, height = 7)

#make ridge plot
seu.obj.sub <- AddMetaData(seu.obj.sub, metadata = ptime.sub, col.name = "ptime.sub")
seu.obj.sub$celltype <- factor(seu.obj.sub$celltype, levels = rev(colz.df$celltype))
p <- RidgePlot(seu.obj.sub, features = "ptime.sub",
               group.by = "celltype",
               cols = colz.df$majCol) + NoLegend() + theme(axis.title = element_blank(),
                                                           plot.title = element_blank())
ggsave(paste0("../output/", outName, "/", outName, "_ridge.png"), width = 7, height = 4)


### Run analysis on Monocytic branch

features <- c("MAFB","LY86","CD300H","SPATS2L") #only works with 1 feat rn
colorBy <- "clusterID2_integrated.harmony" #seu metadata slot to color data by
plotBy <- "ptime" #pseudotime data to use on the x axis -- match with the name of a metadata slot
clus <- "Monocytic"
rmZeros <- F #option to remove zero values - probablaly can keep it at false
bin.size <- 5

seu.obj.sub <- subset(seu.obj, majorID %in% clus)
single.plot <- lapply(features, function(x){
    df <- FetchData(object = seu.obj.sub, vars = c('cellSource', 'name', colorBy, plotBy, x), layer = "data")
    colnames(df)[3] <- "cluster"
    colnames(df)[5] <- "feature"

    if(rmZeros){
        df <- df[df[ ,"feature"] > 0, ]
    }

    df <- df %>% na.omit() %>% mutate(feat = x,
                                      ptime = (ptime - sort(ptime)[3]) / (rev(sort(ptime))[3] - sort(ptime)[3]) 
                                     ) %>% arrange(ptime) %>% mutate(avg = zoo::rollmean(feature, k = bin.size, fill = NA))
    return(df)
})

df <- do.call(rbind, single.plot) %>% na.omit() 

colz.df <- read.csv("./metaData/allCells_ID_disconected_highRes.csv")
colz.df <- colz.df %>% filter(majorID == clus) %>% arrange(desc(order))
df$feat <- factor(df$feat, levels = features)
p <- ggplot(data = df, aes(x = ptime, y = avg, color = cluster)) + theme_classic() + facet_wrap(vars(feat), ncol = 1, scales = "free_y") + geom_point() +
geom_smooth(se = FALSE, colour = "black")  + labs(x = "Pseudotime", y = "Log2 normalized count") + NoLegend() + scale_colour_manual(values = colz.df$majCol) +
scale_x_continuous(breaks = c(0,1))
ggsave(paste0("../output/", outName, "/", outName, "_feats_byTime.png"), width = 4, height = 6)


#create legend
cusLeg(legend = colz.df, colz = 2, rowz = 2, clusLabel = "clusterID2_integrated.harmony", legLabel = "celltype", groupLabel = "majorID", colorz = "majCol", dotSize = 8,
                   groupBy = "", sortBy = "order", labCol = "labCol", headerSize = 6, #add if len labCol == 1 then add col to the df
                   cellHeader = T, bump = 2, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 5, compress_x = 0.75, titleOrder = NULL, spaceBtwnCols = 0.4, breakGroups = T, horizontalWrap = T, returnData = F, overrideData = NULL
                     )
ggsave(paste0("../output/", outName, "/", outName, "_legend.png"), width = 4, height = 4)

#create umap with branch pseudotime
seu.obj$ptime.sub <- NULL
ptime.sub <- df %>% tibble::rownames_to_column(., "rownames") %>% filter(feat == features[1]) %>% pull(ptime, name = rownames)
seu.obj <- AddMetaData(seu.obj, metadata = ptime.sub, col.name = "ptime.sub")
p <- prettyFeats(seu.obj = seu.obj, 
                 reduction = "umap.integrated.harmony",
                 nrow = 1, ncol = 1, 
                 features = "ptime.sub", color = "black", 
                 order = F, titles = NA,
                 noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

p <- FeaturePlot(seu.obj, features = "ptime.sub", reduction = reduction) 
pi <- formatUMAP(p, smallAxes = T) & scale_colour_viridis(na.value = "grey")
ggsave(paste0("../output/", outName, "/", outName, "_pTime.png"), width = 7, height = 7)

#make ridge plot
seu.obj.sub <- AddMetaData(seu.obj.sub, metadata = ptime.sub, col.name = "ptime.sub")
seu.obj.sub$celltype <- factor(seu.obj.sub$celltype, levels = rev(colz.df$celltype))
p <- RidgePlot(seu.obj.sub, features = "ptime.sub",
               group.by = "celltype",
               cols = colz.df$majCol) + NoLegend() + theme(axis.title = element_blank(),
                                                           plot.title = element_blank())
ggsave(paste0("../output/", outName, "/", outName, "_ridge.png"), width = 7, height = 4)


### Run analysis on Lymphocytic branch

features <- c("IL7R","TOX2","PAX5","MS4A1") #only works with 1 feat rn
colorBy <- "clusterID2_integrated.harmony" #seu metadata slot to color data by
plotBy <- "ptime" #pseudotime data to use on the x axis -- match with the name of a metadata slot
clus <- "Lymphocytic"
rmZeros <- F #option to remove zero values - probablaly can keep it at false
bin.size <- 5

seu.obj.sub <- subset(seu.obj, majorID %in% clus)
single.plot <- lapply(features, function(x){
    df <- FetchData(object = seu.obj.sub, vars = c('cellSource', 'name', colorBy, plotBy, x), layer = "data")
    colnames(df)[3] <- "cluster"
    colnames(df)[5] <- "feature"

    if(rmZeros){
        df <- df[df[ ,"feature"] > 0, ]
    }

    df <- df %>% na.omit() %>% mutate(feat = x,
                                      ptime = (ptime - sort(ptime)[3]) / (rev(sort(ptime))[3] - sort(ptime)[3]) 
                                     ) %>% arrange(ptime) %>% mutate(avg = zoo::rollmean(feature, k = bin.size, fill = NA))
    return(df)
})

df <- do.call(rbind, single.plot) %>% na.omit() 

colz.df <- read.csv("./metaData/allCells_ID_disconected_highRes.csv")
colz.df <- colz.df %>% filter(majorID == clus) %>% arrange(desc(order))
df$feat <- factor(df$feat, levels = features)
p <- ggplot(data = df, aes(x = ptime, y = avg, color = cluster)) + theme_classic() + facet_wrap(vars(feat), ncol = 1, scales = "free_y") + geom_point() +
geom_smooth(se = FALSE, colour = "black")  + labs(x = "Pseudotime", y = "Log2 normalized count") + NoLegend() + scale_colour_manual(values = colz.df$majCol) +
scale_x_continuous(breaks = c(0,1))
ggsave(paste0("../output/", outName, "/", outName, "_feats_byTime.png"), width = 4, height = 6)


#create legend
cusLeg(legend = colz.df, colz = 2, rowz = 5, clusLabel = "clusterID2_integrated.harmony", legLabel = "celltype", groupLabel = "majorID", colorz = "majCol", dotSize = 8,
                   groupBy = "", sortBy = "order", labCol = "labCol", headerSize = 6, #add if len labCol == 1 then add col to the df
                   cellHeader = T, bump = 2, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 5, compress_x = 0.75, titleOrder = NULL, spaceBtwnCols = 0.4, breakGroups = T, horizontalWrap = T, returnData = F, overrideData = NULL
                     )
ggsave(paste0("../output/", outName, "/", outName, "_legend.png"), width = 4, height = 4)

#create umap with branch pseudotime
seu.obj$ptime.sub <- NULL
ptime.sub <- df %>% tibble::rownames_to_column(., "rownames") %>% filter(feat == features[1]) %>% pull(ptime, name = rownames)
seu.obj <- AddMetaData(seu.obj, metadata = ptime.sub, col.name = "ptime.sub")
p <- prettyFeats(seu.obj = seu.obj, 
                 reduction = "umap.integrated.harmony",
                 nrow = 1, ncol = 1, 
                 features = "ptime.sub", color = "black", 
                 order = F, titles = NA,
                 noLegend = T) + theme(legend.position = 'bottom') + guides(color = guide_colourbar(barwidth = 1)) + plot_layout(guides = "collect") & scale_colour_viridis(na.value="grey")

p <- FeaturePlot(seu.obj, features = "ptime.sub", reduction = reduction) 
pi <- formatUMAP(p, smallAxes = T) & scale_colour_viridis(na.value = "grey")
ggsave(paste0("../output/", outName, "/", outName, "_pTime.png"), width = 7, height = 7)

#make ridge plot
seu.obj.sub <- AddMetaData(seu.obj.sub, metadata = ptime.sub, col.name = "ptime.sub")
seu.obj.sub$celltype <- factor(seu.obj.sub$celltype, levels = rev(colz.df$celltype))
p <- RidgePlot(seu.obj.sub, features = "ptime.sub",
               group.by = "celltype",
               cols = colz.df$majCol) + NoLegend() + theme(axis.title = element_blank(),
                                                           plot.title = element_blank())
ggsave(paste0("../output/", outName, "/", outName, "_ridge.png"), width = 7, height = 4)

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








