#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
library(monocle3)

##########
# TODO
# 1) modify pseudoDEG to exlcude filterTerm prior to running dgea

###########


############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Begin preliminary analysis   ######## <<<<<<<<<<<<<<
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set some output directives
reduction <- "umap.integrated.harmony"

#load in high resolution clstering and pseudotime
seu.obj <- readRDS("../output/s3/allCells_clean_highRes_integrated.harmony_res1.6_dims45_dist0.15_neigh20_S3.rds")
cds <- readRDS("../output/s3/240220_CDS_bm_cd34_removed_disconnected.rds")
ptime <- pseudotime(cds, reduction_method = "UMAP")
seu.obj$ptime <- ptime

#plot the high-res clustering results
pi <- DimPlot(seu.obj, 
              reduction = reduction, 
              group.by = "clusterID2_integrated.harmony",
              pt.size = 0.1,
              label = T,
              label.box = T,
              repel = F,
)
p <- cusLabels(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_clusID_highRes_UMAP.png"), width = 7, height = 7)

### Figure Xa - plot branch point 1
branch <- "branch1"
base_clus <- "5"
divergent_clus <- c("27","17")
groupBy <- "clusterID2_integrated.harmony"

#dge analysis between divergent_clus[1] and base_clus
btwnClusDEG(seu.obj = seu.obj, groupBy = groupBy, idents.1 = divergent_clus[1], idents.2 = base_clus, bioRep = "orig.ident", padj_cutoff = 0.01, lfcCut = 1, topn=c(20,20),
            minCells = 5, outDir = paste0("../output/", outName, "/", branch, "/"), title = paste0(divergent_clus[1], "v",base_clus), idents.1_NAME = divergent_clus[1], idents.2_NAME = base_clus, strict_lfc = T,
            returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 12, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3, filterTerm = c("^RPS", "^RPL", "^MT-")
                    )

#dge analysis between divergent_clus[2] and base_clus
btwnClusDEG(seu.obj = seu.obj, groupBy = groupBy, idents.1 = divergent_clus[2], idents.2 = base_clus, bioRep = "orig.ident", padj_cutoff = 0.01, lfcCut = 1, topn=c(20,20),
            minCells = 5, outDir = paste0("../output/", outName, "/", branch, "/"), title = paste0(divergent_clus[2], "v",base_clus), idents.1_NAME = divergent_clus[2], idents.2_NAME = base_clus, strict_lfc = T,
            returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 12, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3, filterTerm = c("^RPS", "^RPL", "^MT-")
                    )

#dge analysis between two divergent clusters
btwnClusDEG(seu.obj = seu.obj, groupBy = groupBy, idents.1 = divergent_clus[1], idents.2 = divergent_clus[2], bioRep = "orig.ident", padj_cutoff = 0.01, lfcCut = 1, topn=c(20,20),
            minCells = 5, outDir = paste0("../output/", outName, "/", branch, "/"), title = paste0( divergent_clus[1], "v",divergent_clus[2]), idents.1_NAME = divergent_clus[1], idents.2_NAME = divergent_clus[2], strict_lfc = T,
            returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 12, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3, filterTerm = c("^RPS", "^RPL", "^MT-")
                    )

#collect output
files <- list.files(path = paste0("../output/", outName, "/", branch, "/"), pattern=".csv", all.files=FALSE,
                        full.names=T)
df.list <- lapply(files, read.csv, header = T)
deg_res1 <- do.call(rbind, df.list) %>% mutate(grp = ifelse(log2FoldChange < 0, str_split(gs_base, "_", simplify = T)[ ,3], str_split(gs_base, "_", simplify = T)[ ,1])) %>% filter(grp == base_clus) %>% top_n(n = 20, wt = -padj)
deg_res2 <- do.call(rbind, df.list) %>% mutate(grp = ifelse(log2FoldChange < 0, str_split(gs_base, "_", simplify = T)[ ,3], str_split(gs_base, "_", simplify = T)[ ,1])) %>% filter(gs_base == paste0(divergent_clus[1], "_VS_", divergent_clus[2])) %>% group_by(grp) %>% top_n(n = 20, wt = -padj) %>% ungroup() %>% arrange(grp)

#plot heatmap
feats_forHeat <- c(deg_res2[deg_res2$grp == divergent_clus[2], ]$gene, deg_res1$gene, deg_res2[deg_res2$grp == divergent_clus[1], ]$gene)

#subset the data
Idents(seu.obj) <- groupBy
seu.obj.sub <- subset(seu.obj, idents = c(base_clus, divergent_clus))
seu.obj.sub$type <- factor(paste0(seu.obj.sub@meta.data[ ,groupBy], "--",seu.obj.sub$name))
seu.obj.sub$type <- factor(seu.obj.sub$type, levels = c(levels(seu.obj.sub$type)[grepl(paste0(divergent_clus[2], "--"), levels(seu.obj.sub$type))], 
                                                        levels(seu.obj.sub$type)[grepl(paste0(base_clus, "--"), levels(seu.obj.sub$type))], 
                                                        levels(seu.obj.sub$type)[grepl(paste0(divergent_clus[1], "--"), levels(seu.obj.sub$type))] 
                                                       ))

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
mat <- mat[ ,feats_forHeat]

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
png(file = paste0("../output/", outName, "/", outName, "_", branch, "_Heat_degs.png"), width=3000, height=4000, res=400)
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
seu.obj$pseudotimee <- ifelse(seu.obj@meta.data[ ,groupBy] %in% c(base_clus, divergent_clus), seu.obj$ptime, NA)
p <- FeaturePlot(seu.obj, features = "pseudotimee", reduction = reduction) 
p <- formatUMAP(p) + scale_color_continuous(na.value = "grey80")
ggsave(paste0("../output/", outName, "/", branch, "/", branch, "_pseudo_highlight_UMAP.png"), width = 7, height = 7)

#highlight region by cluster
colArray <- as.data.frame(matrix(levels(seu.obj@meta.data[ ,groupBy]), dimnames = list(NULL, "clusID")))
colArray$color <- gg_color_hue(nrow(colArray))
colArray <- colArray %>% mutate(high_col = ifelse(clusID %in% c(base_clus, divergent_clus), color,"grey"))

pi <- DimPlot(seu.obj, 
        reduction = reduction, 
        group.by = groupBy,
        cols = colArray$high_col,
        pt.size = 0.5,
        label = F,
        label.box = F
 )
umapHighLight <- formatUMAP(pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", branch, "_highlight_UMAP.png"), width = 7, height = 7)


#order cells by time then plot the genes in heatmap
seu.obj.sub <- subset(seu.obj.sub, downsample = min(table(droplevels(seu.obj@meta.data[ ,groupBy]))))
seu.obj.sub <- ScaleData(seu.obj.sub)
mat_scaled <- as.matrix(seu.obj.sub@assays$RNA$scale.data) #use scale data

order.df <- seu.obj.sub@meta.data[ ,c(groupBy, "ptime")] 
colnames(order.df)[1] <- "clus"                                                                              
left.cells <- order.df %>% filter(clus == divergent_clus[2]) %>% arrange(desc(ptime)) %>% rownames()
center.cells <- order.df %>% filter(clus == base_clus)  %>% rownames()
right.cells <- order.df %>% filter(clus == divergent_clus[1]) %>% arrange(desc(ptime)) %>% rownames()

#scale data
mat_scaled <- mat_scaled[rownames(mat_scaled) %in% feats_forHeat,]
mat_scaled <- mat_scaled[ ,match(c(left.cells,center.cells,right.cells), colnames(mat_scaled))]

#prep annotations
samp <- unique(seu.obj.sub$colz)
names(samp) <- unique(seu.obj.sub$name)
ha = HeatmapAnnotation(
    Sample = unlist(lapply(rownames(mat), function(x){strsplit(x,"--")[[1]][2]})),
    Cluster = unlist(lapply(rownames(mat), function(x){strsplit(x,"--")[[1]][1]})),
    border = TRUE,
    col = list(Sample = samp)
)

#plot the data -- this does not work as desired
png(file = paste0("../output/", outName, "/", branch, "_Heat_degs.png"), width=3000, height=4000, res=400)
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

#plot the data using Seurat DoHeatmap
Idents(seu.obj.sub) <- factor(droplevels(seu.obj.sub@meta.data[ ,groupBy]), levels =  c(divergent_clus[2], base_clus, divergent_clus[1]))
p <- DoHeatmap(seu.obj.sub, cells = c(left.cells,center.cells,right.cells), features = feats_forHeat)
ggsave(paste0("../output/", outName, "/", branch, "_seurat_heat.png"), width = 7, height = 7)

# #pick out TFs to highlight
# files <- list.files(path = paste0("../output/", outName, "/"), pattern=".csv", all.files=FALSE,
#                         full.names=T)
# df.list <- lapply(files, read.csv, header = T)
# df <- do.call(rbind, df.list)

# tfs <- read.csv("./metaData/allTFs_hg38.txt", col.names = "tfs")

# df <- df[df$gene %in% tfs$tfs, ]

# #gene expression by pseudotime
# features <- c("PIR","IRF4","TFEC", "NFE2") #only works with 1 feat rn
# colorBy <- "clusterID2_integrated.harmony" #seu metadata slot to color data by
# plotBy <- "ptime" #pseudotime data to use on the x axis -- match with the name of a metadata slot
# clus <- c("17","30")
# rmZeros <- F #option to remove zero values - probablaly can keep it at false
# bin.size <- 5

# seu.obj.sub <- subset(seu.obj, clusterID2_integrated.harmony %in% clus)
# mulit.plot <- lapply(features, function(x){
#     df <- FetchData(object = seu.obj.sub, vars = c('cellSource', 'name', colorBy, plotBy, x), layer = "data")
#     colnames(df)[3] <- "cluster"
#     colnames(df)[5] <- "feature"

#     if(rmZeros){
#         df <- df[df[ ,"feature"] != 0, ]
#     }

#     df <- df %>% na.omit() %>% group_by(cluster) %>% mutate(feat = x,
#                                                             ptime = (ptime - sort(ptime)[3]) / (rev(sort(ptime))[3] - sort(ptime)[3]) 
#                                                            ) %>% arrange(ptime) %>% mutate(avg = zoo::rollmean(feature, k = bin.size, fill = NA)) %>% ungroup()
#     return(df)
# })

# df <- do.call(rbind, mulit.plot) %>% na.omit() 

# p <- ggplot(data = df, aes(x = ptime, y = avg, color = feat, linetype = cluster)) + geom_smooth(se = FALSE)  + labs(title = "TFs by time", x = "Pseudotime", y = "Log2 normalized count") + guides(colour = guide_legend(ncol = 3)) 
# ggsave(paste0("../output/", outName, "/", outName, "_feats_byTime.png"), width = 7, height = 7)
