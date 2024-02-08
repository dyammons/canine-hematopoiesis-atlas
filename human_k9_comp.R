#!/usr/bin/Rscript

#load custom functions & packages
message(paste0(Sys.time(), " INFO: loading packages and custom functions."))
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
library(ape)
library(ggtree)

#load in processed data
message(paste0(Sys.time(), " INFO: initating code."))
seu.obj <- readRDS("../output/s3/integrated.harmony_res0.8_dims45_dist0.2_neigh20_S3.rds")
outName <- "human_k9_comp"

#set metadata levels to compare
hu.meta <- "mantonID"
dog.meta <- "minorIdent"

# ### Fig - plot the umap with cell annotations - no label
# pi <- DimPlot(seu.obj, 
#               reduction = "umap.integrated", 
#               group.by = "mantonID",
#               pt.size = 0.1,
#               label = F,
#               label.box = F,
#               raster = F,
#               repel = F
# )
# p <- formatUMAP(plot = pi) + NoLegend()
# ggsave(paste0("../output/", outName, "/", outName, "_UMAP_hu_noLabs.png"), width = 7, height = 7)


# ### Fig - plot the umap with cell annotations - labeled
# pi <- DimPlot(seu.obj, 
#               reduction = "umap.integrated", 
#               group.by = "mantonID",
#               pt.size = 0.1,
#               label = T,
#               label.box = T,
#               raster = F,
#               repel = F
# )
# p <- formatUMAP(plot = pi) + NoLegend()
# ggsave(paste0("../output/", outName, "/", outName, "_UMAP_hu.png"), width = 7, height = 7)


# ### Data supplemental - generate violin plots of defining features
# vilnPlots(seu.obj = seu.obj, groupBy = "mantonID", numOfFeats = 24, outName = "manton_bm",
#                       outDir = paste0("../output/viln/", outName, "/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
#                       min.pct = 0.25, only.pos = T)


### Integrate cross-species

#read in processed k9 data
seu.obj.k9 <- readRDS("../output/s3/allCellslabTransfer_S3.rds")

#split then merge objects
message(paste0(Sys.time(), " INFO: splitting data from k9 and human."))
seu.list <- c(SplitObject(seu.obj.k9, split.by = "orig.ident"), SplitObject(seu.obj, split.by = "orig.ident"))
samNames <- names(seu.list)
seu.list <- lapply(1:length(seu.list), function(i){
    cnts <- seu.list[[i]]@assays$RNA@layers$counts
    rownames(cnts) <- rownames(seu.list[[i]])
    colnames(cnts) <- colnames(seu.list[[i]])
    CreateSeuratObject(cnts, project = samNames[i], assay = "RNA", meta.data = seu.list[[i]]@meta.data)
})

message(paste0(Sys.time(), " INFO: merging data from k9 and human."))
seu.merge <- merge(seu.list[1][[1]], y = seu.list[2:length(seu.list)],
                  add.cell.ids = samNames, 
                  project = "hu_k9_comp"
                 )
rm(seu.list)
gc()

#integrate the data
message(paste0(Sys.time(), " INFO: integrating data from k9 and human."))
seu.obj <- integrateData(din = NULL, pattern = NULL,
                          saveRDS = T, 
                          outName = outName,  dout = "../output/s2/",
                          orig.reduction = "pca",
                          normalization.method = "LogNormalize", 
                          method = "HarmonyIntegration",
                          indReClus = TRUE, seu.obj = seu.merge,
                          runAllMethods = FALSE
                        )
gc()

#complete data visualization & save the RDS file
message(paste0(Sys.time(), " INFO: data integration complete. compeleting dimension reduction and saving integrated object as a .rds file in ../s3/."))
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "integrated.harmony", 
                        final.dims = 45, final.res = 0.8, stashID = "clusterID", algorithm = 3, min.dist = 0.2, n.neighbors = 20,
                        prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated",
                        saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "HLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)
gc()

#update user
message(paste0(Sys.time(), " INFO: file saved. moving to plot hierchical clustering."))

#tag cell type labels with species labels to provide contrast
seu.obj$cellSource <- ifelse(grepl("Manton", seu.obj$orig.ident), "Human", "Canine")    
seu.obj$type <- ifelse(grepl("Canine", seu.obj$cellSource), paste0("c_", seu.obj[[dog.meta]]), paste0("hu_", seu.obj[[hu.meta]]))    

#extract data to plot
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$integrated@data))
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#get average expression by cluster
avg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(avg_expression) <- avg_expression$anno_merge
avg_expression$anno_merge <- NULL

#complete heirchical clustering and plot the data using two different methods
M <- (1- cor(t(avg_expression), method="pearson"))/2

outfile <- paste("./output/", outName, "/", outName, "_phylo.png", sep = "")
png(file = outfile, width=4000, height=4000, res=400)
par(mfcol=c(1,1))
hc <- hclust(as.dist(M),method="complete")
p <- plot.phylo(as.phylo(hc), type = "fan", cex = 0.8, no.margin = FALSE, tip.col = "black")
dev.off()

ggtree(as.phylo(hc), layout = "fan") + geom_tiplab() + NoLegend() + xlim(NA, 0.75)
ggsave(paste("./output/", outName, "/", outName, "_ggTree.png", sep = ""), width = 9, height = 9) 


### Use defining features to evaluate cell type gene signature overlap between species
message(paste0(Sys.time(), " INFO: evlauting gene expression intersections between species."))

#load in preprocessed data
seu.obj.hu <- readRDS("../output/s3/integrated.harmony_res0.8_dims45_dist0.2_neigh20_S3.rds")
seu.obj.dog <- seu.obj.k9
rm(seu.obj.k9)

seu.obj.hu[[hu.meta]] <- droplevels(seu.obj.hu[[hu.meta]])

#find overlapping features
interSECT <- rownames(seu.obj.hu)[rownames(seu.obj.hu) %in% rownames(seu.obj.dog)]

#run FindAllMarkers only including 1:1 orthologoues - on human
vilnPlots(seu.obj = seu.obj.hu, inFile = NULL, groupBy = hu.meta, numOfFeats = 24, outName = "human",
                      outDir = "", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL, returnViln = F, features = interSECT
                     )

#run FindAllMarkers only including 1:1 orthologoues - on dog
vilnPlots(seu.obj = seu.obj.dog, inFile = NULL, groupBy = dog.meta, numOfFeats = 24, outName = "dog",
                      outDir = "", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL, returnViln = F, features = interSECT
                     )


############ CODE BELOW IS UNTESTED ############


#read in the cell type gene lists
dog.df <- read.csv("../output/viln/dog_gene_list.csv")
human.df <- read.csv("../output/viln/human_gene_list.csv")

#only include features that were defining in each speices - this might be cheating
dog.df <- dog.df[dog.df$gene %in% interSECT, ]
human.df <- human.df[human.df$gene %in% interSECT, ]

#alternative approach - use stringent filtering to increase the liklihood of only using truely defining features
# dog.df <- dog.df[dog.df$pct.2 < 0.5 & dog.df$pct.1 > 0.6, ] #alternate filtering para
# human.df <- human.df[human.df$pct.2 < 0.5 & human.df$pct.1 > 0.6, ] #alternate filtering para

#calculate the number of overlaping DEGs
res <- lapply(unique(dog.df$cluster), function(x){
    
    dog.list <- dog.df[dog.df$cluster == x, ] %>% .$gene
    
    res_pre <- lapply(unique(human.df$cluster), function(y){
        human.list <- human.df[human.df$cluster == y, ] %>% .$gene
#         degList <- length(unique(c(dog.list,human.list)))
        intERsct <- length(dog.list[dog.list %in% human.list])
        names(intERsct) <- x
        
        return(intERsct)
    })
    
    res_pre <- do.call(rbind, res_pre)
    rownames(res_pre) <- unique(human.df$cluster)
    return(res_pre)
    
})
    
res1 <- do.call(cbind, res)


#calculate the percent overlaping DEGs
res <- lapply(unique(dog.df$cluster), function(x){
    
    dog.list <- dog.df[dog.df$cluster == x, ] %>% .$gene
    
    res_pre <- lapply(unique(human.df$cluster), function(y){
        human.list <- human.df[human.df$cluster == y, ] %>% .$gene
        intERsct <- length(dog.list[dog.list %in% human.list])/length(unique(c(dog.list,human.list)))
        names(intERsct) <- x
        
        return(intERsct)
    })
    
    res_pre <- do.call(rbind, res_pre)
    rownames(res_pre) <- unique(human.df$cluster)
    return(res_pre)
    
})
    
res2 <- do.call(cbind, res)

#use this code block to reorder when the time comes
# rowTarg <- c("Endothelial","Fibroblast" ,"Osteoblast","Cycling osteoblast",
#        "CD14_monocyte","NR4A3_Macrophage","TXNIP_Macrophage","FABP5_Macrophage",
#        "IFN-TAM","Pre-OC","Mature-OC","pDC","cDC2","Mast",
#        "CD8 T cell","CD4 T cell","IFN-sig T cell","NK cell",
#        "Plasma cell","B cell")


# colTarg <- c("Endothelial cell","Fibroblast","Osteoblast_1","Osteoblast_2","Osteoblast_3","Hypoxic_osteoblast",
#          "IFN-osteoblast","Osteoblast_cycling","M-MDSC","TIM","ANGIO_TAM","TAM_INT","TAM_ACT","LA-TAM_SPP2_hi",
#          "LA-TAM_C1QC_hi","IFN-TAM","Cycling_OC","CD320_OC","Mature_OC","pDC","cDC1","cDC2","mregDC",
#          "Mast cell","Neutrophil","CD8 T cell","CD4 T cell","T_IFN","T_cycling","NK","Plasma cell","B cell")

# res1 <- res1[match(rowTarg, rownames(res1)),]        
# res1 <- res1[ ,match(colTarg, colnames(res1))]

# res2 <- res2[match(rowTarg, rownames(res2)),]        
# res2 <- res2[ ,match(colTarg, colnames(res2))] 
       

#basic plot using pheatmap
png(file = paste0("../output/", outName, "/", outName, "_celltype_comp_pheatmap.png"), width=4000, height=4000, res=400)
par(mfcol=c(1,1))                    
pheatmap(t(res2), display_numbers = T, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


#pretty plot using ComplexHeatmap
small_mat <-  as.matrix(t(res1))
png(file = paste0("../output/", outName, "/", outName, "_celltype_comp_heatmap.png"), width=4000, height=4000, res=400)
par(mfcol=c(1,1))         

ht <- Heatmap(t(res2), #name = "mat", #col = col_fun,
              name = "% overlapping DEGs",
              cluster_rows = F,
              row_title = "Dog cell types",
              col=viridis(100),
              cluster_columns = F,
              column_title = "Human cell types",
              column_title_side = "bottom",
              column_names_rot = 45,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
              
              cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.0f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
              })

draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),heatmap_legend_side = "top")

dev.off()


#update user of completion
message(paste0(Sys.time(), " INFO: end of script!"))
