#!/usr/bin/Rscript

#load custom functions & packages
message(paste0(Sys.time(), " INFO: loading packages and custom functions."))
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
library(ape)
library(ggtree)

#load in processed data
message(paste0(Sys.time(), " INFO: initating code."))
seu.obj <- readRDS("../output/s3/manton_ref_res0.8_dims45_dist0.2_neigh20_S3.rds")
exclude <- unlist(lapply(1:8, function(x){
    namez <- lapply(c(1,3,5:8), function(y){
        namez <- paste0("MantonBM", x, "_HiSeq_",y)
    })
    namez <- unlist(namez)
}))
seu.obj <- subset(seu.obj, invert = T, subset = orig.ident %in% exclude)
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
cnts <- seu.obj.k9@assays$RNA$counts
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "dog",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))
seu.obj.k9 <- CreateSeuratObject(cnts, project = "humanConvert", assay = "RNA",
                                  min.cells = 0, min.features = 0, names.field = 1,
                                  names.delim = "_", meta.data = seu.obj.k9@meta.data)

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

rm(seu.obj.k9)
rm(seu.obj)
gc()

#ensure all objects have percent.mt stored in metadata as this will be regressed during integration
seu.list <- lapply(seu.list, PercentageFeatureSet, pattern = "^MT-", col.name = "percent.mt")

message(paste0(Sys.time(), " INFO: integrating data from k9 and human."))
seu.obj <- integrateData_v4(seu.list = seu.list, outDir = "../output/s2/", subName = "human_k9",
                            vars.to.regress = "percent.mt", saveRDS = T
                           )
rm(seu.list)
gc()


# seu.obj <- readRDS("../output/s2/human_k9_S2.rds")
#complete data visualization & save the RDS file
message(paste0(Sys.time(), " INFO: data integration complete. compeleting dimension reduction and saving integrated object as a .rds file in ../s3/."))
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "integrated_v4", 
                        final.dims = 45, final.res = 0.8, stashID = "clusterID", algorithm = 3, min.dist = 0.2, n.neighbors = 20,
                        prefix = "integrated_snn_res.", assay = "integrated", reduction = "pca",
                        saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "HLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)
gc()

#update user
message(paste0(Sys.time(), " INFO: file saved. moving to plot hierchical clustering."))


#reload in the integrated object
# seu.obj <- readRDS("../output/s3/integrated_v4_res0.8_dims45_dist0.2_neigh20_S3.rds")

#tag cell type labels with species labels to provide contrast
seu.obj$cellSource <- ifelse(grepl("Manton", seu.obj$orig.ident), "Human", "Canine")

#load in the processed data for further work up
seu.obj.k9_anno <- readRDS("../output/s3/allCells_clean_highRes_integrated.harmony_res1.6_dims45_dist0.15_neigh20_S3.rds")
seu.obj.k9_anno <- loadMeta(seu.obj = seu.obj.k9_anno, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "celltype")
df <- seu.obj.k9_anno@meta.data %>% rownames_to_column() %>% select(rowname, celltype)
df <- seu.obj@meta.data %>% rownames_to_column() %>% left_join(df, by = "rowname")
table(is.na(df$celltype))
df <- df %>% mutate(celltype = ifelse(cellSource == "Canine", 
                                      ifelse(is.na(celltype), as.character(minorIdent), as.character(celltype)), NA))
table(is.na(df$celltype))
table(seu.obj$cellSource)
rownames(df) <- df$rowname
df$rowname <- NULL
seu.obj@meta.data <- df

seu.obj$type <- as.factor(ifelse(grepl("Canine", seu.obj$cellSource), paste0("c_", seu.obj$celltype), paste0("hu_", seu.obj$mantonID)))
table(seu.obj$type)

#remove unwanted cell types (aka those that were not annotated)
Idents(seu.obj) <- "type"
seu.obj <- subset(seu.obj, invert = T, idents = levels(seu.obj$type)[table(seu.obj$type) < 23][1:16])
seu.obj$type <- droplevels(seu.obj$type)

#extract data to plot
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(as.matrix(seu.obj@assays$integrated@data)))
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#get average expression by cluster
avg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), sum)) %>% as.data.frame()
rownames(avg_expression) <- avg_expression$anno_merge
avg_expression$anno_merge <- NULL

#complete heirchical clustering and plot the data using two different methods
M <- (1- cor(t(avg_expression), method="pearson"))/2

outfile <- paste("../output/", outName, "/", outName, "_phylo.png", sep = "")
png(file = outfile, width=4000, height=4000, res=400)
par(mfcol=c(1,1))
hc <- hclust(as.dist(M),method="complete")
p <- plot.phylo(as.phylo(hc), type = "fan", cex = 0.8, no.margin = FALSE, tip.col = "black")
dev.off()

ggtree(as.phylo(hc), layout = "fan") + geom_tiplab() + NoLegend() + xlim(NA, 0.75)
ggsave(paste("../output/", outName, "/", outName, "_ggTree.png", sep = ""), width = 9, height = 9) 


### Use defining features to evaluate cell type gene signature overlap between species
message(paste0(Sys.time(), " INFO: evlauting gene expression intersections between species."))

#load in preprocessed data
seu.obj.hu <- subset(seu.obj, subset = cellSource == "Human")
seu.obj.hu <- JoinLayers(seu.obj.hu)
seu.obj.dog <- subset(seu.obj, subset = cellSource == "Canine")
seu.obj.dog <- JoinLayers(seu.obj.dog)

dog.meta <- "celltype"

seu.obj.hu[[hu.meta]] <- droplevels(seu.obj.hu[[hu.meta]])
seu.obj.dog[[dog.meta]] <- droplevels(seu.obj.dog[[dog.meta]])

#find overlapping features to only use features expressed in both species
interSECT <- rownames(seu.obj.hu)[rownames(seu.obj.hu) %in% rownames(seu.obj.dog)]

#run FindAllMarkers - on human
vilnPlots(seu.obj = seu.obj.hu, inFile = NULL, groupBy = hu.meta, numOfFeats = 24, outName = "human",
                      outDir = "../output/viln/human/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), 
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL, returnViln = F, features = interSECT
                     )

#run FindAllMarkers - on dog
vilnPlots(seu.obj = seu.obj.dog, inFile = NULL, groupBy = dog.meta, numOfFeats = 24, outName = "dog",
                      outDir = "../output/viln/dog/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"),
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL, returnViln = F, features = interSECT
                     )


#read in the cell type gene lists
dog.df <- read.csv("../output/viln/dog/dog_celltype_gene_list.csv") # need to convert these to human gene sybmols
human.df <- read.csv("../output/viln/human/human_mantonID_gene_list.csv")

#check the data quality
dog.df %>% group_by(cluster) %>% summarize(nn = n()) %>% summarize(min = min(nn))
human.df %>% group_by(cluster) %>% summarize(nn = n()) %>% summarize(min = min(nn))

#calc the jaccard similarity index
res <- lapply(unique(dog.df$cluster), function(x){
    dog.list <- dog.df[dog.df$cluster == x, ] %>% .$gene
    
    res_pre <- lapply(unique(human.df$cluster), function(y){
        human.list <- human.df[human.df$cluster == y, ] %>% .$gene
        
        interSect <- length(intersect(human.list, dog.list)) 
        uni <- length(human.list) + length(dog.list) - interSect 
        JaccardIndex <- interSect/uni
        names(JaccardIndex) <- x
        
        return(JaccardIndex)
    })
    
    res_pre <- do.call(rbind, res_pre)
    rownames(res_pre) <- unique(human.df$cluster)
    return(res_pre)
    
})
    
res1 <- do.call(cbind, res)

# #order the celltypes
# rowTarg <- c("Endothelial","Fibroblast" ,"Osteoblast","Cycling osteoblast",
#        "CD14_monocyte","NR4A3_Macrophage","TXNIP_Macrophage","FABP5_Macrophage",
#        "IFN-TAM","Pre-OC","Mature-OC","pDC","cDC2","Mast",
#        "CD8 T cell","CD4 T cell","IFN-sig T cell","NK cell",
#        "Plasma cell","B cell")
# colTarg <- c("Endothelial cell","Fibroblast","Osteoblast_1","Osteoblast_2","Osteoblast_3","Hypoxic_osteoblast",
#          "IFN-osteoblast","Osteoblast_cycling","CD4+_TIM","CD4-_TIM","ANGIO_TAM","TAM_INT","TAM_ACT","LA-TAM_SPP2_hi",
#          "LA-TAM_C1QC_hi","IFN-TAM","Cycling_OC","CD320_OC","Mature_OC","pDC","cDC1","cDC2","mregDC",
#          "Mast cell","Neutrophil","CD8 T cell","CD4 T cell","T_IFN","T_cycling","NK","Plasma cell","B cell")
# res1 <- res1[match(rowTarg, rownames(res1)),]        
# res1 <- res1[ ,match(colTarg, colnames(res1))]   
       
#plot the data
png(file = paste0("../output/", outName, "/", outName, "_jaccard.png"), width=4000, height=4000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(t(res1), #name = "mat", #col = col_fun,
              name = "Jaccard similarity index",
              cluster_rows = T,
              row_title = "Canine cell types",
              row_title_gp = gpar(fontsize = 24),
              col=viridis(option = "magma",100),
              cluster_columns = T,
              column_title = "Human cell types",
              column_title_gp = gpar(fontsize = 24),
              column_title_side = "bottom",
              column_names_rot = 45,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
             )
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),heatmap_legend_side = "top")
dev.off()




# ### Transfer labels over from 
# hu.reference <- readRDS("../output/s3/manton_ref_res0.8_dims45_dist0.2_neigh20_S3.rds")

# #transfer scArches_Cluster annotations
# ref.anchors <- FindTransferAnchors(reference = hu.reference, query = seu.obj, dims = 1:30,
#     reference.reduction = "pca", features = rownames(seu.obj)[rownames(seu.obj) %in% rownames(hu.reference)])
# predictions <- TransferData(anchorset = ref.anchors, refdata = hu.reference$mantonID, dims = 1:30)
# seu.obj <- AddMetaData(seu.obj, metadata = predictions)


# #plot the previously annotated labels (only for healthy cells)
# pi <- DimPlot(seu.obj, 
#               reduction = "umap.integrated.harmony", 
#               group.by = "predicted.id",
#               split.by = "predicted.id",
#               pt.size = 0.1,
#               ncol = 5,
#               repel = F
# )
# p <- formatUMAP(plot = pi) + NoLegend()
# ggsave(paste0("../output/", outName, "/", outName, "_mantonID_split_UMAP.png"), width = 14, height = 14)


# #plot transfered labels
# pi <- DimPlot(seu.obj, 
#               reduction = "umap.integrated.harmony", 
#               group.by = "predicted.id",
#               pt.size = 0.1,
#               label = T,
#               label.box = T,
#               repel = T
# )
# p <- formatUMAP(plot = pi) + NoLegend()
# ggsave(paste0("../output/", outName, "/", outName, "_mantonID_transfer_UMAP.png"), width = 7, height = 7)


# #transfer lsc_class annotations
# predictions <- TransferData(anchorset = ref.anchors, refdata = hu.reference$lsc_class, dims = 1:30)
# seu.obj <- AddMetaData(seu.obj, metadata = predictions)


# seu.obj$predicted.id <- factor(seu.obj$predicted.id)
# #inspect data for proper import -- looks good
# pi <- DimPlot(seu.obj, 
#               reduction = "umap.integrated", 
#               group.by = "predicted.id",
#               split.by = "predicted.id",
#               pt.size = 0.1,
#               ncol = 3,
#               label = F,
#               label.box = F,
#               repel = F
# ) + NoLegend()
# ggsave(paste0("../output/", outName, "/", outName, "_lsc_class_transfer_UMAP.png"), width = 12, height = 4)

#update user of completion
message(paste0(Sys.time(), " INFO: end of script!"))
