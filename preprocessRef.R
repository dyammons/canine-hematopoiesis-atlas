#!/usr/bin/Rscript

#load custom functions & packages
message(paste0(Sys.time(), " INFO: loading packages and custom functions."))
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
outName <- "human_k9_comp"

#get file names
message(paste0(Sys.time(), " INFO: initiating script."))
files <- list.files(path = "../external_data/manton/" , pattern = "*.h5", all.files = FALSE,
                    full.names = F)

#load in the barcodes that have annotations to stash cell identities and prefilter
cellMeta.df <- read.csv("./metaData/manton_cell_IDs.csv", header = T)
cellMeta.df <- cellMeta.df %>% mutate(cleanBar = str_split(barcode, "_", simplify = T)[,4],
                                      samName = substr(barcode, 1, (nchar(barcode)-19))
                                     )

#get names so orig.ident can be stored in metadata
sampleNames <-  unlist(lapply(files, function(x){paste(strsplit(x, "_")[[1]][1:3], collapse = "_")}))

if(length(files[!sampleNames %in% cellMeta.df$samName]) > 0){
    message("WARN: The following samples do not have associated metadata and will be skipped: \n", paste(files[!sampleNames %in% cellMeta.df$samName], collapse = "\n"))
    files <- files[sampleNames %in% cellMeta.df$samName]
    sampleNames <- sampleNames[sampleNames %in% cellMeta.df$samName]
}


#load each in as a seperate seurat object
message(paste0(Sys.time(), " INFO: loading in matricies and creating Seurat object list."))
seu.list <- lapply(files, function(inFile){
 
    #get key paths
    inFile_pwd <- paste0("../external_data/manton/", inFile)
    projName <- paste(strsplit(inFile, "_")[[1]][1:3], collapse = "_")

    #load the data
    spar.matrix <- Read10X_h5(inFile_pwd, use.names = TRUE, unique.features = TRUE)
    seu.obj <- CreateSeuratObject(spar.matrix, project = projName)

    #load cell metadata
    cellMeta.df.sub <- cellMeta.df[cellMeta.df$samName == projName, ]
    cellNames <- cellMeta.df.sub$celltype
    names(cellNames) <- cellMeta.df.sub$cleanBar
    seu.obj <- AddMetaData(seu.obj, cellNames, col.name = "mantonID")

    #remove any cells not in the list
    seu.obj <- subset(seu.obj, cells = names(cellNames))
    message(paste0(Sys.time(), " INFO: Successfully loaded: ", inFile))
    return(seu.obj)
})

#merge all objects into one
message(paste0(Sys.time(), " INFO: samples loaded. Merging objects."))
seu.obj <- merge(seu.list[1][[1]], y = seu.list[2:length(seu.list)],
                  add.cell.ids = sampleNames, 
                  project = "manton_bm_ref"
                 )

#integrate the data using harmony
message(paste0(Sys.time(), " INFO: files loaded in as Seurat objects and integration is begining."))
seu.obj <- integrateData(outName = "manton_ref", orig.reduction = "pca", saveRDS = F,
                          normalization.method = "LogNormalize", method = "HarmonyIntegration",
                          indReClus = TRUE, seu.obj = seu.obj,
                          runAllMethods = FALSE
                        )

#complete data visualization & save the RDS file
message(paste0(Sys.time(), " INFO: data integration complete. generating some figures and saving integrated object as a .rds file in ../s3/."))
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "manton_ref", 
                        final.dims = 45, final.res = 0.8, stashID = "clusterID", algorithm = 3, min.dist = 0.2, n.neighbors = 20,
                        prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated",
                        saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "HLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)

#plot the umap with cell annotations
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated", 
              group.by = "mantonID",
              pt.size = 0.1,
              label = T,
              label.box = T,
              repel = F
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_minorIdent_UMAP.png"), width = 7, height = 7)

#update user of completion
message(paste0(Sys.time(), " INFO: end of script!"))
