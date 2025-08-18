args <- commandArgs(T)
main_dir<- args[1]
resolution<-as.numeric(args[2])#50 
pythondir<-args[3]#"C:/Users/denghaoke/AppData/Local/Programs/Python/Python311"
binsize<-1
input=file.path(main_dir,"knn_interpolation.zarr")
input_id=file.path(main_dir,"id_to_name.tsv")
input_i=file.path(main_dir,"metab_identi.xls")

outdir=main_dir
mz_table="mz_table"
st_table="gene_table"


suppressMessages(library(reticulate))
use_python(pythondir)
sd <- import("spatialdata")
options(Seurat.object.assay.version = 'v3')
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(rjson))
suppressMessages(library(scales))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(readr))

generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE){
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }
  
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius / max(dim(x = image))
  
  return(new(Class = 'VisiumV1', 
             image = image, 
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                          fiducial = scale.factors$fiducial_diameter_fullres, 
                                          hires = scale.factors$tissue_hires_scalef, 
                                          lowres = scale.factors$tissue_lowres_scalef), 
             coordinates = tissue.positions, 
             spot.radius = spot.radius))
}

create_seurat <- function(mat,cell_coords,binsize){
  rownames(cell_coords) <- colnames(mat)
  colnames(cell_coords) <- c("x","y")
  seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'Stereo', assay = 'Spatial',names.delim = ':', meta.data = cell_coords)
  
  cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
  cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1
  tissue_positions_list <- data.frame(row.names = rownames(cell_coords),
                                      tissue = 1,
                                      row = cell_coords$y, col = cell_coords$x,
                                      imagerow = cell_coords$y, imagecol = cell_coords$x)
  
  tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))
  scalefactors_json <- toJSON(list(fiducial_diameter_fullres = binsize,
                                   tissue_hires_scalef = binsize,
                                   tissue_lowres_scalef = binsize))
  
  spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                    scale.factors = fromJSON(scalefactors_json), 
                                    tissue.positions = tissue_positions_list)
  
  spatialObj <- spatialObj[Cells(x = seurat_spatialObj)]
  DefaultAssay(spatialObj) <- 'Spatial'
  
  seurat_spatialObj[['slice1']] <- spatialObj
  return(seurat_spatialObj)
}



identi<-fread(input_i,sep ="\t", header = TRUE)
identi$mz<-as.character(identi$mz)
sdata <- sd$read_zarr(input)
mz_coords <- as.data.frame(sdata[mz_table]$obsm['spatial'])
str(mz_coords)
mz_coords[,1]<-(mz_coords[,1])/(resolution*2)
mz_coords[,2]<-(mz_coords[,2])/(resolution*2)
mz_mat <- as(t(sdata[mz_table]$X), "CsparseMatrix")
mz<- rownames(sdata[mz_table]$var)

rownames(mz_mat) <- mz
colnames(mz_mat) <- paste0("sample",":",mz_coords[,1],"_",mz_coords[,2])#
mz_mat <- mz_mat[rownames(mz_mat) %in% identi$mz, ]

mz_obj <- create_seurat(mz_mat,mz_coords,binsize)

saveRDS(mz_obj, file = file.path(outdir,"metab.rds"))

identi<-identi |>
dplyr::select(mz,Name,KEGG.ID) |>
dplyr::filter(mz %in% rownames(mz_mat))
write_delim(identi,file.path(outdir,"metabolite_annotation.txt"),delim = "\t")

gc()

genetype="geneid"
sdata <- sd$read_zarr(input)
st_coords <- as.data.frame(sdata[st_table]$obsm['spatial'])
st_coords[,1]<-(st_coords[,1])/(resolution*2)
st_coords[,2]<-(st_coords[,2])/(resolution*2)
st_mat <- as(t(sdata[st_table]$X), "CsparseMatrix")
gc()
geneid<-rownames(sdata[st_table]$var)
positionid<-paste0("sample",":",st_coords[,1],"_",st_coords[,2])#rownames(sdata[st_table]$obs)

rownames(st_mat) <- geneid
colnames(st_mat) <- positionid
#id to name
if(genetype=="geneid"){
  id_to_name<-read_delim(input_id)
  geneid <- rownames(st_mat)
  st_mat <- as.data.table(as.matrix(st_mat))
  gc()
  setDT(st_mat)
  set(st_mat, j = "geneid", value = geneid)
  st_mat <- merge(st_mat, id_to_name, by = "geneid", all.x = TRUE)
  gc()
  numeric_cols <- sapply(st_mat, is.numeric)
  st_mat <- st_mat[, lapply(.SD, sum, na.rm = TRUE), by = symbol, .SDcols = numeric_cols]
  gc()
  genename<-st_mat$symbol
  st_mat <- st_mat[, symbol := NULL]
  st_mat <- st_mat[, geneid := NULL]
  st_mat<-as.matrix(st_mat)
  st_mat<-as(st_mat, "CsparseMatrix")
  rownames(st_mat) <- genename
  colnames(st_mat) <- positionid
  gc()
}

st_obj <- create_seurat(st_mat,st_coords,binsize)

vst_rds=FindVariableFeatures(st_obj,selection.method = "vst", nfeatures = 10000)
gene_vector <- VariableFeatures(object = vst_rds)
st_obj <- subset(st_obj, features = gene_vector)
st_mat <- st_mat[rownames(st_mat) %in% gene_vector, ]

gc()
saveRDS(st_obj, file = file.path(outdir,"trans.rds"))




