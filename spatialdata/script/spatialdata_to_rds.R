# Parse command line arguments
args <- commandArgs(T)
main_dir <- args[1]  # Main directory path
resolution <- as.numeric(args[2])  # Resolution parameter (e.g., 50)
pythondir <- args[3]  # Python directory path
binsize <- 1  # Fixed bin size parameter

# Define input file paths
input <- file.path(main_dir, "knn_interpolation.zarr")
input_id <- file.path(main_dir, "id_to_name.tsv")
input_i <- file.path(main_dir, "metab_identi.xls")

# Define output directory (same as main_dir)
outdir <- main_dir
mz_table <- "mz_table"  # Metabolite table name in Zarr
st_table <- "gene_table"  # Gene table name in Zarr

# Load required R packages
suppressMessages(library(reticulate))
use_python(pythondir)  # Set Python environment
sd <- import("spatialdata")  # Import Python spatialdata module

# Set Seurat assay version and load R packages
options(Seurat.object.assay.version = 'v3')
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(rjson))
suppressMessages(library(scales))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(readr))

# Function to generate a spatial object compatible with Seurat
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
  # Filter matrix to include only tissue positions where tissue == 1
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }
  
  # Calculate spot radius for visualization
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius / max(dim(x = image))
  
  # Create and return VisiumV1 spatial object
  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      lowres = scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}

# Function to create Seurat object from spatial data
create_seurat <- function(mat, cell_coords, binsize) {
  # Set row names of coordinates to match matrix column names
  rownames(cell_coords) <- colnames(mat)
  colnames(cell_coords) <- c("x", "y")
  
  # Create Seurat object with spatial coordinates as metadata
  seurat_spatialObj <- CreateSeuratObject(
    counts = mat,
    project = 'Stereo',
    assay = 'Spatial',
    names.delim = ':',
    meta.data = cell_coords
  )
  
  # Normalize coordinates to start from 1
  cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
  cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1
  
  # Create tissue positions data frame for spatial object
  tissue_positions_list <- data.frame(
    row.names = rownames(cell_coords),
    tissue = 1,
    row = cell_coords$y,
    col = cell_coords$x,
    imagerow = cell_coords$y,
    imagecol = cell_coords$x
  )
  
  # Create low-resolution image matrix
  tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))
  
  # Create JSON string for scale factors
  scalefactors_json <- toJSON(list(
    fiducial_diameter_fullres = binsize,
    tissue_hires_scalef = binsize,
    tissue_lowres_scalef = binsize
  ))
  
  # Generate spatial object
  spatialObj <- generate_spatialObj(
    image = tissue_lowres_image,
    scale.factors = fromJSON(scalefactors_json),
    tissue.positions = tissue_positions_list
  )
  
  # Subset spatial object to match cells in Seurat object
  spatialObj <- spatialObj[Cells(x = seurat_spatialObj)]
  DefaultAssay(spatialObj) <- 'Spatial'
  
  # Add spatial object to Seurat object as a new slice
  seurat_spatialObj[['slice1']] <- spatialObj
  return(seurat_spatialObj)
}

# Read metabolite identification data
identi <- fread(input_i, sep = "\t", header = TRUE)
identi$mz <- as.character(identi$mz)  # Convert mz column to character

# Read metabolite data from Zarr file
sdata <- sd$read_zarr(input)

# Extract metabolite spatial coordinates and adjust by resolution
mz_coords <- as.data.frame(sdata[mz_table]$obsm['spatial'])
mz_coords[, 1] <- (mz_coords[, 1]) / (resolution * 2)
mz_coords[, 2] <- (mz_coords[, 2]) / (resolution * 2)

# Extract metabolite matrix and convert to sparse matrix format
mz_mat <- as(t(sdata[mz_table]$X), "CsparseMatrix")
mz <- rownames(sdata[mz_table]$var)  # Get metabolite IDs

# Set row and column names for metabolite matrix
rownames(mz_mat) <- mz
colnames(mz_mat) <- paste0("sample", ":", mz_coords[, 1], "_", mz_coords[, 2])

# Filter metabolite matrix to include only identified metabolites
mz_mat <- mz_mat[rownames(mz_mat) %in% identi$mz, ]

# Create Seurat object for metabolite data
mz_obj <- create_seurat(mz_mat, mz_coords, binsize)

# Save metabolite Seurat object as RDS file
saveRDS(mz_obj, file = file.path(outdir, "metab.rds"))

# Create metabolite annotation file
identi <- identi |>
  dplyr::select(mz, Name, KEGG.ID) |>
  dplyr::filter(mz %in% rownames(mz_mat))
write_delim(identi, file.path(outdir, "metabolite_annotation.txt"), delim = "\t")

# Garbage collection to free memory
gc()

# Process gene expression data
genetype <- "geneid"  # Gene identifier type

# Read spatial data again for gene expression
sdata <- sd$read_zarr(input)

# Extract gene spatial coordinates and adjust by resolution
st_coords <- as.data.frame(sdata[st_table]$obsm['spatial'])
st_coords[, 1] <- (st_coords[, 1]) / (resolution * 2)
st_coords[, 2] <- (st_coords[, 2]) / (resolution * 2)

# Extract gene expression matrix and convert to sparse matrix
st_mat <- as(t(sdata[st_table]$X), "CsparseMatrix")
gc()  # Garbage collection

# Get gene IDs and position IDs
geneid <- rownames(sdata[st_table]$var)
positionid <- paste0("sample", ":", st_coords[, 1], "_", st_coords[, 2])

# Set row and column names for gene expression matrix
rownames(st_mat) <- geneid
colnames(st_mat) <- positionid

# Map gene IDs to gene names if needed
if (genetype == "geneid") {
  id_to_name <- read_delim(input_id)  # Read gene ID to name mapping
  geneid <- rownames(st_mat)
  
  # Convert to data.table for efficient merging and aggregation
  st_mat <- as.data.table(as.matrix(st_mat))
  gc()
  setDT(st_mat)
  set(st_mat, j = "geneid", value = geneid)
  
  # Merge with gene name mapping
  st_mat <- merge(st_mat, id_to_name, by = "geneid", all.x = TRUE)
  gc()
  
  # Sum expression values by gene symbol
  numeric_cols <- sapply(st_mat, is.numeric)
  st_mat <- st_mat[, lapply(.SD, sum, na.rm = TRUE), by = symbol, .SDcols = numeric_cols]
  gc()
  
  # Extract gene symbols and remove metadata columns
  genename <- st_mat$symbol
  st_mat <- st_mat[, symbol := NULL]
  st_mat <- st_mat[, geneid := NULL]
  
  # Convert back to sparse matrix
  st_mat <- as.matrix(st_mat)
  st_mat <- as(st_mat, "CsparseMatrix")
  rownames(st_mat) <- genename
  colnames(st_mat) <- positionid
  gc()
}

# Create Seurat object for gene expression data
st_obj <- create_seurat(st_mat, st_coords, binsize)

# Find variable features for gene expression data
vst_rds <- FindVariableFeatures(st_obj, selection.method = "vst", nfeatures = 10000)
gene_vector <- VariableFeatures(object = vst_rds)

# Subset gene expression matrix to include only variable features
st_obj <- subset(st_obj, features = gene_vector)
st_mat <- st_mat[rownames(st_mat) %in% gene_vector, ]

# Garbage collection
gc()

# Save gene expression Seurat object as RDS file
saveRDS(st_obj, file = file.path(outdir, "trans.rds"))