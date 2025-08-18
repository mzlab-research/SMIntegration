# Spatial Metabolomics and Transcriptomics Data Integration Tutorial
## Environment Configuration
Install SpatialData following the official documentation:
```bash
# Create and activate conda environment
conda create -n spatialdata python=3.10
conda activate spatialdata

# Install SpatialData and visualization tools
pip install spatialdata napari napari-spatialdata
```
## Data Preparation
### **1. Input Data Formats**
  - Metabolite feature matrix: ./spatialdata/data/metab_C02928A2_pos.txt 
  - Gene feature matrix: ./spatialdata/data/C02928A2.tissue.gef
  - ssDNA-stained image: ./spatialdata/data/C02928A2_ssDNA_regist.tif

### **2. Convert Metabolite Matrix to Zarr**
```bash
python spatialmetab.py \
  --metab_file "metab_C02928A2_pos.txt" \
  --out_zarr "metab_C02928A2_pos.zarr" \
  --resolution 50
```

### **3. Convert Gene Matrix to Zarr**
```bash
python stereoseq.py \
  --bins 100 \
  --out_zarr "trans_C02928A2.zarr" \
  --regist_image "C02928A2_ssDNA_regist.tif" \
  --squarebin_gef "C02928A2.tissue.gef"
```
## Spatial Registration of Metabolite Matrix
### **1. Registration Workflow**
```python
import spatialdata as sd
from napari_spatialdata import Interactive
from spatialdata.transformations import (
    get_transformation_between_landmarks,
    set_transformation,
    Identity,
    Sequence
)

# Load datasets
mt = sd.read_zarr("metab_C02928A2_pos.zarr")
st = sd.read_zarr("trans_C02928A2.zarr")

# Step 1: Initial landmark selection
Interactive([st, mt])
```
### **2. Napari Interface Operations**
Create Landmark Layers:
  - Add two new layers named st and mt
  - Use point tool to add corresponding landmarks (maintain identical order)
  - Press Shift+L to link landmarks to their respective data layers
  - Press Shift+E to save landmarks to Zarr files
  - Close window after completion
### **3. Apply Affine Transformations**
```python
# First transformation
affine1 = get_transformation_between_landmarks(
    references_coords=st["st"],
    moving_coords=mt["mt"]
)
set_transformation(
    mt["mz_point"],
    affine1,
    to_coordinate_system="demo2"
)

# Second transformation (refinement)
Interactive([st, mt])  # Add new landmarks (st2, mt2)

affine2 = get_transformation_between_landmarks(
    references_coords=st["st2"],
    moving_coords=mt["mt2"]
)
set_transformation(
    mt["mz_point"],
    Sequence([affine1, affine2]),
    to_coordinate_system="demo2"
)

# Save aligned data
mt.write("mz_aligned.zarr")
```
Interactive landmark selection in Napari (spatialdata/landmark_alignment.png)


## KNN Interpolation
```bash
python knn_interpolation.py \
  -s "trans_C02928A2.zarr" \
  -m "mz_aligned.zarr" \
  -o "knn_interpolation.zarr" \
  --stPoint "bin100_genes" \
  --mzPoint "mz_point" \
  --stTable "bin100_table" \
  --mzTable "mz_table" \
  --stCoord "global" \
  --mzCoord "demo2"
```

## Convert to RDS/TXT Format
### **Required metadata files**:
  - metab_identi.xls with columns: mz, Name, KEGG.ID
  - id_to_name.tsv with columns: geneid, symbol

```bash
Rscript spatialdata_to_rds.R \
  "./spatialdata/data/" \
  50 \ 
  "C:/Users/denghaoke/AppData/Local/Programs/Python/Python311"
```


```bash
Rscript rds_to_txt.R \
  "./spatialdata/data/"
```
