from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import Any, Union

import anndata as ad # type: ignore
import h5py # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
from dask_image.imread import imread # type: ignore
from imageio import imread as imread2 # type: ignore
from xarray import DataArray # type: ignore
from scipy.sparse import coo_matrix
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel

from spatialdata_io.readers._utils._utils import _initialize_raster_models_kwargs # type: ignore
from typing import Any, Union, Optional, List

# Path keys for GEF file structure
GENE_EXP = "geneExp"
EXPRESSION = "expression"
# Column names
GENE = "gene"
GENE_ID = "geneID"
GENE_NAME = "geneName"
COUNT = "count"
COORD_X = "x"
COORD_Y = "y"
# Key words for spatial data
SPATIAL_KEY = "spatial"
REGION_KEY = "region"
INSTANCE_KEY = "instance_id"


def stereoseq(
    squarebin_gef: Optional[Union[str, Path]] = None,
    bins: Union[List[int], int] = [1, 20, 40, 50, 100, 150, 200],
    regist_image: Optional[Union[str, Path]] = None,
    cellbin_gef: Optional[Union[str, Path]] = None,
    tissue_mask: Optional[Union[str, Path]] = None,
    cell_mask: Optional[Union[str, Path]] = None,
    imread_kwargs: Mapping[str, Any] = {},
    image_models_kwargs: Mapping[str, Any] = {},
) -> SpatialData:
    """
    Read Stereo-seq data and convert it into SpatialData format.
    
    Args:
        squarebin_gef: Path to the square bin GEF file (tissue.gef for v8 or ssDNA.gef for v7)
        bins: List of bin sizes to process or a single integer bin size
        regist_image: Path to the registered image (HE or ssDNA image)
        cellbin_gef: Path to the cell bin GEF file (not yet supported)
        tissue_mask: Path to tissue mask (not yet supported)
        cell_mask: Path to cell mask (not yet supported)
        imread_kwargs: Keyword arguments for image reading
        image_models_kwargs: Keyword arguments for image model initialization
    
    Returns:
        SpatialData: A SpatialData object containing the processed data
    """
    
    # Initialize image model keyword arguments
    image_models_kwargs, _ = _initialize_raster_models_kwargs(image_models_kwargs, {})
    
    # Initialize data containers
    points = {}
    tables = {}
    images = {}
    
    # Process registered image if provided
    if regist_image:
        # Handle two types of images: HE and ssDNA
        data = imread(regist_image, **imread_kwargs)
        if len(data.shape) != 4:
            # Standard image format
            images[Path(regist_image).stem] = Image2DModel.parse(
                    imread(regist_image, **imread_kwargs), dims=("c", "y", "x"), **image_models_kwargs
                )
        else:
            # Special format requiring reshaping
            data = imread2(regist_image, **imread_kwargs).squeeze().transpose(2, 0, 1)
            image = DataArray(data, dims=("c", "y", "x"))
            images[Path(regist_image).stem] = Image2DModel.parse(image, **image_models_kwargs)
    
    # Process square bin GEF file if provided
    if squarebin_gef:
        # Create points model using SquareBin.gef file
        path_squarebin = squarebin_gef
        squarebin_gef = h5py.File(str(path_squarebin), "r")
        
        # Determine GEF version based on SAW version (v7 or v8)
        gef_version = squarebin_gef.attrs['version'][0]
        if gef_version == 2:
            FEATURE_KEY = GENE
            DF_GENE_COL = [FEATURE_KEY, COUNT]
        elif gef_version == 4:
            FEATURE_KEY = GENE_ID
            DF_GENE_COL = [FEATURE_KEY, GENE_NAME, COUNT]
        else:
            raise ValueError(
                f"GEF version can only be 2 or 4, which depends on the SAW version 7 or 8."
            )
        
        # Convert bins to string format for GEF file keys
        bins = [f"bin{i}" for i in bins]

        # Iterate through each bin in the GEF file
        for i in squarebin_gef[GENE_EXP].keys():
            if i not in bins: 
                continue
            
            # Get gene information
            arr = squarebin_gef[GENE_EXP][i][GENE][:]
            df_gene = pd.DataFrame(arr, columns=DF_GENE_COL)
            for j in DF_GENE_COL[:-1]:
                df_gene[j] = df_gene[j].str.decode("utf-8")

            # Create dataframe for points model
            arr = squarebin_gef[GENE_EXP][i][EXPRESSION][:]
            df_points = pd.DataFrame(arr, columns=[COORD_X, COORD_Y, COUNT])

            # Unroll gene names by count to create mapping between coordinate counts and gene names
            df_points[FEATURE_KEY] = [
                name
                for _, (name, *symbol, cell_count) in df_gene[DF_GENE_COL].iterrows()
                for _ in range(cell_count)
            ]
            df_points[FEATURE_KEY] = df_points[FEATURE_KEY].astype("category")

            # Calculate expression index
            points_coords = df_points[[COORD_X, COORD_Y]].copy()
            points_coords.drop_duplicates(inplace=True)
            points_coords.reset_index(inplace=True, drop=True)
            points_coords["bin_id"] = points_coords.index

            # Create names for points and table elements
            name_points_element = f"{i}_genes"
            name_table_element = f"{i}_table"
            
            # Create mapping from index to bin_id
            index_to_bin_id = pd.merge(
                df_points[[COORD_X, COORD_Y]],
                points_coords,
                on=[COORD_X, COORD_Y],
                how="left",
                validate="many_to_one",
            )

            # Create sparse expression matrix
            expression = coo_matrix(
                (
                    df_points[COUNT],
                    (index_to_bin_id.loc[df_points.index]["bin_id"].to_numpy(), df_points[FEATURE_KEY].cat.codes),
                ),
                shape=(len(points_coords), len(df_points[FEATURE_KEY].cat.categories)),
            ).tocsr().astype("int64")
            
            # Create spatial coordinates array
            obsm = dict()
            obsm[SPATIAL_KEY] = points_coords[["x","y"]].to_numpy()
            
            # Create points element (using points instead of shapes for performance)
            points_element = PointsModel.parse(
                obsm[SPATIAL_KEY],
                pd.DataFrame(index=points_coords['bin_id'])
            )

            # Create observation dataframe
            obs = pd.DataFrame({INSTANCE_KEY: points_coords["bin_id"], REGION_KEY: name_points_element})
            obs.index = points_coords.apply(lambda row: row[COORD_X] << 32 | row[COORD_Y], axis = 1)
            obs.index.name = None
            obs[REGION_KEY] = obs[REGION_KEY].astype("category")
            obs["total_counts"] = expression.sum(axis=1)

            # Add more gene information to variable dataframe
            df_gene = df_gene.set_index(FEATURE_KEY)
            df_gene.index.name = None
            df_gene = df_gene.loc[df_points[FEATURE_KEY].cat.categories, :]
            
            # Create AnnData object
            adata = ad.AnnData(expression, obs=obs, var=df_gene, obsm=obsm)

            # Create table model
            table = TableModel.parse(
                adata,
                region=name_points_element,
                region_key=REGION_KEY,
                instance_key=INSTANCE_KEY,
            )

            # Store elements in dictionaries
            tables[name_table_element] = table
            points[name_points_element] = points_element

    # Check for unsupported features
    if cellbin_gef:
        raise ValueError(
            f"Cell bin GEF processing is not supported yet."
        )
    
    # Create and return SpatialData object
    sdata = SpatialData(images=images, tables=tables, points=points)
    return sdata


def main(out_zarr, **kwargs):
    """Main function to convert Stereo-seq data to SpatialData and write to Zarr format."""
    # Convert comma-separated bins string to list
    kwargs["bins"] = kwargs["bins"].split(",")
    
    # Convert data to SpatialData format
    sto = stereoseq(**kwargs)
    
    # Write SpatialData to Zarr format
    sto.write(out_zarr)


if __name__ == "__main__":
    # Command-line interface
    import argparse
    parser = argparse.ArgumentParser(description="Convert Stereo-seq data to SpatialData format")
    parser.add_argument(
        "--squarebin_gef",
        help="The tissue cut GEF file: `.ssDNA.gef` for SAW v7 or `.tissue.gef` for SAW v8",
        required=True,
    )
    parser.add_argument(
        "--bins",
        help="The bin sizes to create (comma-separated list)",
        required=True,
    )
    parser.add_argument(
        "--regist_image",
        help="Path to the registered image",
    )
    parser.add_argument(
        "--out_zarr",
        help="Output directory path for Zarr storage",
    )
    args = parser.parse_args()
    main(**args.__dict__)