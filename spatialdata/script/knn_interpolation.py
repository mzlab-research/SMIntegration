__author__ = "zhengshanqiao"
__version__ = "0.0.1"

import spatialdata as sd
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel
from spatialdata.transformations import set_transformation, Identity
import anndata as ad
from scipy.spatial import KDTree
from scipy.sparse import coo_matrix
import numpy as np
import pandas as pd
from math import sqrt
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from typing import Union
from typing_extensions import Self


class Spatial:
    """Spatial data handler for spatial transcriptomics and metabolomics data."""
    
    def __init__(self, sdata, point, table, coord, resolution=None) -> Self:
        """
        Initialize Spatial object.
        
        Args:
            sdata: Path to spatial data zarr file
            point: Name of point data in spatial data
            table: Name of table data in spatial data
            coord: Coordinate system to use
            resolution: Optional resolution value, will be estimated if not provided
        """
        self.sdata = sd.read_zarr(sdata)
        self.point = self.sdata[point]
        self.table = self.sdata[table]
        self.coord = coord
        self.point_name = point
        self.table_name = table
        self._pre_porcess(resolution)

    def _pre_porcess(self, resolution):
        """Preprocess spatial data and build KDTree for nearest neighbor search."""
        spatial_coord = self._get_spatial_coordinate()
        self.kdt = KDTree(spatial_coord, 2)  # 2-dimensional KDTree
        self.resolution = resolution if resolution else self._guess_resolution()

    def _get_spatial_coordinate(self) -> pd.DataFrame:
        """
        Calculate transformed coordinates from registered zarr data for network construction.
        
        Returns:
            DataFrame with spatial coordinates
        """
        trans = sd.get_centroids(
            self.point, coordinate_system=self.coord
        ).compute()
        return trans

    def _guess_resolution(self) -> float:
        """
        Estimate spatial resolution based on nearest neighbor distances.
        
        Returns:
            Estimated resolution as minimum distance to nearest neighbor
        """
        kdt_radius = self.kdt.query(self.kdt.data, k=[2])[0].min()
        return kdt_radius


def _calculate_weight(x):
    """
    Calculate weights for metabolite points based on distances.
    
    Args:
        x: Array of distances
        
    Returns:
        Array of normalized weights
    """
    if len(x[np.isnan(x) == False]) == 1:
        ratio = np.full(9, np.nan)
        ratio[0] = 1
    else:
        ratio = 1 - x / np.nansum(x)
        ratio = ratio / np.nansum(ratio)
    return ratio


def get_point_pair(
    st: Spatial, mz: Spatial
) -> Union[np.ndarray, np.ndarray, np.ndarray]:
    """
    Find k nearest metabolite data points for each spatial transcriptomics spot.
    
    Metabolite points closer to the spot have higher weights.
    Registration deviation is defined as the scaling ratio of metabolite points.
    A spatial transcriptomics spot should have at most 9 metabolite points 
    (otherwise should be merged rather than weighted average).
    
    Args:
        st: Spatial transcriptomics data handler
        mz: Metabolomics data handler
        
    Returns:
        Tuple of:
        - Boolean array indicating which spatial transcriptomics spots have neighbors
        - Array of indices for metabolite neighbors
        - Array of weights for metabolite neighbors
    """

    # Calculate appropriate neighbor search range
    distance = st.resolution * (st.resolution / mz.resolution)
    neighbors_num = 9  # Maximum number of neighbors to consider

    # Find nearest neighbors
    dd, ii = mz.kdt.query(st.kdt.data, k=neighbors_num, distance_upper_bound=distance)
    
    # Filter spots without neighbors
    st_index = [np.all(np.isinf(d)) == False for d in dd]
    mz_distance = dd[st_index]
    mz_distance[np.isinf(mz_distance)] = np.nan
    
    # Calculate weights based on distances
    mz_weight = np.apply_along_axis(_calculate_weight, 1, mz_distance)

    mz_max = len(mz.kdt.data)
    mz_index = [i[i != mz_max] for i in ii[st_index]]

    return st_index, mz_index, mz_weight


def _imputation(index, weight, adata, row_num):
    """
    Perform weighted imputation of metabolite data.
    
    Args:
        index: Array of metabolite indices for each spot
        weight: Array of weights for each metabolite
        adata: AnnData object with metabolite data
        row_num: Number of rows in output matrix
        
    Returns:
        Imputed metabolite expression matrix
    """
    metabolite_X = np.zeros([row_num, adata.X.shape[1]], dtype=np.int64)
    index_max = adata.obs[adata.uns["spatialdata_attrs"]["instance_key"]].max()
    
    for e, (i, w) in enumerate(zip(index, weight)):
        n = i[i <= index_max]
        w = w[np.isnan(w) == False]
        metabolite_X[e] = adata.X[n].T.multiply(w).sum(axis=1).round().astype(int).T
        
    metabolite_X = coo_matrix(metabolite_X).tocsr().astype("int64")
    return metabolite_X


def create_spatialdata(
    st: Spatial,
    paired_points: pd.DataFrame,
    gene_table: ad.AnnData,
    mz_table: ad.AnnData,
) -> SpatialData:
    """
    Create a new SpatialData object containing only paired points, images, and annotation tables.
    
    Args:
        st: Spatial transcriptomics data handler
        paired_points: DataFrame of paired spatial coordinates
        gene_table: Gene expression AnnData table
        mz_table: Metabolite expression AnnData table
        
    Returns:
        New SpatialData object
    """
    points = {
        "paired_points": PointsModel.parse(
            paired_points
        )
    }
    tables = {"gene_table": gene_table, "mz_table": mz_table}
    images = {i: j.copy() for i, j in st.sdata.images.items()}
    sdata = SpatialData(images=images, tables=tables, points=points)
    return sdata


def get_paired_element(
    st: Spatial,
    mz: Spatial,
    st_index: np.ndarray,
    mz_index: np.ndarray,
    mz_weight: np.ndarray,
) -> Union[pd.DataFrame, ad.AnnData, ad.AnnData]:
    """
    Calculate new X matrix and related information based on nearest neighbors and weights.
    
    Args:
        st: Spatial transcriptomics data handler
        mz: Metabolomics data handler
        st_index: Boolean array indicating which spots have neighbors
        mz_index: Array of metabolite neighbor indices
        mz_weight: Array of metabolite neighbor weights
        
    Returns:
        Tuple of:
        - DataFrame of paired point coordinates
        - Gene expression AnnData table
        - Metabolite expression AnnData table
    """
    region = "paired_points"
    
    # Get coordinates for paired points
    paired_points = sd.get_centroids(
            st.point, coordinate_system=st.coord
        ).compute()[st_index]
    
    # Prepare gene expression table
    gene_table = st.table[st_index].copy()
    gene_table.X = gene_table.X.astype("int64")
    region_key = gene_table.uns["spatialdata_attrs"]["region_key"]
    gene_table.obs[region_key] = region
    gene_table.obs[region_key] = gene_table.obs[region_key].astype("category")
    gene_table.uns["spatialdata_attrs"]["region"] = region

    # Perform metabolite imputation
    metabolite_X = _imputation(mz_index, mz_weight, mz.table, gene_table.X.shape[0])
    
    # Create metabolite table
    var = pd.DataFrame(
        metabolite_X.sum(axis=0).T,
        index=mz.table.var.index.values,
        columns=["MIDCount"],
    )
    obsm = gene_table.obsm.copy()
    uns = gene_table.uns.copy()
    obs = gene_table.obs.copy()
    obs["total_counts"] = metabolite_X.sum(axis=1)

    mz_table = ad.AnnData(metabolite_X, obs=obs, var=var, obsm=obsm, uns=uns)
    
    return paired_points, gene_table, mz_table


def main(
    spatial, metab, st_point, st_table, st_coord, mz_point, mz_table, mz_coord, outdir
) -> None:
    """
    Main function to find pairing relationships between metabolomics and spatial transcriptomics spots.
    
    Returns a merged spatialdata object.
    
    Args:
        spatial: Path to spatial transcriptomics data
        metab: Path to metabolomics data
        st_point: Name of point data in spatial transcriptomics
        st_table: Name of table data in spatial transcriptomics
        st_coord: Coordinate system for spatial transcriptomics
        mz_point: Name of point data in metabolomics
        mz_table: Name of table data in metabolomics
        mz_coord: Coordinate system for metabolomics
        outdir: Output directory for results
    """
    st = Spatial(spatial, st_point, st_table, st_coord)
    mz = Spatial(metab, mz_point, mz_table, mz_coord)
    st_index, mz_index, mz_weight = get_point_pair(st, mz)
    paired_points, gene_table, mz_table = get_paired_element(
        st, mz, st_index, mz_index, mz_weight
    )
    sdata = create_spatialdata(st, paired_points, gene_table, mz_table)
    sdata.write(outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--spatial", required=True)
    parser.add_argument("-m", "--metab", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--stPoint", dest="st_point", required=True)  # point
    parser.add_argument("--mzPoint", dest="mz_point", required=True)  # point
    parser.add_argument("--stTable", dest="st_table", required=True)  # point
    parser.add_argument("--mzTable", dest="mz_table", required=True)  # point
    parser.add_argument("--stCoord", dest="st_coord", required=True)  # point
    parser.add_argument("--mzCoord", dest="mz_coord", required=True)  # point
    args = parser.parse_args()
    main(**args.__dict__)