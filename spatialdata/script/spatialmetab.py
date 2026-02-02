__author__ = "zhengshanqiao"
__version__ = "0.0.1"
import pandas as pd # type: ignore
import anndata as ad # type: ignore
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel
from scipy.sparse import coo_matrix
from spatialdata import SpatialData

# mz  x       y       intensity        
# 103.9562616     1       1       0             
# 103.9562616     1       2       0             
# 103.9562616     1       3       0            
# 103.9562616     1       4       0            
# 103.9562616     1       5       0            
# 103.9562616     1       6       0             
# 103.9562616     1       7       0       

def spatialmetab(file,resolution):
    """
    Process spatial metabolomics data from meta_group_cluster.txt format
    """
    COORD_X = "x"
    COORD_Y = "y"
    MID_COUNT = "intensity"  # "MIDCount" in original format
    GENE_ID = "mz"  # "geneID" in original format
    SPATIAL_KEY = "spatial"
    INSTANCE_KEY = "instance_id"
    REGION_KEY = "region"

    name_points_element = "mz_point"
    name_table_element = "mz_table"

    # Read and filter data
    df_points = pd.read_table(file,sep="\t")
    df_points = df_points[df_points.intensity != 0]  # Filter out zero intensity points
    df_points.reset_index(inplace=True, drop=True)
    df_points[GENE_ID] = df_points[GENE_ID].astype("category")
    
    # Convert pixel coordinates to physical coordinates (1 px = 0.5 μm)
    df_points[[COORD_X]] = df_points[[COORD_X]] * 2 * resolution
    df_points[[COORD_Y]] = df_points[[COORD_Y]] * 2 * resolution

    # Create unique coordinate bins
    points_coords = df_points[[COORD_X,COORD_Y]].copy()
    points_coords.drop_duplicates(inplace=True)
    points_coords.reset_index(inplace=True, drop=True)
    points_coords["bin_id"] = points_coords.index

    # Map each data point to its bin ID
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
            df_points[MID_COUNT],
            (index_to_bin_id.loc[df_points.index]["bin_id"].to_numpy(), df_points[GENE_ID].cat.codes),
        ),
        shape=(len(points_coords), len(df_points[GENE_ID].cat.categories)),
    ).tocsr().astype("int64")
    

    # Prepare spatial coordinates for AnnData
    obsm = dict()
    obsm[SPATIAL_KEY] = points_coords[["x","y"]].to_numpy()
    
    # Create SpatialData points element
    points_element = PointsModel.parse(
        obsm[SPATIAL_KEY],
        pd.DataFrame(index=points_coords['bin_id'])
    )

    # Create observation metadata
    obs = pd.DataFrame({INSTANCE_KEY: points_coords["bin_id"], REGION_KEY: name_points_element})
    # Create unique index by combining x and y coordinates
    obs.index = points_coords.apply(lambda row: row[COORD_X] << 32 | row[COORD_Y], axis = 1)
    obs.index.name = None
    obs[REGION_KEY] = obs[REGION_KEY].astype("category")
    obs["total_counts"] = expression.sum(axis=1)

    # Add more gene information to variable metadata
    df_gene = df_points[[GENE_ID,MID_COUNT]].groupby(GENE_ID,observed=True).sum()
    df_gene.index.name = None
    adata = ad.AnnData(expression, obs=obs, var=df_gene, obsm=obsm)

    # Create SpatialData table element
    table = TableModel.parse(
        adata,
        region=name_points_element,
        region_key=REGION_KEY,
        instance_key=INSTANCE_KEY,
    )
    
    # Assemble SpatialData object
    tables = {name_table_element:table}
    points = {name_points_element:points_element}
    sdata = SpatialData(tables=tables, points=points)
    return sdata

def main(out_zarr,**kwargs):
    """
    Main function to process data and write to Zarr format
    """
    sto = spatialmetab(**kwargs)
    sto.write(out_zarr)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert spatial metabolomics data to SpatialData format")
    parser.add_argument(
        "--metab_file",
        help="the gem format file of spatial metabolome",
        dest="file",
        required=True,
    )
    parser.add_argument(
        "--resolution",
        help="the resolution of spatial metabolome, micron(μm)",
        required=True,
        type=int
    )
    parser.add_argument(
        "--out_zarr",
        help="empty directory path for output Zarr file",
    )
    args = parser.parse_args()
    main(**args.__dict__)