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
    meta_group_cluster.txt
    """
    COORD_X = "x"
    COORD_Y = "y"
    MID_COUNT = "intensity"#"MIDCount"
    GENE_ID = "mz"#"geneID"
    SPATIAL_KEY = "spatial"
    INSTANCE_KEY = "instance_id"
    REGION_KEY = "region"

    name_points_element = "mz_point"
    name_table_element = "mz_table"

    df_points = pd.read_table(file,sep="\t")
    df_points = df_points[df_points.intensity != 0]#df_points[df_points.MIDCount != 0]
    df_points.reset_index(inplace=True, drop=True)
    df_points[GENE_ID] = df_points[GENE_ID].astype("category")
    df_points[[COORD_X]] = df_points[[COORD_X]] * 2 * resolution # convert to physical coordinate, 1 px = 0.5 μm
    df_points[[COORD_Y]] = df_points[[COORD_Y]] * 2 * resolution # convert to physical coordinate, 1 px = 0.5 μm

    points_coords = df_points[[COORD_X,COORD_Y]].copy()
    points_coords.drop_duplicates(inplace=True)
    points_coords.reset_index(inplace=True, drop=True)
    points_coords["bin_id"] = points_coords.index

    index_to_bin_id = pd.merge(
        df_points[[COORD_X, COORD_Y]],
        points_coords,
        on=[COORD_X, COORD_Y],
        how="left",
        validate="many_to_one",
    )
    
    expression = coo_matrix(
        (
            df_points[MID_COUNT],
            (index_to_bin_id.loc[df_points.index]["bin_id"].to_numpy(), df_points[GENE_ID].cat.codes),
        ),
        shape=(len(points_coords), len(df_points[GENE_ID].cat.categories)),
    ).tocsr().astype("int64")
    

    obsm = dict()
    obsm[SPATIAL_KEY] = points_coords[["x","y"]].to_numpy()
    
    points_element = PointsModel.parse(
        obsm[SPATIAL_KEY],
        pd.DataFrame(index=points_coords['bin_id'])
    )

    obs = pd.DataFrame({INSTANCE_KEY: points_coords["bin_id"], REGION_KEY: name_points_element})
    obs.index = points_coords.apply(lambda row: row[COORD_X] << 32 | row[COORD_Y], axis = 1)
    obs.index.name = None
    obs[REGION_KEY] = obs[REGION_KEY].astype("category")
    obs["total_counts"] = expression.sum(axis=1)

    # add more gene info to var
    df_gene = df_points[[GENE_ID,MID_COUNT]].groupby(GENE_ID,observed=True).sum()
    df_gene.index.name = None
    adata = ad.AnnData(expression, obs=obs, var=df_gene,obsm=obsm)

    table = TableModel.parse(
        adata,
        region=name_points_element,
        region_key=REGION_KEY,
        instance_key=INSTANCE_KEY,
    )
    tables = {name_table_element:table}
    points = {name_points_element:points_element}
    sdata = SpatialData(tables=tables, points=points)
    return sdata

def main(out_zarr,**kwargs):
    sto = spatialmetab(**kwargs)
    sto.write(out_zarr)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="")
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
        help="empty directory path",
    )
    args = parser.parse_args()
    main(**args.__dict__)

