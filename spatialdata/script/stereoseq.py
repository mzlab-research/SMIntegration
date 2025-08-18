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
# gef path
GENE_EXP = "geneExp"
EXPRESSION = "expression"
# colnames
GENE = "gene"
GENE_ID = "geneID"
GENE_NAME = "geneName"
COUNT = "count"
COORD_X = "x"
COORD_Y = "y"
# key word
SPATIAL_KEY = "spatial"
REGION_KEY = "region"
INSTANCE_KEY = "instance_id"

##def stereoseq(
#    squarebin_gef: str | Path = None,
#    bins: list | int = [1,20,40,50,100,150,200],
#    regist_image: str | Path = None,
#    cellbin_gef: str | Path = None,
#    tissue_mask: str | Path = None,
#    cell_mask: str | Path = None,
#    imread_kwargs: Mapping[str, Any] = MappingProxyType({}),
#    image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
#) -> SpatialData:
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

    image_models_kwargs, _ = _initialize_raster_models_kwargs(image_models_kwargs, {})
    points = {}
    tables = {}
    images = {}
    if regist_image:# 处理HE和ssDNA两种图像类型
        data = imread(regist_image, **imread_kwargs)
        if len(data.shape) != 4:
            images[Path(regist_image).stem] = Image2DModel.parse(
                    imread(regist_image, **imread_kwargs), dims=("c", "y", "x"), **image_models_kwargs
                )
        else:
            data = imread2(regist_image, **imread_kwargs).squeeze().transpose(2, 0, 1)
            image = DataArray(data, dims=("c", "y", "x"))
            images[Path(regist_image).stem] = Image2DModel.parse(image, **image_models_kwargs)
    if squarebin_gef:
        # create points model using SquareBin.gef
        path_squarebin = squarebin_gef
        squarebin_gef = h5py.File(str(path_squarebin), "r")
        
        gef_version = squarebin_gef.attrs['version'][0]
        if gef_version == 2:
            FEATURE_KEY = GENE
            DF_GENE_COL = [FEATURE_KEY, COUNT]
        elif gef_version == 4:
            FEATURE_KEY = GENE_ID
            DF_GENE_COL = [FEATURE_KEY, GENE_NAME, COUNT]
        else:
            raise ValueError(
                f"gef version only can be 2 or 4, which depends on the SAW version 7 or 8."
            )
        bins = [f"bin{i}" for i in bins]

        for i in squarebin_gef[GENE_EXP].keys():
            if i not in bins: continue
            # get gene info
            arr = squarebin_gef[GENE_EXP][i][GENE][:]
            df_gene = pd.DataFrame(arr, columns=DF_GENE_COL)
            for j in DF_GENE_COL[:-1]:
                df_gene[j] = df_gene[j].str.decode("utf-8")

            # create df for points model
            arr = squarebin_gef[GENE_EXP][i][EXPRESSION][:]
            df_points = pd.DataFrame(arr, columns=[COORD_X, COORD_Y, COUNT])

            # unroll gene names by count such that there exists a mapping between coordinate counts and gene names
            df_points[FEATURE_KEY] = [
                name
                for _, (name, *symbol, cell_count) in df_gene[DF_GENE_COL].iterrows()
                for _ in range(cell_count)
            ]
            df_points[FEATURE_KEY] = df_points[FEATURE_KEY].astype("category")

            # cauculate expression index
            points_coords = df_points[[COORD_X, COORD_Y]].copy()
            points_coords.drop_duplicates(inplace=True)
            points_coords.reset_index(inplace=True, drop=True)
            points_coords["bin_id"] = points_coords.index

            name_points_element = f"{i}_genes"
            name_table_element = f"{i}_table"
            index_to_bin_id = pd.merge(
                df_points[[COORD_X, COORD_Y]],
                points_coords,
                on=[COORD_X, COORD_Y],
                how="left",
                validate="many_to_one",
            )

            expression = coo_matrix(
                (
                    df_points[COUNT],
                    (index_to_bin_id.loc[df_points.index]["bin_id"].to_numpy(), df_points[FEATURE_KEY].cat.codes),
                ),
                shape=(len(points_coords), len(df_points[FEATURE_KEY].cat.categories)),
            ).tocsr().astype("int64")
            obsm = dict()
            obsm[SPATIAL_KEY] = points_coords[["x","y"]].to_numpy()
            
            # may use this bin_id for point is better:
            # it would be more natural to use shapes, but for performance reasons we use points
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
            df_gene = df_gene.set_index(FEATURE_KEY)
            df_gene.index.name = None
            df_gene = df_gene.loc[df_points[FEATURE_KEY].cat.categories, :]
            adata = ad.AnnData(expression, obs=obs, var=df_gene,obsm=obsm)

            table = TableModel.parse(
                adata,
                region=name_points_element,
                region_key=REGION_KEY,
                instance_key=INSTANCE_KEY,
            )

            tables[name_table_element] = table
            points[name_points_element] = points_element

    if cellbin_gef:
        raise ValueError(
            f"Not supported yet."
        )
    sdata = SpatialData(images=images, tables=tables, points=points)
    return sdata

def main(out_zarr,**kwargs):
    kwargs["bins"] = kwargs["bins"].split(",")
    sto = stereoseq(**kwargs)
    sto.write(out_zarr)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "--squarebin_gef",
        help="the tissucut gef file, `.ssDNA.gef` for v7 or `.tissue.gef` for v8",
        required=True,
    )
    parser.add_argument(
        "--bins",
        help="the bins you want to create",
        required=True,
    )
    parser.add_argument(
        "--regist_image",
        help="regisited image",
    )
    parser.add_argument(
        "--out_zarr",
        help="empty directory path",
    )
    args = parser.parse_args()
    main(**args.__dict__)

