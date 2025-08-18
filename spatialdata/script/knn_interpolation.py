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
    def __init__(self, sdata, point, table, coord, resolution=None) -> Self:
        self.sdata = sd.read_zarr(sdata)
        self.point = self.sdata[point]
        self.table = self.sdata[table]
        self.coord = coord
        self.point_name = point
        self.table_name = table
        self._pre_porcess(resolution)

    def _pre_porcess(self, resolution):
        spatial_coord = self._get_spatial_coordinate()
        self.kdt = KDTree(spatial_coord, 2)
        self.resolution = resolution if resolution else self._guess_resolution()

    def _get_spatial_coordinate(self) -> pd.DataFrame:
        """
        从配准后的zarr中计算transfrom的坐标，用于构建网络
        不要依赖anndata中的spatail坐标，通过instance_id去匹配
        """
        trans = sd.get_centroids(
            self.point, coordinate_system=self.coord
        ).compute()
        return trans

    def _guess_resolution(self) -> float:
        kdt_radius = self.kdt.query(self.kdt.data, k=[2])[0].min()
        return kdt_radius


def _calculate_weight(x):
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
    为每个空转数据寻找k个代谢数据。
    距离更近的代谢点有更高的权重。
    配准偏差定义为代谢点的缩放比例
    一个空转点周围不应超过 9 个代谢点（否则应该合并而不是加权平均），实际个数受限于distence。
    简单起见，视为圆形状，因此重叠面积与角度无关，仅与距离有关
    """

    # 计算合适的邻居范围
    distence = st.resolution * (st.resolution / mz.resolution)
    neighbors_num = 9  #

    dd, ii = mz.kdt.query(st.kdt.data, k=neighbors_num, distance_upper_bound=distence)
    # 过滤未找到邻居的点
    st_index = [np.all(np.isinf(d)) == False for d in dd]
    mz_distence = dd[st_index]
    mz_distence[np.isinf(mz_distence)] = np.nan
    mz_weight = np.apply_along_axis(_calculate_weight, 1, mz_distence)

    mz_max = len(mz.kdt.data)
    mz_index = [i[i != mz_max] for i in ii[st_index]]

    return st_index, mz_index, mz_weight


def _imputation(index, weight, adata, row_num):
    metabolite_X = np.zeros([row_num, adata.X.shape[1]], dtype=np.int64)
    index_max = adata.obs[adata.uns["spatialdata_attrs"]["instance_key"]].max()
    for e, (i, w) in enumerate(zip(index, weight)):
        n = i[i <= index_max]
        w = w[np.isnan(w) == False]
        metabolite_X[e] = adata.X[n].T.multiply(w).sum(axis=1).round().astype(int).T
        adata.X[n].T.multiply(w).toarray()
    metabolite_X = coo_matrix(metabolite_X).tocsr().astype("int64")
    return metabolite_X


def create_spatialdata(
    st: Spatial,
    paired_points: pd.DataFrame,
    gene_table: ad.AnnData,
    mz_table: ad.AnnData,
) -> SpatialData:
    """
    创建一个新的spatialdata，仅包含配对点、图像和注释表。
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
    根据最近邻居和权重计算出新的X矩阵以及相关信息
    """
    region = "paired_points"
    paired_points = sd.get_centroids(
            st.point, coordinate_system=st.coord
        ).compute()[st_index]
    gene_table = st.table[st_index].copy()
    gene_table.X = gene_table.X.astype("int64")
    region_key = gene_table.uns["spatialdata_attrs"]["region_key"]
    gene_table.obs[region_key] = region
    gene_table.obs[region_key] = gene_table.obs[region_key].astype("category")
    gene_table.uns["spatialdata_attrs"]["region"] = region

    metabolite_X = _imputation(mz_index, mz_weight, mz.table, gene_table.X.shape[0])
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
    寻找代谢和空转数据spot的配对关系，返回一个合并后的sptaildata
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
