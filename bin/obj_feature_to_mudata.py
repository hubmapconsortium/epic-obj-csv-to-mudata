#!/usr/bin/env python3
import csv
import json
import re
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

import anndata
import mudata
import numpy as np
import pandas as pd

mudata.set_options(pull_on_update=False)

obj_file_pattern = re.compile(r"^(.+)-objects\.tsv$")
known_column_classes = {
    "mask": "mask",
    "spatial": "spatial",
}
# TODO: expand as more are known:
X_cols = {"antibody", "IMS", "gene", "default"}
known_col_sets = X_cols | set(known_column_classes.values())
type_mapping = {
    "integer": int,
    "float": float,
}
nullable_boolean_values = frozenset({"", False, True})


def subset_df_cols(
    data: pd.DataFrame,
    header: pd.DataFrame,
    col_type: str,
) -> tuple[pd.DataFrame, dict[str, str]]:
    df = data.loc[:, header.columns[header.loc["Feature class", :] == col_type]]
    new_columns = [c.casefold() for c in df.columns]
    col_mapping = dict(zip(new_columns, df.columns))
    df.columns = [c.casefold() for c in df.columns]
    return df, col_mapping


def df_cols_to_anndata(
    data: pd.DataFrame,
    header: pd.DataFrame,
    col_type: str,
) -> tuple[anndata.AnnData, dict[str, str]]:
    X_df, mapping = subset_df_cols(data, header, col_type)

    var = pd.DataFrame(
        {
            "alternative_id": header.loc["Alternative ID", X_df.columns],
            "description": header.loc["Description", X_df.columns],
            "protocol_doi": header.loc[
                "Protocol used to derive this feature value (DOI)", X_df.columns
            ],
        }
    )

    return anndata.AnnData(X=np.array(X_df), var=var), mapping


def check_duplicate_objects(data: pd.DataFrame):
    if len(set(data.index)) == data.shape[0]:
        return
    counts = data.index.value_counts()
    duplicates = counts[counts > 1]
    message_pieces = [
        "Found duplicate object IDs:",
        *(f"\t{i}\t({count} occurrences)" for i, count in duplicates.items()),
    ]
    raise ValueError("\n".join(message_pieces))


def get_header_row_count(tsv_path: Path) -> int:
    """
    :param tsv_path:
    :return: The value to use for `skiprows` in `pandas.read_table
    """
    with open(tsv_path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        for i, row in enumerate(r):
            if all(col in {"", "[ THIS ROW MUST BE PRESENT AND EMPTY ]"} for col in row):
                return i
    raise ValueError("Couldn't find empty row between header rows and content")


def read_tsv(tsv_path: Path) -> mudata.MuData:
    print("Reading", tsv_path)
    header_row_count = get_header_row_count(tsv_path)
    header = pd.read_table(tsv_path, nrows=header_row_count, index_col=0, header=None)
    data = pd.read_table(tsv_path, skiprows=header_row_count + 1, index_col=0)
    header.columns = data.columns
    header.loc["Feature class", :] = header.loc["Feature class", :].apply(str.casefold)
    orig_object_ids = list(data.index)
    # Use types in header in case Pandas is wrong, or data is malformed.
    # Coerce boolean to float, to allow NaNs if concatenating data frames
    # without the same set of columns.
    for column in header:
        if (column_type := header.loc["Type", column]) in type_mapping:
            data[column] = data[column].astype(type_mapping[column_type])
    check_duplicate_objects(data)
    data.index = data.index.astype(str)
    filename_piece = obj_file_pattern.match(tsv_path.name).group(1)
    data.index = [f"{filename_piece}-{i}" for i in data.index]
    data.index.name = "Object"
    other_col_types = set(header.loc["Feature class", :]) - known_col_sets

    # Special handling for known mask, data (X), and spatial keys:
    # mask data goes in overall .obs, spatial information goes in
    # the X_spatial key in .obsm
    obs, obs_col_mapping = subset_df_cols(data, header, known_column_classes["mask"])
    obs["object id"] = orig_object_ids
    orig_col_mapping = {"obs": obs_col_mapping, "mod": {}}
    spatial, spatial_col_mapping = subset_df_cols(data, header, known_column_classes["spatial"])
    orig_col_mapping["obsm"] = {"X_spatial": spatial_col_mapping}
    obsm = {"X_spatial": spatial}

    # Each feature class of primary measurement gets its own AnnData
    # (maybe empty here if there are no columns of a particular class)
    adatas = {}
    for col_class in X_cols:
        ad, mapping = df_cols_to_anndata(data, header, col_class)
        adatas[col_class] = ad
        if mapping:
            orig_col_mapping["mod"][col_class] = mapping
    adatas_to_use = {modality: adata for modality, adata in adatas.items() if adata.shape[1]}
    if not adatas_to_use:
        adatas_to_use["default"] = anndata.AnnData(
            X=None,
            shape=(obs.shape[0], 0),
            obs=pd.DataFrame(index=obs.index),
        )

    # Load everything else into its own key in .obsm (neighborhood, annotations)
    for remaining in other_col_types:
        df, col_mapping = subset_df_cols(data, header, remaining)
        obsm[remaining] = df
        orig_col_mapping["obsm"][remaining] = col_mapping

    mdata = mudata.MuData(
        adatas_to_use,
        obs=obs,
        obsm=obsm,
        uns={"column_orig_name_mapping": orig_col_mapping},
    )

    print(mdata)
    return mdata


def read_convert_tsv(input_dir: Path):
    tsvs = input_dir.glob("**/*-objects.tsv")
    mudatas = [read_tsv(tsv) for tsv in tsvs]
    obs = pd.concat([md.obs for md in mudatas])
    obsm_pieces = defaultdict(list)
    mod_pieces = defaultdict(list)
    for md in mudatas:
        for key, ad in md.mod.items():
            mod_pieces[key].append(ad)
        obsm_keys = set(md.obsm) - {"default"}
        for key in obsm_keys:
            item = md.obsm[key]
            obsm_pieces[key].append(item)
    obsm = {}
    for key, values in obsm_pieces.items():
        if all(isinstance(item, np.ndarray) for item in values):
            obsm[key] = np.vstack(values)
        else:
            df = pd.concat(values).copy()
            for col in df.columns:
                if df.dtypes[col] == "object":
                    df[col] = df[col].fillna("")
            obsm[key] = df

    mod = {}
    for key, pieces in mod_pieces.items():
        mod[key] = anndata.concat(pieces)

    uns = defaultdict(dict)
    for md in mudatas:
        for key, value in md.uns.items():
            uns[key] |= value

    mdata = mudata.MuData(mod, obs=obs, obsm=obsm, uns=dict(uns))
    for key, item in mdata.obsm.items():
        if isinstance(item, np.ndarray):
            continue
        for column in item.columns:
            if set(item.loc[:, column]) == nullable_boolean_values:
                # found nullable boolean
                new_arr = [item if isinstance(item, bool) else None for item in item[column]]
                item[column] = pd.array(new_arr)
    return mdata


def extract_metadata_write_json(mdata: mudata.MuData, output_json: Path):
    data = {}
    if "ontology" in mdata.obsm:
        if "object type" in mdata.obsm["ontology"]:
            data["object_types"] = sorted(set(mdata.obsm["ontology"]["object type"]))
    if "annotation tool" in mdata.obs:
        data["annotation_tools"] = sorted(set(mdata.obs["annotation tool"]))
    if "mask name" in mdata.obs:
        data["mask_names"] = sorted(set(mdata.obs["mask name"]))
    image_dimension = "2D"
    if "X_spatial" in mdata.obsm:
        # might be a DataFrame
        spatial_array = np.array(mdata.obsm["X_spatial"])
        if spatial_array.shape[1] >= 3:
            z_value_count = len(np.unique(spatial_array[:, 2]))
            if z_value_count > 1:
                image_dimension = "3D"
    data["image_dimension"] = image_dimension
    with open(output_json, "w") as f:
        json.dump(data, f)


def main(input_dir: Path, output_h5mu: Path, output_json: Path):
    mdata = read_convert_tsv(input_dir)
    print("Overall MuData:")
    print(mdata)

    mdata.write_h5mu(output_h5mu)
    extract_metadata_write_json(mdata, output_json)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("input_dir", type=Path)
    p.add_argument("output_h5mu", type=Path, nargs="?", default=Path("objects.h5mu"))
    p.add_argument("output_json", type=Path, nargs="?", default=Path("calculated_metadata.json"))
    args = p.parse_args()

    mdata = main(args.input_dir, args.output_h5mu, args.output_json)
