import argparse

import json
from anndata import read_h5ad
import pyarrow as pa
import pandas as pd

from constants import *
from utils import *

def generate_json_files(
    input_cells_arrow_file, input_annotations_csv_file,
    output_cells_json_file, output_factors_json_file, output_flat_cell_sets_json_file
):

    cells_df = pa.ipc.open_file(input_cells_arrow_file).read_pandas()
    annotation_df = pd.read_csv(input_annotations_csv_file).set_index(COLUMNS.CELL_ID.value)

    df = cells_df.join(annotation_df)
    df.index = df.index.rename(COLUMNS.CELL_ID.value)

    df['leiden'] = df['leiden'].apply(lambda i: f"Cluster {i}")
    df[COLUMNS.ANNOTATION.value] = df[COLUMNS.ANNOTATION.value].astype(str)

    cells_df_items = df.T.to_dict().items()

    # Generate .cells.json
    cells_json = {
        k: {
            "mappings": {"UMAP": [v['umap_x'], v['umap_y']]},
            "factors": {
                "Leiden Clustering": v['leiden'],
                "Cell Type Annotation": v[COLUMNS.ANNOTATION.value]
            }
        }
        for (k,v) in cells_df_items
    }
    with open(output_cells_json_file, 'w') as f:
        json.dump(cells_json, f)

    # Generate .factors.json
    leiden_clusters = sorted(df['leiden'].unique().tolist())
    cell_types = sorted(df[COLUMNS.ANNOTATION.value].unique().tolist())
    factors_json = {
        "Leiden Clustering": {
            "map": leiden_clusters,
            "cells": { k: leiden_clusters.index(v['leiden']) for (k,v) in cells_df_items }
        },
        "Cell Type Annotation": {
            "map": leiden_clusters,
            "cells": { k: cell_types.index(v[COLUMNS.ANNOTATION.value]) for (k,v) in cells_df_items }
        }
    }
    with open(output_factors_json_file, 'w') as f:
        json.dump(factors_json, f)
    

    # Generate .flat.cell_sets.json
    df = df.reset_index()
    tree = init_tree()

    leiden_clusters_children = []
    for cluster_name, cluster_df in df.groupby("leiden"):
        leiden_clusters_children.append({
            "name": cluster_name,
            "set": cluster_df[COLUMNS.CELL_ID.value].unique().tolist()
        })
    
    tree["tree"].append({
        "name": "Leiden Clustering",
        "children": leiden_clusters_children
    })

    cell_type_annotation_children = []
    for cell_type, cell_type_df in df.groupby(COLUMNS.ANNOTATION.value):
        cell_type_annotation_children.append({
            "name": cell_type,
            "set": cell_type_df[COLUMNS.CELL_ID.value].unique().tolist()
        })
    
    tree["tree"].append({
        "name": "Cell Type Annotations",
        "children": cell_type_annotation_children
    })
    
    with open(output_flat_cell_sets_json_file, 'w') as f:
        json.dump(tree, f)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-ic', '--input_cells_arrow_file',
        type=str,
        required=True,
        help='Input Arrow file'
    )
    parser.add_argument('-ia', '--input_annotations_csv_file',
        type=str,
        required=True,
        help='Input Arrow file'
    )
    parser.add_argument('-oc', '--output_cells_json_file',
        required=True,
        help='Output cells.json file'
    )
    parser.add_argument('-of', '--output_factors_json_file',
        required=True,
        help='Output factors.json file'
    )
    parser.add_argument('-ofcs', '--output_flat_cell_sets_json_file',
        required=True,
        help='Output flat.cell_sets.json file'
    )
    args = parser.parse_args()
    generate_json_files(
        args.input_cells_arrow_file,
        args.input_annotations_csv_file,
        args.output_cells_json_file,
        args.output_factors_json_file,
        args.output_flat_cell_sets_json_file
    )
