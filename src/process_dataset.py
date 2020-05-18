import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ
import json

from anndata import read_h5ad
import pyarrow as pa
from pandas import DataFrame



def arrow_to_json(arrow_file, cells_json_file, factors_json_file):
    df = pa.ipc.open_file(arrow_file).read_pandas()
    df_items = df.T.to_dict().items()

    id_to_umap = {
        k: {
            "mappings": {"UMAP": [v['umap_x'], v['umap_y']]},
            "factors": {"Leiden Cluster": str(int(v['leiden']))}
        }
        for (k,v) in df_items
    }
    pretty_json_umap = json.dumps(id_to_umap).replace('}},', '}},\n')
    with open(cells_json_file, 'w') as f:
        f.write(pretty_json_umap)

    leiden_clusters = sorted(df['leiden'].unique().astype('uint8'))
    id_to_factors = {
        'Leiden Cluster': {
            'map': [str(cluster) for cluster in leiden_clusters],
            'cells': { k: v['leiden'] for (k,v) in df_items }
        }
    }
    pretty_json_factors = json.dumps(id_to_factors).replace('}},', '}},\n')
    with open(factors_json_file, 'w') as f:
        f.write(pretty_json_factors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_arrow_file',
        required=True,
        help='directory containing h5ad files to read'
    )
    parser.add_argument('-oc', '--output_cells_file',
        required=True,
        help='directory where arrow files should be written'
    )
    parser.add_argument('-of', '--output_factors_file',
        required=True,
        help='directory where arrow files should be written'
    )
    args = parser.parse_args()
    arrow_to_json(
        args.input_arrow_file,
        args.output_cells_file,
        args.output_factors_file
    )
