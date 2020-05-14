import pandas as pd
import numpy as np
import json
import networkx
import obonet

from constants import *
from utils import *

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0])
    df = df.loc[df[COLUMNS.DATASET_ID.value] == snakemake.wildcards[COLUMNS.DATASET_ID.value]]
    
    
    tree = {
        "datatype": "cell",
        "version": "0.1.2",
        "tree": [
            {
                "name": "Cell Type Annotations",
                "children": [
                    
                ]
            }
        ]
    }

    for cell_type, cell_type_df in df.groupby(COLUMNS.ANNOTATION.value):
        tree["tree"][0]["children"].append({
            "name": cell_type,
            "set": cell_type_df[COLUMNS.CELL_ID.value].unique().tolist()
        })
    
    with open(snakemake.output[0], 'w') as f:
        json.dump(tree, f)