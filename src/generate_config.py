import pandas as pd
import numpy as np
import json
import yaml

from constants import COLUMNS

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0])
    dataset_ids = df[COLUMNS.DATASET_ID.value].unique().tolist()
    # TODO
    config = { "[instance_here]": dataset_ids }
    with open(snakemake.output[0], 'w') as f:
        yaml.dump(config, f, default_flow_style=False)