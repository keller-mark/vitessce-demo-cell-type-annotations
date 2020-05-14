# %%
import pandas as pd
import numpy as np
from os.path import join

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0])

    for dataset_id, dataset_df in df.groupby('dataset_id'):
        dataset_df.to_csv(join(snakemake.output[0], f"{dataset_id}.csv"))
