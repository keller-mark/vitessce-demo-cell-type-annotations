import pandas as pd
import numpy as np
import json

from constants import *

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0])
    df = df.loc[df[COLUMNS.DATASET_ID.value] == snakemake.wildcards[COLUMNS.DATASET_ID.value]]

    annotation_values = df[COLUMNS.ANNOTATION.value].unique().tolist()

    df["annotation_index"] = df[COLUMNS.ANNOTATION.value].apply(lambda val: annotation_values.index(val))


    
    factors = {
        "Cell Type Annotations": {
            "map": annotation_values,
            "cells": dict(zip(df[COLUMNS.CELL_ID.value].tolist(), df["annotation_index"].tolist()))
        }
    }
    
    with open(snakemake.output[0], 'w') as f:
        json.dump(factors, f)