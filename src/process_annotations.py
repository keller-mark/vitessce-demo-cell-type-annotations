# %%
import pandas as pd
import numpy as np
import json
import os
from os.path import join
import requests

# %%
# TODO: load data from real file

DATASET_ID_1 = "3683b49e27133c0.1"
DATASET_ID_2 = "3683b49e27133c0.2"

df = pd.DataFrame(data=[
    {
        "cell_id": "TTTGTTGAGCAGCGAT",
        "annotation": "B cell",
        "prediction_score": 0.9885534,
        "dataset_id": DATASET_ID_1
    },
    {
        "cell_id": "CGATGCGTCGATAACC",
        "annotation": "B cell",
        "prediction_score": 1.00000,
        "dataset_id": DATASET_ID_1
    },
    {
        "cell_id": "ATTGGGTTCACCGACG",
        "annotation": "B cell",
        "prediction_score": 1.00000,
        "dataset_id": DATASET_ID_1
    },
    {
        "cell_id": "TTGCGTCGTAGCGTCC",
        "annotation": "CD14-positive monocyte",
        "prediction_score": 0.8830666,
        "dataset_id": DATASET_ID_1
    },
    {
        "cell_id": "GCGATCGGTACGGGAT",
        "annotation": "B cell",
        "prediction_score": 0.9359248,
        "dataset_id": DATASET_ID_1
    },
    {
        "cell_id": "GAGATGGTCCTGCTAC",
        "annotation": "alpha-beta T cell",
        "prediction_score": 0.9934960,
        "dataset_id": DATASET_ID_1
    },
    {
        "cell_id": "TCGTGGGCAACGACAG",
        "annotation": "alpha-beta T cell",
        "prediction_score": 0.8950038,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "ATGGATCTCGGACAAG",
        "annotation": "alpha-beta T cell",
        "prediction_score": 0.9203139,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "ATGCGATGTACACGTT",
        "annotation": "B cell",
        "prediction_score": 1.000000,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "TCTCCGACACCAAATC",
        "annotation": "splenic macrophage",
        "prediction_score": 0.8988188,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "TTGTGTTAGTAAGAGG",
        "annotation": "B cell",
        "prediction_score": 1.00000,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "GTCTACCAGATTGAGT",
        "annotation": "B cell",
        "prediction_score": 1.00000,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "GATGATCCACTGCGAC",
        "annotation": "B cell",
        "prediction_score": 0.9589355,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "TTACCGCGTCGGCACT",
        "annotation": "alpha-beta T cell",
        "prediction_score": 0.9748914,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "ATGGTTGGTCATTCCC",
        "annotation": "B cell",
        "prediction_score": 1.0000,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "GTTCCTTCATGGCGCT",
        "annotation": "alpha-beta T cell",
        "prediction_score": 1.0000,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "GTTGTAGAGTTACTCG",
        "annotation": "B cell",
        "prediction_score": 0.9886951,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "GCAACATAGGAGTACC",
        "annotation": "B cell",
        "prediction_score": 1.0000,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "GCGATCGTCGGAACTT",
        "annotation": "B cell",
        "prediction_score": 0.9737594,
        "dataset_id": DATASET_ID_2
    },
    {
        "cell_id": "GAGGGTATCAGCGCGT",
        "annotation": "B cell",
        "prediction_score": 0.9890959,
        "dataset_id": DATASET_ID_2
    }
])
df

# %%

