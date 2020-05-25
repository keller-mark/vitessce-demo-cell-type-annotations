# Processing Cell Type Annotations

columns:
- `cell_id`: The cell barcode label from each experiment
- `annotation`: Our predicted annotation. Each value is a term from the EBI Cell Ontology
- `prediction_score`: Confidence level in prediction, ranging from 0..1
- `dataset_id`: Links to the dataset IDs from the Loaded Dataset Information HuBMAP google sheet. An example entry is HBM336.FWTN.636
- `globus_id`: corresponds to the directory on Globus where the original input data was uploaded


## Setup

```sh
conda env create -f environment.yml
conda activate cell-type-annotation-for-vitessce
```

## Run

```sh
snakemake --cores 2
```
