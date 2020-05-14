# Processing HuBMAP Cell Set Annotations

columns:
- cell name: The cell barcode label from each experiment
- annotation: Our predicted annotation. Each value is a term from the EBI Cell Ontology
- prediction score: Confidence level in prediction, ranging from 0..1
- dataset ID: Links to the dataset IDs from the Loaded Dataset Information HuBMAP google sheet. An example entry is HBM336.FWTN.636
- Globus ID: corresponds to the directory on Globus where the original input data was uploaded (for example, 2dca1bf5832a4102ba780e9e54f6c350)


## Setup

```sh
conda env create -f environment.yml
conda activate hubmap-cell-sets
```