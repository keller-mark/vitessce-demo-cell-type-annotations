from enum import Enum

# URLs
CL_OBO_URL = "https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl.obo"
CL_ROOT_ID = 'CL:0000000'

# Column name constants
class COLUMNS(Enum):
    CELL_ID = "cell_id"
    ANNOTATION = "annotation"
    PREDICTION_SCORE = "prediction_score"
    DATASET_ID = "dataset_id"
    GLOBUS_ID = "globus_id"
