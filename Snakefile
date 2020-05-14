from os.path import join

SRC_DIR = "src"
DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
INTER_DIR = join(DATA_DIR, "intermediate")
PROC_DIR = join(DATA_DIR, "processed")

rule all:
    input:
        join(RAW_DIR, "annotations_spleen_0510.csv")

checkpoint split:
    input:
        join(RAW_DIR, "annotations_{instance}.csv")
    output:
        directory(join(INTER_DIR, "{instance}"))
    script:
        join(SRC_DIR, "split_by_dataset_id.py")

def list_annotation_files_after_split(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    return expand(
        join(INTER_DIR, "{instance}", "{dataset_id}.csv"),
        instance=wildcards.instance,
        dataset_id=glob_wildcards(join(checkpoint_output, "{dataset_id}.csv")).dataset_id
    )
    
rule process_all:
    input:
        list_annotation_files_after_split

rule process_dataset:
    input:
        join(INTER_DIR, "{instance}", "{dataset_id}.csv")
    output:
        join(PROC_DIR, "{instance}", "{dataset_id}.json")
    script:
        join(SRC_DIR, "process_dataset.py")
