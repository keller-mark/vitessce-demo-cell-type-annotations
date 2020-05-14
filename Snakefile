from os.path import join

configfile: 'config.yml'

SRC_DIR = "src"
DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
PROC_DIR = join(DATA_DIR, "processed")

print(config.items())

rule all:
    input:
        [ join(PROC_DIR, "{}", "{}.json").format(instance, dataset_id)
            for instance, instance_datasets in config.items()
            for dataset_id in instance_datasets ]


rule process_dataset:
    input:
        join(RAW_DIR, "{instance}.csv")
    output:
        join(PROC_DIR, "{instance}", "{dataset_id}.json")
    script:
        join(SRC_DIR, "process_dataset.py")
