
import pandas as pd
from snakemake.utils import validate


###### Load configuration file
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

ss = pd.read_csv(config["input"]["sum_stats"], na_filter=False)

data_dir = config["out"]["data_dir"] #where the data is
out_dir = config["out"]["output_dir"] #where results will go

kmax = str(config["analysis"]["kmax"])
seed = str(config["analysis"]["seed"])
norm_type = str(config["analysis"]["norm_type"])

prefix = config["out"]["prefix"]

method = config["analysis"]["method"]

rule all:
  input: expand("file_{m}.txt", m = method)
    
rule test:
  output: out = "file_{m}.txt"
  shell: "touch {output.out}"
