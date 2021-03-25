# Snakemake pipeline for analyzing gwas summary statistic data using flash
#
#
# LICENSE: CC0. Do what you want with the code, but it has no guarantees.
#          https://creativecommons.org/share-your-work/public-domain/cc0/
#
#
# source activate cause_large
#
# ./run_snakemake.sh
#
# don't forget to update cluster.yaml


import pandas as pd
from snakemake.utils import validate

localrules: all, summ_to_cor
###### Load configuration file
configfile: "config_sub.yaml"
#validate(config, schema="schemas/config.schema.yaml")

ss = pd.read_csv(config["input"]["sum_stats"], na_filter=False)

kmax = str(config["analysis"]["kmax"])
seeds = config["analysis"]["seeds"]

min_nonmissing = config["analysis"]["min_nonmissing"]
norm_type = config["analysis"]["norm_type"]
method = config["analysis"]["method"]
p_thresh = config["analysis"]["flash_pthresh"]

data_dir = config["out"]["data_dir"] #where the data is
out_dir = config["out"]["output_dir"] #where results will go
formatted_gwas_dir = config["out"]["formatted_gwas_dir"]
data_prefix = config["out"]["data_prefix"]
results_prefix = config["out"]["results_prefix"]

est_L = config["analysis"]["est_L"]

n_sub = config["analysis"]["n_sub"]
seeds = config["analysis"]["seeds"]
tsvd_thresh = config["analysis"]["tsvd_thresh"]
tsvd_k = config["analysis"]["tsvd_k"]


if est_L:
    #inp = expand(out_dir + results_prefix + "estL_z_ff-1_{mnm}nomiss_{pt}_seed{s}_sub{ns}.tsvd_{thresh}.{chr}.RDS",chr = range(1, 23), mnm = min_nonmissing, pt = p_thresh, s = seeds, ns = n_sub, thresh = tsvd_thresh)
    #inp = expand(out_dir + results_prefix + "estL_snp_list_123nomiss_0.01.tsvd_{thresh}.{chrom}.RDS", chrom = range(1, 23), thresh = tsvd_thresh)
    inp = expand(out_dir + results_prefix + "estL_snp_list_123nomiss_0.01.tsvd_{thresh}.gwas.tsv", thresh = tsvd_thresh)
else:
    inp = expand(out_dir + results_prefix + "fit_z_ff-1_{mnm}nomiss_{pt}_seed{s}_sub{ns}.tsvd_{thresh}.RDS", mnm = min_nonmissing, pt = p_thresh, s = seeds, ns = n_sub, thresh = tsvd_thresh)

    #inp2 = expand(out_dir + results_prefix + "snp_list_123nomiss_0.01.tsvd_{thresh}.RDS",  thresh = tsvd_thresh)

rule all:
    input: inp



rule tsvd1:
    input: flash_fit = out_dir + results_prefix + "fit_{flash_prefix}.RDS"
    output:  out = out_dir + results_prefix + "fit_{flash_prefix}.tsvd_{thresh}.RDS"
    params: k = tsvd_k, 
    shell: 'Rscript R/7_run_tsvd.R {output.out} NA \
            {wildcards.thresh} {params.k} {input.flash_fit}'

rule tsvd2:
    input: NB = expand(data_dir + data_prefix + "normbeta_ldpruned_plink.z.{chrom}.RDS", chrom = range(1, 23)),
           snps = out_dir + "{snp_list_prefix}.RDS"
    output:  out = out_dir + results_prefix + "fit_{snp_list_prefix}.tsvd_{thresh}.RDS"
    params: k = tsvd_k 
    shell: 'Rscript R/7_run_tsvd.R {output.out} 0 \
            {wildcards.thresh} {params.k} NA {input.snps} {input.NB}'

rule est_L:
    input:  fit  = out_dir + results_prefix + "fit_{tsvd_prefix}.RDS",
            data =  data_dir + data_prefix + "normbeta.z.{chrom}.RDS",
            R = data_dir + data_prefix + "R_estimate.z.RDS"
    output: out = out_dir + results_prefix + "estL_{tsvd_prefix}.{chrom}.RDS"
    wildcard_constraints: chrom = "\d+"
    shell: 'Rscript R/6_est_L_svd.R  {input.fit} {input.data} {input.R} {output.out}'

rule write_gwas:
    input: estL = expand(out_dir + results_prefix + "estL_{{prefix}}.{chrom}.RDS", chrom=range(1, 23))
    output: out = out_dir + results_prefix + "estL_{prefix}.gwas.tsv"
    shell: 'Rscript R/8_format_estL.R  {output.out} {input.estL} '
