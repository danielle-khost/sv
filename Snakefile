import yaml
import glob, os, pathlib


#Config file
configfile: "config/config.yaml"
res_config = yaml.load(open("config/resources.yaml"),Loader=yaml.FullLoader)
samsheet="config/samples.tsv"

samdict = {}
with open(samsheet) as fin:
    for line in fin:
        line = line.rstrip()
        data = line.split("\t")
        samdict[data[0]] = {'bamfile' : data[1]}

samples = list(samdict.keys())

#Global variables from config
ref_org = config["ref_org"]
int_files = config["intermediate_files"]

#Functions
def get_bam(wildcards):
    bfile = samdict[wildcards.sample]['bamfile']
    return bfile


#Global rule
rule all:
    input:
        expand("results/{sample}_finalCall.vcf", sample=samples)


#Rules
rule snifflesCall:
    input:
        ref = config["reference"],
        bamfile = get_bam
    output:
        vcf = int_files + "{sample}_vs_" + ref_org + "_sniffles.vcf",
        snf = int_files +  "{sample}_vs_" + ref_org + "_sniffles.snf"
    conda:
        "workflow/envs/sniffles.yaml"
    resources:
        mem_mb = res_config["sniffles_call"]["mem_mb"],
        time = res_config["sniffles_call"]["time"]
    threads:
        res_config['sniffles_call']['threads']
    shell:
        "sniffles --threads {threads} --reference {input.ref} --input {input.bamfile} --vcf {output.vcf} --snf {output.snf}"

rule cutesvCall:
    input:
        ref = config["reference"],
        bamfile = get_bam
    output:
        vcf = int_files + "{sample}_vs_" + ref_org + "_cutesv.vcf"
    conda:
        "workflow/envs/cutesv.yaml"
    resources:
        mem_mb = res_config["cutesv_call"]["mem_mb"],
        time = res_config["cutesv_call"]["time"]
    threads:
        res_config['cutesv_call']['threads']
    params:
        tmpdir = int_files + "{sample}"
    shell:
        """
        mkdir -p {params.tmpdir}
        cuteSV -t {threads} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DE 0.3 {input.bamfile} {input.ref} {output.vcf} {params.tmpdir}
        """

rule svimCall:
    input:
        ref = config["reference"],
        bamfile = get_bam
    output:
        vcf = int_files + "{sample}_vs_" + ref_org + "_svim.vcf"
    conda:
        "workflow/envs/svim.yaml"
    resources:
        mem_mb = res_config["svim_call"]["mem_mb"],
        time = res_config["svim_call"]["time"]
    params:
        outdir = "workflow/{sample}_svim"
    shell:
        """
        svim alignment {params.outdir} {input.bamfile} {input.ref}
        mv {params.outdir}/variants.vcf {output.vcf}
        """

rule vcfMerge:
    input:
        snif = int_files + "{sample}_vs_" + ref_org + "_sniffles.vcf",
        cute = int_files + "{sample}_vs_" + ref_org + "_cutesv.vcf",
        svim = int_files + "{sample}_vs_" + ref_org + "_svim.vcf"
    output:
        int_files + "{sample}_merged.vcf"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * res_config["vcf_merge"]["mem_mb"],
        time = res_config["vcf_merge"]["time"]
    params:
        surv_path = config["survivor_path"],
        text = "{sample}_vcfs.txt"
    shell:
        """
        ls {input.snif} {input.cute} {input.svim} > {params.text}
        {params.surv_path} merge {params.text} 1000 3 1 1 0 50 {output}
        """

rule vcfFilter:
    input:
        int_files + "{sample}_merged.vcf"
    output:
        filt = int_files + "{sample}_merged_filtered_nogaps.vcf"
    resources:
        mem_mb = res_config["vcf_filter"]["mem_mb"],
        time = res_config["vcf_filter"]["time"]
    params:
        gaps = config["gap_bed"],
        surv_path = config["survivor_path"],
        inter = int_files + "{sample}_merged_filtered.vcf"
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        """
        grep -v '<TRA>' {input} | bcftools view -i 'SVLEN<100000' - > {params.inter}
        {params.surv_path} filter {params.inter} {params.gaps} 50 -1 0 -1 {output.filt}
        """

rule forceCall:
    input:
        ref = config["reference"],
        bamfile = get_bam,
        knownsv = int_files + "{sample}_merged_filtered_nogaps.vcf"
    output:
        vcf = temp("results/{sample}_finalCall.vcf")
    resources:
        mem_mb = res_config["force_call"]["mem_mb"],
        time = res_config["force_call"]["time"]
    threads:
        res_config['force_call']['threads']
    conda:
        "workflow/envs/sniffles.yaml"
    shell:
        "sniffles --threads {threads} --input {input.bamfile} --genotype-vcf {input.knownsv} --vcf {output.vcf}"


rule fixHeader:
    input:
        "results/{sample}_finalCall.vcf"
    output:
        vcf =  "results/{sample}_finalCall.sample.vcf"
    resources:
        mem_mb = res_config["vcf_filter"]["mem_mb"],
        time = res_config["vcf_filter"]["time"]
    conda:
        "workflow/envs/bcftools.yaml"
    params:
        rename = temp("{sample}_rename.txt"),
        trimvcf = temp("{sample}_trimmed.vcf")
    shell:
        """
        echo "{sample}" > {params.rename}
        cut -f1-10 {input} > {params.trimvcf}
        bcftools reheader -s {params.rename} {params.trimvcf} -o {output.vcf}
        """