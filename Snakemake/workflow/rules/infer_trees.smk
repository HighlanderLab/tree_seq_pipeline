

# ---------------------------------------------------------------------------- #
#                                TREE                                          #
# ---------------------------------------------------------------------------- #

rule prepare_sample_file:
    input:
        vcf="Tsinfer/{chromosome}_ancestral.vcf.gz",
        meta=config['meta']
    conda:
        "tsinfer"
    output:
        "Tsinfer/samples/{chromosome}.samples"
    params:
        chrLength= lambda wildcards:
            config['chromosome_length'][int(wildcards.chromosome.split("_")[1])],
        ploidy=config['ploidy'],
        snakemakedir=config['snakemakedir']
    shell:
        "python Snakemake/scripts/PrepareTsinferSampleFile.py "
        "{input.vcf} {input.meta} {output} {params.ploidy} {params.chrLength}"

rule infer:
    input:
        "Tsinfer/samples/{chromosome}.samples"
    conda:  "tsinfer"
    output:
        "Tsinfer/trees/{chromosome}.trees"
    params:
        snakemakedir=config['snakemakedir']
    shell:
        "python Snakemake/scripts/InferTrees.py {input} {output}"
