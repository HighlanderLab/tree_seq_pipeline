configfile: "../config/tsinfer.yaml"
workdir: config['workdir']

rule all:
    input:
        expand("Tsinfer/Chr{chromosome}.trees", chromosome = range(1, config['noChromosomes'] + 1))

rule prepare_sample_file:
    input:
        vcf="Tsinfer/Chr{chromosome}_ancestral.vcf.gz",
        meta=config['meta']
    conda:
        "tsinfer"
    output:
        "Tsinfer/Chr{chromosome}.samples"
    params:
        chrLength= lambda wildcards: config['chromosome_length'][int(wildcards.chromosome)],
        ploidy=config['ploidy'],
        snakemakedir=config['snakemakedir']
    shell:
        "python {params.snakemakedir}/scripts/PrepareTsinferSampleFile.py {wildcards.chromosome} {input.vcf} {input.meta} {output} {params.ploidy} {params.chrLength}"

rule infer:
    input:
        "Tsinfer/Chr{chromosome}.samples"
    conda:  "tsinfer"
    output:
        "Tsinfer/Chr{chromosome}.trees"
    params:
        snakemakedir=config['snakemakedir']
    shell:
        "python {params.snakemakedir}/scripts/InferTrees.py {input} {output}"
