configfile: "../config/tsinfer.yaml"
workdir: config['workdir']

rule all:
    input:
        expand("Tsinfer/Chr{chromosome}.trees", chromosome = range(1, config['noChromosomes'] + 1))

rule prepare_sample_file:
    input:
        "Tsinfer/Chr{chromosome}_ancestral.vcf.gz"
    output:
        "Tsinfer/Chr{chromosome}.samples"
    script:
        "../scripts/PrepareTsinferSampleFile.py"

rule infer:
    input:
        "Tsinfer/Chr{chromosome}.samples"
    output:
        "Tsinfer/Chr{chromosome}.trees"
    script:
        "../scripts/InferTrees.py"
