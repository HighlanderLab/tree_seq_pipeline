import re

# ---------------------------------------------------------------------------- #
#                                TREE                                          #
# ---------------------------------------------------------------------------- #

rule prepare_sample_file:
    input:
        vcf="../Project/Tsinfer/{chromosome}_ancestral.vcf.gz",
        meta=config['meta']
    conda:
        "tsinfer"
    output:
        "../Project/Tsinfer/samples/{chromosome}.samples"
    params:
        chrLength= lambda wildcards:  config['chromosome_length'][int(re.findall(r'\d+', wildcards.chromosome)[0])],
        ploidy=config['ploidy']
    shell:
        "python ../scripts/PrepareTsinferSampleFile.py "
        "{input.vcf} {input.meta} {output} {params.ploidy} {params.chrLength}"

rule infer:
    input:
        "../Project/Tsinfer/samples/{chromosome}.samples"
    conda:  "tsinfer"
    output:
        "../Project/Tsinfer/trees/{chromosome}.trees"
    shell:
        "python ../scripts/InferTrees.py {input} {output}"
