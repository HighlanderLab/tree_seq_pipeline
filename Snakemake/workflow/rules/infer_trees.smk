import re

# ---------------------------------------------------------------------------- #
#                                TREE                                          #
# ---------------------------------------------------------------------------- #

rule prepare_sample_file:
    input:
        vcf=f'{vcfdir}/{{chromosome}}_ancestral.vcf',
        meta=config['meta']
    output:
        "../Project/Tsinfer/samples/{chromosome}.samples"
    params:
        chrLength= lambda wildcards:  config['chromosome_length'][int(re.findall(r'\d+', wildcards.chromosome)[0])],
        ploidy=config['ploidy']
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    shell:
        "python scripts/PrepareTsinferSampleFile.py "
        "{input.vcf} {input.meta} {output} {params.ploidy} {params.chrLength}"

rule infer:
    input:
        "../Project/Tsinfer/samples/{chromosome}.samples"
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    output:
        "../Project/Tsinfer/trees/{chromosome}.trees"
    shell:
        "python scripts/InferTrees.py {input} {output}"
