import re

rule prepare_sample_file:
    input:
        vcf=f'{vcfdir}/{{chromosome}}_ancestral.vcf.gz',
        meta=config['meta']
    output:
        f"../{Project}/Tsinfer/samples/{{chromosome}}.samples"
    params:
        chrLength= lambda wildcards:  config['chromosome_length'][int(re.findall(r'\d+', wildcards.chromosome)[0])],
        ploidy=config['ploidy']
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/prepare_sample_file_{chromosome}.log'
    shell:
        "python scripts/PrepareTsinferSampleFile.py "
        "{input.vcf} {input.meta} {output} {params.ploidy} {params.chrLength}"

rule infer:
    input:
        f"../{Project}/Tsinfer/samples/{{chromosome}}.samples"
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/infer_{chromosome}.log'
    output:
        f"../{Project}/Tsinfer/trees/{{chromosome}}.trees"
    shell:
        "python scripts/InferTrees.py {input} {output}"
