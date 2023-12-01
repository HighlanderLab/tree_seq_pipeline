import re

# if multiallelic:
#     vcf_4_inference = f'{vcfOut}/{{chromosome}}_allbi.vcf.gz'
# else:
#     vcf_4_inference = f'{vcfOut}/{{chromosome}}_ancestral.vcf.gz'

rule prepare_sample_file:
    input:
        #vcf=vcf_4_inference,
        vcf = f'{vcfOut}/{{chromosome}}_ancestral.vcf.gz',
        meta = Path(vcfdir, config['meta'])
    output:
        f"{oDir}/Tsinfer/samples/{{chromosome}}.samples"
    params:
        chrLength= lambda wildcards:  config['chromosome_length'][int(re.findall(r'\d+', wildcards.chromosome)[0])],
        ploidy=config['ploidy']
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=128000, time_min=200
    log: 'logs/Prepare_sample_file_{chromosome}.log'
    shell:
        "python scripts/Prepare_tsinfer_sample_file.py "
        "{input.vcf} {input.meta} {output} {params.ploidy} {params.chrLength}"

rule infer:
    input:
        f"{oDir}/Tsinfer/samples/{{chromosome}}.samples"
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=128000, time_min=300
    log: 'logs/Infer_{chromosome}.log'
    output:
        f"{oDir}/Tsinfer/trees/{{chromosome}}.trees"
    script:
        "../scripts/Infer_trees.py" # path to script is relative to the file containing this rule. Hence the prepended "../"
