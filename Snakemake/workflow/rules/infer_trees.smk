import re

rule prepare_sample_file:
    input:
        vcf = f'{vcfdir}/{{chromosome}}_ancestral.vcf.gz',
        meta = config['meta']
    output:
        f'{projdir}/Tsinfer/samples/{{chromosome}}.samples'
    params:
        chrLength= lambda wildcards:  config['chromosome_length'][int(re.findall(r'\d+', wildcards.chromosome)[0])],
        ploidy=config['ploidy']
    conda: "HLab_tsinfer"
    threads: 32
    resources: cpus=32, mem_mb=768000, time_min=200
    benchmark:
            'benchmarks/{chromosome}.prep_sample.benchmark.txt'
    #log: 'logs/Prepare_sample_file_{chromosome}.log'
    shell:
        "python scripts/Prepare_tsinfer_sample_file.py "
        "{input.vcf} {input.meta} {output} {params.ploidy} {params.chrLength}"

rule infer:
    input:
        f'{projdir}/Tsinfer/samples/{{chromosome}}.samples',
    output:
        f'{projdir}/Tsinfer/trees/{{chromosome}}.trees',
    conda: "HLab_tsinfer"
    threads: 16
    resources: cpus=16, mem_mb=1024000, time_min=10080
    benchmark:
            'benchmarks/{chromosome}.infer.benchmark.txt'
    #log: 'logs/Infer_{chromosome}.log'
    shell:
        r"""
        # Check amount of memory (in kbytes) as seen by the job
        ulimit -v
        python scripts/Infer_trees.py {input} {output}
        """
