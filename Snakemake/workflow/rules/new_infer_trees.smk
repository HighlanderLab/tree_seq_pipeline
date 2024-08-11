# Convert vcf to samples file
# generate ancestos samples based only on sites with ancestral allele
# information
# truncate ancestors to remove those spanning too many sites
# generate ancestors trees
# match ancestors and samples to infer trees
# final output is a simplified tree sequence

rule prepare_sample_file:
    input:
        vcf = input_ts,
        # f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf.gz',
        meta = config['meta']
    output: f'{projdir}/Tsinfer/samples/{{chromosome}}.samples'
    params:
        chrLength = lambda wildcards:  config['chromosome_length'][int(re.findall(r'\d+', wildcards.chromosome)[0])]
    conda: "HLab_tsinfer"
    threads: 32
    resources: cpus = 32, mem_mb = 768000, time_min = 200
    benchmark: 'benchmarks/{chromosome}_new.prep_sample.benchmark.txt'
    shell:
        r"""
        ulimit -v
        
        python scripts/generate_samples.py \
                {input.vcf} {output} {input.meta} {params.chrLength} {threads}
        
        """

rule generate_ancestors:
    input: 
        samples = rules.prepare_sample_file.output,
        sites_to_exclude = f'{vcfdir}/{{chromosome}}/{{chromosome}}.exclude',
    output:
        ancestors_samples = f'{projdir}/Tsinfer/samples/{{chromosome}}.ancestors',
    conda: 'HLab_tsinfer'
    threads: 64
    resources: cpus = 64, mem_mb = 2048000, time_min = 1680
    benchmark: 'benchmarks/{chromosome}_generate_ancestors.benchmark.txt'
    shell:
        r"""
        # Check amount of memory (in kbytes) as seen by the job
        ulimit -v
        python scripts/generate_ancestors.py {input.samples} {input.sites_to_exclude} {output.ancestors_samples} {threads}
        """

rule truncate_ancestors:
    input:
        ancestors_samples = rules.generate_ancestors.output.ancestors_samples,
    output: 
        truncated_ancestors = f'{projdir}/Tsinfer/samples/{{chromosome}}.ancestors.truncated',
    params:
       upper = config['truncate_upper'],
       lower = config['truncate_lower'],
    conda: 'HLab_tsinfer'
    threads: 32
    resources: cpus=32, mem_mb=2048000, time_min=1680
    benchmark: 'benchmarks/{chromosome}_truncate.benchmark.txt'
    shell:
        r"""
        ulimit -v
        python scripts/truncate.py {input} {params.upper} {params.lower} {output}
        """

rule ancestors_trees:
    input:
        samples = rules.prepare_sample_file.output,
        truncated_ancestors = rules.truncate_ancestors.output.truncated_ancestors,
    output: 
        ancestors_trees = f'{projdir}/Tsinfer/trees/{{chromosome}}.atrees',
    params:
        re = f"--recombination-rate {config['recombination_rate_anc']}" if config['recombination_rate_anc'] else '',
        misma = f"--mismatch-ratio {config['mismatch_ratio_anc']}" if config['mismatch_ratio_anc'] else '',
    conda: 'HLab_tsinfer'
    threads: 32
    resources: cpus = 32, mem_mb = 2048000, time_min = 1680
    benchmark: 'benchmarks/{chromosome}_ancestors_trees.benchmark.txt'
    shell:
        r"""
        tsinfer match-ancestors {input.samples} \
                -A {output} \
                {params.re} \
                {params.misma} \
                --log-section tsinfer.inference \
                --num-threads {threads} \
                --progress \
                --verbosity \
                --ancestors {input.truncated_ancestors} 
        """

rule infer_trees:
    input:
        samples = rules.prepare_sample_file.output,
        ancestors_trees = rules.ancestors_trees.output.ancestors_trees,
    output:
        trees = f'{projdir}/Tsinfer/trees/{{chromosome}}.trees',
    params:
        re = f"--recombination-rate {config['recombination_rate_tree']}" 
            if config['recombination_rate_tree'] else '',
        misma = f"--mismatch-ratio {config['mismatch_ratio_tree']}" 
            if config['mismatch_ratio_tree'] else '',
    conda: 'HLab_tsinfer'
    threads: 32
    resources: cpus=32, mem_mb=2048000, time_min=1680
    benchmark: 'benchmarks/{chromosome}_infer.benchmark.txt'
    shell:
        r"""
        # Check amount of memory (in kbytes) as seen by the job
        ulimit -v
        tsinfer match-samples {input.samples} \
                -A {input.ancestors_trees} \
                --num-threads {threads} \
                {params.re} \
                {params.misma} \
                --progress \
                --log-section tsinfer.inference \
                --verbosity \
                -O {output}
        """

rule simplify_trees:
    input: 
        trees = rules.infer_trees.output.trees,
    output: 
        f'{projdir}/Tsinfer/trees/{{chromosome}}.dated.trees',
    params:
        genint = config['generation_interval'],
        mutrate = config['mutation_rate'],
    conda: 'HLab_tsinfer'
    threads: 32
    resources: cpus=32, mem_mb=2048000, time_min=1680
    benchmark: 'benchmarks/{chromosome}_simplify.benchmark.txt'
    shell: 
        r"""
        ulimit -v
        python scripts/simplify.py {input} {params.genint} {params.mutrate} {output}
        """             

# ,[^.]        