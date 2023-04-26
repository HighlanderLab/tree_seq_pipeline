# ---------------------------------------------------------------------------- #
#                               MERGE                                          #
# ---------------------------------------------------------------------------- #

# Starts by extracting a list of samples ids from each vcf.
# Executes for all chromosomes and all files at once. Outputs temporary .txt
# for each list.

rule get_samples:
    input:
        files=f'{vcfdir}/{{chromosome}}/{{file}}.vcf.gz'
    output: temp(f'{vcfdir}/{{chromosome}}/{{file}}.txt')
    conda:
        'envs/vcfEdit.yaml'
    log: 'logs/{chromosome}_{file}.log'
    shell:
        """
        bcftools query -l {input.files} > {output}
        """

# Executes per chromosome! Takes all .txt files at once, compare and remove
# duplicates, and write a new temporary .txt for each one.
rule compare:
    input: [f'{vcfdir}/{{chromosome}}/{{file}}.txt'
    output: temp([f'{vcfdir}/{{chromosome}}/{{file}}.filtered')
    log: expand('logs/{{chromosome}}_{file}.log', file=allfiles)
    run:
        from itertools import combinations
        samplelist=[]
        for ifile in input:
            with open(ifile, 'r') as f:
                samplelist.append(f.read().splitlines())

        for a, b in combinations(samplelist, 2):
            [b.remove(element) for element in a if element in b]

        for i, ofile in enumerate(output):
            with open(ofile, 'w') as f:
                [f.write(f'{line}\n') for line in samplelist[i]]
                f.close()

# Executes for all chromosomes all files at once. Takes .txt and .vcf as input
# and filters the vcfs based on the new filtered samples list. index the output.
rule filter:
    input:
        vcf = f'{vcfdir}/{{chromosome}}/{{file}}.vcf.gz',
        ids = f'{vcfdir}/{{chromosome}}/{{file}}.filtered'
    output: f'{vcfdir}/{{chromosome}}/{{file}}.filtered.vcf.gz'
    conda:
        'envs/vcfEdit.yaml'
    log: 'logs/{chromosome}_{file}.log'
    shell:
        """
        bcftools view -S {input.ids} --force-samples {input.vcf} -O z -o {output}
        bcftools index {output}
        """

# Again takes all files for each chromosome at once.
# Merges all vcfs and indexes them.
# Outputs the final merged vcf for each chromosome directly in the vcfDir.
rule merge:
    input: [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{file}}.filtered.vcf.gz']
    output: f'{vcfdir}/{{chromosome}}.vcf.gz'
    conda:
        'env/vcfEdit.yaml'
    log: 'logs/{chromosome}_merged.log'
    shell:
        """
        bcftools merge -m all {input} -O z -o {output}
        bcftools index {output}
        """
