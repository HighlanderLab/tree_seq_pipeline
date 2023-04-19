# import some python modules to define variables
from os import listdir
from os.path import isfile, join

configfile: '../config/tsinfer.yaml'
# files are expected to be stored in the vcfDir folder.
# files for each chromosome should be in a separate folder (for clarity).

# defining variables that will be used throughout the script
vcfdir=config['vcfDir'] # from config file, the path to vcfs
mapdir=config['genMap'] # from config file, the path to recombination maps
samples = [f'Chr{n}' for n in range(1, config['noChromosomes'] + 1)]
# set a list of 'Chr' + chr# combinations to be used as wildcards

allfiles = [f.split('.')[0].split('_')[1] for f in listdir(f'{vcfdir}/{samples[0]}')
if isfile(join(f'{vcfdir}/{samples[0]}', f)) and 'vcf' in f]
# defines a list of names for all different files
#groups = {s : [file for file in allfiles if s in file] for s in samples}

# At the end it should generate a phased vcf file that combines all files for
# each chromosome. These files are stored in the vcf folder declaired in the
# configuration file. It also checks that vcfs have been filtered before merged
rule all:
    input:
        expand(f'{vcfdir}/{{sample}}/{{sample}}_{{file}}.filtered.vcf.gz',
        sample=samples, file=allfiles),
        expand(f'{vcfdir}/{{sample}}_phased.vcf.gz', sample=samples)

# Starts by extracting a list of samples ids from each vcf.
# Executes for all chromosomes and all files at once. Outputs temporary .txt
# for each list.
rule get_samples:
    input: f'{vcfdir}/{{sample}}/{{sample}}_{{file}}.vcf.gz'
    output: temp(f'{vcfdir}/{{sample}}/{{sample}}_{{file}}.txt')
    conda:
        'envs/vcfEdit.yaml'
    log: 'logs/{sample}_{file}.log'
    shell:
        "bcftools query -l {input} > {output}"

# Executes per chromosome! Takes all .txt files at once, compare and remove
# duplicates, and write a new temporary .txt for each one.
rule compare:
    input: [f'{vcfdir}' + x for x in expand('/{{sample}}/{{sample}}_{file}.txt', file=allfiles)]
    output: temp([f'{vcfdir}' + x for x in expand('/{{sample}}/{{sample}}_{file}.filtered', file=allfiles)])
    log: expand('logs/{{sample}}_{file}.log', file=allfiles)
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
        vcf = f'{vcfdir}/{{sample}}/{{sample}}_{{file}}.vcf.gz',
        ids = f'{vcfdir}/{{sample}}/{{sample}}_{{file}}.filtered'
    output: f'{vcfdir}/{{sample}}/{{sample}}_{{file}}.filtered.vcf.gz'
    conda:
        'envs/vcfEdit.yaml'
    log: 'logs/{sample}_{file}.log'
    shell:
        """
        bcftools view -S {input.ids} --force-samples {input.vcf} -O z -o {output}
        bcftools index {output}
        """

# Again takes all files for each chromosome at once.
# Merges all vcfs and indexes them.
# Outputs the final merged vcf for each chromosome directly in the vcfDir.
rule merge:
    input: [f'{vcfdir}' + x for x in expand('/{{sample}}/{{sample}}_{file}.filtered.vcf.gz', file=allfiles)]
    output: f'{vcfdir}/{{sample}}_merged.vcf.gz'
    conda:
        'env/vcfEdit.yaml'
    log: 'logs/{sample}_merged.log'
    shell:
        """
        bcftools merge -m all {input} -O z -o {output}
        bcftools index {output}
        """

# Takes the final vcf and the recombination map as inputs and run shapeit4.
# using a specific environment file so installs shapeit4 if needed.
# Outputs the phased vcf.gz with index for all chromosomes. 
rule phase:
    input:
        vcf = rules.merge.output,
        map = f'{mapdir}/{{sample}}.gmap'

    output: f'{vcfdir}/{{sample}}_phased.vcf.gz'
    conda:
        'env/shapeit.yaml'
    log: 'logs/{sample}_phased.log'
    shell:
        """
        str='{wildcards.sample}'
        chr=$(echo ${{str:3}})

        shapeit4 --input {input.vcf} --map {input.map} --region ${{chr}} --output {output} --sequencing
        bcftools index {output}
        """
