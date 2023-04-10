#
#
#
#
import itertools
import numpy as np
configfile: 'config.yaml'

dir = config['vcfDir']
# single file patterns
(PREFIX, SUFFIX) = glob_wildcards(dir+'/{prefix}_{suffix}.vcf')
#print(PREFIX)


rule all:
    input:
        expand('{prefix}_{suffix}', zip, prefix=(PREFIX), suffix=SUFFIX)
        #expand('Chr{chromosome}.vcf.gz', chromosome=range(1, config['noChromosomes'] + 1))


rule extract_samples_list:
    input: dir + '/{prefix}_{suffix}.vcf'
    output: '{prefix}_{suffix}.txt'
    shell:
        """
        bcftools query -l {input} > {output}
        """

#rule extract_samples_list:
#    input:
#        expand(dir + '/{prefix}_{suffix}.vcf', zip, prefix=PREFIX, suffix=SUFFIX)
#        files = lambda wildcards: expand(
#            ['{files}'], files = file_dic[int(wildcards.chromosome)]
#        )
#    output:
#        expand('{prefix}_{suffix}.txt', zip, prefix=PREFIX, suffix=SUFFIX)
        #expand('ids{{chromosome}}_{i}', i = range(1, count + 1))
#    shell:
#        """
#        vcfs=( {input} )
#        outfiles=( {output} )

#        for ((i=0; i<${{#vcfs[@]}}; i++)); do
#            bcftools query -l ${{vcfs[$i]}} > ${{outfiles[$i]}}
#        done

#        """

rule compare_samples:
    input:
#    # takes input from previous rule
        expand('{prefix}_{suffix}.txt', zip, prefix=PREFIX, suffix=SUFFIX)
    output:
    # create output based on combination of wildcards
        #expand('{{prefix}}_{suffix}_filtered.txt', zip, prefix=PREFIX, suffix=SUFFIX)
        expand('{prefix}_{suffix}_filtered.txt', zip, prefix=PREFIX, suffix=SUFFIX)
    run:
#    # run in python
        prefix = list(set([file.split('_')[0] for file in input]))
        # group input and output files by chromosome

        for name in prefix:
        # execute for each chromosome group
            infiles = [file for file in input if name in file] # define input files
            outfiles = [file for file in output if name in file] # define output files

            samplelist=[] # creates holder to read in content of input files (this is updated every loop)

            for ifile in infiles:  # read the content of the files and append to list (each file content is a list within list)
                with open(ifile, 'r') as f:
                    samplelist.append(f.read().splitlines())

            for a, b in itertools.combinations(samplelist,2): # compare the lists in samplelist by pair
                [b.remove(element) for element in a if element in b] # if id exists already in firs list, remove from second

            for i, ofile in enumerate(outfiles): # write out filtered lists
                with open(ofile, 'w') as f:
                    [f.write(f'{line}\n') for line in samplelist[i]]
                    f.close()

rule filter_samples_in_vcf:
    input:
        file = dir + '/{prefix}_{suffix}.vcf',
        id = '{prefix}_{suffix}_filtered.txt'
    output: '{prefix}_{suffix}.vcf.gz'
    shell:
        """
        bcftools view -S {input.id} --force-samples {input.file} -O z -o {output}
        bcftools index {output}
        """

rule merge_vcfs:
    input:
        files=expand('{{prefix}}_{suffix}.vcf.gz', suffix=set(SUFFIX))
        #direct={dir}
    output: directory(expand('{{prefix}}_{suffix}', suffix=set(SUFFIX)))
    shell:
        """
        files='{output}'
        dir='{config[vcfDir]}'
        echo ${{files}}

        bcftools merge -m all {input.files} -O z -o ${{dir}}/{wildcards.prefix}.vcf.gz
        bcftools index ${{dir}}/{wildcards.prefix}.vcf.gz

        mkdir {output}

        for f in ${{files}}; do
            mv ${{f}}*.* ${{dir}}/${{f}}*.* ${{f}}
        done
        """
