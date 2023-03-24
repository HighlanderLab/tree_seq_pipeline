#
#
#
#
import itertools
import numpy as np
configfile: 'test.yaml'

def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]

def dic_files_per_chromosome(direc_lst, noChromosomes):
    file_dic={}
    for i in range(1, noChromosomes + 1):
        file_dic[i] = [name for name in direc_lst if str(i) in name]

    return file_dic

rule all:
    input:
        expand('Chr{chromosome}.vcf.gz', chromosome=range(1, config['noChromosomes'] + 1))

#workdir: config['vcfDir']

VCFs =[x for x in list_full_paths(config['vcfDir']) if x.endswith(".vcf")]
print('here are vcfs', VCFs)

file_dic = dic_files_per_chromosome(VCFs, config['noChromosomes'])
count = list(set(map(len, file_dic.values())))[0]
#print(file_dic)
#print(count)

rule extract_samples_list:
    input:
        files = lambda wildcards: expand(
            ['{files}'], files = file_dic[int(wildcards.chromosome)]
        )
    output:
        expand('ids{{chromosome}}_{i}', i = range(1, count + 1))
    shell:
        """
        vcfs=( {input.files} )
        outfiles=( {output} )

        for ((i=0; i<${{#vcfs[@]}}; i++)); do
            bcftools query -l ${{vcfs[$i]}} > ${{outfiles[$i]}}
        done
        """

rule compare_samples:
    input:
        rules.extract_samples_list.output
    output:
        expand('ids{{chromosome}}_filtered_{i}', i = range(1, count + 1))
    run:
        samplelist=[]
        for file in input:
            with open(file, 'r') as f:
                samplelist.append(f.read().splitlines())

        #samplelist.append(list(np.loadtxt(lst, dtype = str))) for lst in input]
        print(len(samplelist))

        for a, b in itertools.combinations(samplelist, 2):
            [b.remove(element) for element in a if element in b]


        for i, file in enumerate(output):
            with open(file, 'w') as f:
                for line in samplelist[i]:
                    f.write(f'{line}\n')
            f.close()

rule filter_samples_in_vcf:
    input:
        files = lambda wildcards: expand(
            ['{files}'], files = file_dic[int(wildcards.chromosome)]
        ),
        ids  = rules.compare_samples.output
    output:
        expand('Chr{{chromosome}}_{i}.vcf.gz', i = range(1, count + 1))
    shell:
        """
        vcfs=( {input.files} )
        ids=( {input.ids} )
        outfiles=( {output} )

        for ((i=0; i<${{#vcfs[@]}}; i++)); do
            bcftools view -S ${{ids[$i]}} --force-samples ${{vcfs[$i]}} -O z -o ${{outfiles[$i]}}
            bcftools index ${{outfiles[$i]}}
        done
        """

# Problem here
rule merge_vcfs:
    input:
        rules.filter_samples_in_vcf.output
    output:
        'Chr{chromosme}.vcf.gz'
    shell:
        """
        bcftools merge -m all {input} -O z -o {output}
        bcftools index {output}
        """

# This is just for tests so it runs without complaining
#rule testing:
#    input:
#        rules.merge_vcfs.output
        #ids = rules.extract_samples_list.output
#    output:
#        'Chr{chromosome}.vcf.gz'
#    shell:
#        """
#        echo {output}
#        """
