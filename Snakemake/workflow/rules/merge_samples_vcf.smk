import yaml
#vcfdir = "vcfDir"
#chromosomes = [f'Chr{n}' for n in range(1, 3)]

mergeList = f'{vcfdir}/Chr1/MergeList.yaml'
with open(mergeList) as f:
    files2merge = yaml.safe_load(f)

files=[f for f in files2merge.values()][0]
nfiles2merge = len(files)
names = [f.strip('.vcf.gz').split('/')[-1].split('_')[-1] for f in files]

# rule all:
#  input: expand(f'{vcfdir}/{{chromosome}}_final.vcf.gz', chromosome=chromosomes)

if nfiles2merge != 1:
    print(f'{nfiles2merge} files to merge')
    rule get_samples:
        input:
            files=f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{name}}.vcf.gz'
        output: temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{name}}.txt')
        shell:
            """
            bcftools query -l {input.files} > {output}
            """

#     # Executes per chromosome! Takes all .txt files at once, compare and remove
#     # duplicates, and write a new temporary .txt for each one.
    rule compare:
        input:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{name}.txt',
                    name=names)]
        output:
            temp([f'{vcfdir}' + x for x in
                expand('/{{chromosome}}/{{chromosome}}_{name}.filtered',
                    name=names)])
# #         log: expand('logs/{{chromosome}}_{file}.log', file=allfiles)
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
#
# #     # Executes for all chromosomes all files at once. Takes .txt and .vcf as input
# #     # and filters the vcfs based on the new filtered samples list. index the output.
    rule filter:
        input:
            vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{name}}.vcf.gz',
            ids = f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{name}}.filtered'
        output: f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{name}}.filtered.vcf.gz'
#         conda:
#             'envs/vcfEdit.yaml'
#         log: 'logs/{chromosome}_{file}.log'
        shell:
            """
            bcftools view -S {input.ids} --force-samples {input.vcf} -O z -o {output}
            bcftools index {output}
            """
#
# #     # Again takes all files for each chromosome at once.
# #     # Merges all vcfs and indexes them.
# #     # Outputs the final merged vcf for each chromosome directly in the vcfDir.
    rule merge:
        input: [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{name}.filtered.vcf.gz', name=names)]
        output: f'{vcfdir}/{{chromosome}}_final.vcf.gz'
#         conda:
#             'env/vcfEdit.yaml'
#         log: 'logs/{chromosome}_merged.log'
        shell:
            """
            bcftools merge -m all {input} -O z -o {output}
            bcftools index {output}
            """

else:
    rule rename:
        input:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/MergeList.yaml',
                    chromosome = chromosomes)]
        output: f'{vcfdir}/{{chromosome}}_final.vcf.gz'
        run:
            with open(mergeList) as f:
                dict = yaml.safe_load(f)

            dict_values = [f for f in files2merge.values()]
            file = [item for sublist in dict_values for item in sublist]

            shell("mv {file} {output}")
