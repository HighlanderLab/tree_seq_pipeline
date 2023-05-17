import yaml
#vcfdir = "vcfDir"
#chromosomes = [f'Chr{n}' for n in range(1, 3)]

if len(allFiles) != 1:
    rule get_samples:
        input:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixOne}.vcf.gz',
                    suffixOne = splitFiles)] if len(splitFiles) != 0 else [],
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixTwo}.vcf.gz',
                    suffixTwo = combinedFiles)] if len(combinedFiles) != 0 else []
        output: temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.txt')
        shell:
            """
            bcftools query -l {input} > {output}
            """

#     # Executes per chromosome! Takes all .txt files at once, compare and remove
#     # duplicates, and write a new temporary .txt for each one.
    rule compare:
        input:
            samples = [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{file}.txt',
                    file=allFiles)],
            vcfs = [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{file}.vcf.gz',
                    file=allFiles)]
        output:
            temp([f'{vcfdir}' + x for x in
                expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz',
                    file=allFiles)])
# #         log: expand('logs/{{chromosome}}_{file}.log', file=allfiles)
        run:
            from itertools import combinations
            vcflist_input=input.vcfs
            vcflist_output=output
            samplelist=[]
            track=[]

            for ifile in input.samples:
                with open(ifile, 'r') as f:
                    samplelist.append(f.read().splitlines())


            # Take sublist in pairs (return index [0] and sublist [1]) and compare
            for a, b in combinations(enumerate(samplelist), 2):
                    # if there are overlapping entries, remove duplicates
                    if not set(a[1]).isdisjoint(b[1]) == True:
                        [b[1].remove(element) for element in a[1] if element in b[1]]
                        # keeps track of the sublists that changed using index
                        track.append(b[0])

            # filter files if sample list has changed, otherwise only rename
            for i in range(len(samplelist)):
                vcf=vcflist_input[i]
                print("VCFLIST")
                print(vcf)
                print("SAMPLES")
                print(samples)
                ovcf=vcflist_output[i]
                samples=samplelist[i]
                if i in set(track):
                    shell('bcftools view -S {samples} --force-samples {vcf} -O z -o {ovcf} \
                            bcftoools index {ovcf}')
                else:
                    shell('bcftools view {vcf} -O z -o {ovcf} \
                            bcftools index {ovcf}')

            # for i, ofile in enumerate(output):
            #     with open(ofile, 'w') as f:
            #         [f.write(f'{line}\n') for line in samplelist[i]]
            #         f.close()
#
# #     # Executes for all chromosomes all files at once. Takes .txt and .vcf as input
# #     # and filters the vcfs based on the new filtered samples list. index the output.
#     rule filter:
#         input:
#             vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.vcf.gz',
#             ids = f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.filtered'
#         output: f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.filtered.vcf.gz'
# #         conda:
# #             'envs/vcfEdit.yaml'
# #         log: 'logs/{chromosome}_{file}.log'
#         shell:
#             """
#             bcftools view -S {input.ids} --force-samples {input.vcf} -O z -o {output}
#             bcftools index {output}
#             """
#
# #     # Again takes all files for each chromosome at once.
# #     # Merges all vcfs and indexes them.
# #     # Outputs the final merged vcf for each chromosome directly in the vcfDir.
    rule merge:
        input: [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz', file=allFiles)]
        output: f'{vcfdir}/{{chromosome}}_final.vcf.gz'
#         conda:
#             'env/vcfEdit.yaml'
#         log: 'logs/{chromosome}_merged.log'
        shell:
            """
            bcftools merge {input} -O z -o {output}
            bcftools index {output}
            """

else:
    rule rename:
        input:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{file}.vcf.gz',
                    suffixOne = splitFiles)] if len(splitFiles) != 0 else [],
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{file}.vcf.gz',
                    suffixTwo = combinedFiles)] if len(combinedFiles) != 0 else []
        output:
            f'{vcfdir}/{{chromosome}}_final.vcf.gz'
        shell:
            """
            bcftools view {file} -O z -o {output}
            bcftools index {output}
            """
