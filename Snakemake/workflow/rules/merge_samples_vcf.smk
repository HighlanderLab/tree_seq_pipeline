import yaml
#vcfdir = "vcfDir"
#chromosomes = [f'Chr{n}' for n in range(1, 3)]

if len(allFiles) != 1:
    rule get_samples:
        input:
            f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.vcf.gz'
        output: temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.txt')
        # envmodules:
        #     config['bcftoolsModule']
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
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
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        shell:
            "python scripts/CompareVCFs.py"
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
        input:
            [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz', file=allFiles)]
        output:
            vcf=f'{vcfdir}/{{chromosome}}_final.vcf.gz',
            index=f'{vcfdir}/{{chromosome}}_final.vcf.gz.csi'
#         conda:
#             'env/vcfEdit.yaml'
#         log: 'logs/{chromosome}_merged.log'
        # envmodules:
        #     config['bcftoolsModule']
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        shell:
            """
            bcftools merge {input} -O z -o {output.vcf}
            bcftools index {output.vcf}
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
        # envmodules:
        #     config['bcftoolsModule']
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        shell:
            """
            bcftools view {input} -O z -o {output} #GABRIELA: SHOULD THIS BE INPUT????
            bcftools index {output}
            """
