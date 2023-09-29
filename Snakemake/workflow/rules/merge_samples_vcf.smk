import yaml


if len(allFiles) != 1:
    rule get_samples:
        input:
            f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.vcf.gz'
        output: temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.txt')
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        log: 'logs/get_samples_{chromosome}_{file}.log'
        shell:
            """
            bcftools query -l {input} > {output}
            """

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
#        log: 'logs/compare_{chromosome}_{file}.log'
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        shell:
            """
            python scripts/CompareVCFs.py {input.samples}  {input.vcfs}  {output}
            """

    rule merge:
        input:
            [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz', file=allFiles)]
        output:
            vcf=f'{vcfdir}/{{chromosome}}_final.vcf.gz',
            index=f'{vcfdir}/{{chromosome}}_final.vcf.gz.csi'
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        log: 'logs/merge_{chromosome}.log'
        shell:
            """
            bcftools merge {input} -O z -o {output.vcf}
            bcftools index {output.vcf}
            """

else:
    rule rename:
        input:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixOne}.vcf.gz',
                    suffixOne = splitFiles)] if len(splitFiles) != 0 else [],
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixTwo}.vcf.gz',
                    suffixTwo = combinedFiles)] if len(combinedFiles) != 0 else []
        output:
            f'{vcfdir}/{{chromosome}}_final.vcf.gz'
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        log: 'logs/rename_{chromosome}.log'
        shell:
            """
            bcftools view {input} -O z -o {output}
            bcftools index {output}
            """
