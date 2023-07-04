
splitFiles = list(set([f.strip("vcf.gz").split("_")[1]
    for f in os.listdir(f'{vcfdir}/RawVCF')
        if (f.endswith("vcf.gz") and f.startswith("Chr"))]))
combinedFiles = list(set([f.strip("vcf.gz").split("_")[1]
    for f in os.listdir(f'{vcfdir}/RawVCF')
        if (f.endswith("vcf.gz") and f.startswith("Combined"))]))
allFiles = splitFiles + combinedFiles

if len(splitFiles) > 0:
    rule move_vcf:
        input:
            f'{vcfdir}/RawVCF/{{chromosome}}_{{suffixOne}}.vcf.gz'
        output:
            f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{suffixOne}}.vcf.gz'
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        log: 'logs/move_vcf_{chromosome}.log'
        shell:
            """
            echo {input}
            cp {input} {output}
            bcftools index {output}
            """

if len(combinedFiles) > 0:
    rule split_and_move_vcfs:
        input:
            f'{vcfdir}/RawVCF/Combined_{{suffixTwo}}.vcf.gz'
        output:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{{suffixTwo}}.vcf.gz',
                    chromosome = chromosomes, suffixTwo = combinedFiles)]
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=4000, time_min=5
        log: expand('logs/split_and_move_vcfs_{{chromosome}}_{{suffixTwo}}.log')
        shell:
            """
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})

            bcftools view -r ${{chr}} {input} -O z -o {output}
            bcftools index {output}
            """
