
splitFiles = list(set([f.strip("vcf.gz").split("_")[1]
    for f in os.listdir(f'{vcfdir}/RawVCF')
        if (f.endswith("vcf.gz") and (f.startswith("Chr") | f.startswith("chr")) )]))
combinedFiles = list(set([f.strip("vcf.gz").split("_")[1]
    for f in os.listdir(f'{vcfdir}/RawVCF')
        if (f.endswith("vcf.gz") and f.startswith("Combined"))]))
allFiles = splitFiles + combinedFiles

if len(splitFiles) > 0:
    rule move_vcf:
        input:
            vcf = f'{vcfdir}/RawVCF/{{chromosome}}_{{suffixOne}}.vcf.gz',
            idx = f'{vcfdir}/RawVCF/{{chromosome}}_{{suffixOne}}.vcf.gz.csi'
        output:
            vcf = f'{vcfOut}/{{chromosome}}/{{chromosome}}_{{suffixOne}}.vcf.gz',
            idx = f'{vcfOut}/{{chromosome}}/{{chromosome}}_{{suffixOne}}.vcf.gz.csi'
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=60
        log: 'logs/Move_vcf_{chromosome}_{suffixOne}.log'
        shell:
            """
            if [ -h {input.vcf} ]; then
                ln -s $( realpath {input.vcf} ) {output.vcf}
                ln -s $( realpath {input.vcf} ).csi {output.vcf}
            else
                ln -s {input.vcf} {output.vcf}
                ln -s {input.idx} {output.idx}
            fi
            """

if len(combinedFiles) > 0:
    rule split_and_move_vcfs:
        input:
            f'{vcfdir}/RawVCF/Combined_{{suffixTwo}}.vcf.gz'
        output:
            vcf = [f'{vcfOut}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{{suffixTwo}}.vcf.gz',
                    chromosome = chromosomes, suffixTwo = combinedFiles)],
            idx = [f'{vcfOut}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{{suffixTwo}}.vcf.gz.csi',
                    chromosome = chromosomes, suffixTwo = combinedFiles)]
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=64000, time_min=60
        log: expand('logs/Split_and_move_vcfs_{{chromosome}}_{{suffixTwo}}.log')
        shell:
            """
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})

            bcftools view -r ${{chr}} {input} -O z -o {output.vcf}
            bcftools index -f {output.vcf}
            """
