import yaml


if len(allFiles) != 1:
    rule get_samples:
        input:
            f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.vcf.gz'
        output: temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{file}}.txt')
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=30
        log: 'logs/Get_samples_{chromosome}_{file}.log'
        shell:
            """
            bcftools query -l {input} > {output}
            """

    rule filter:
        input:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{file}.txt',
                    file=allFiles)],
        output:
            temp([f'{vcfdir}' + x for x in
                expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.{ext}', file=allFiles, ext=['gz', 'gz.csi'])])
        params:
            duplicated = '{chromosome}.ids'
       # conda: 'bcftools'
        resources: cpus=1, mem_mb=64000, time_min=60
        run:
            import os
            filtered = []

            for i in range(len(input)):
                print(i)
                for j, file in enumerate(input):
                    if file != input[i] and j > i:
                        file1 = input[i]
                        file2 = input[j]
                        print(file1, file2)

                        shell("comm -12 <(sort {file1}) <(sort {file2}) > {params.duplicated}")

                        remove = params.duplicated
                        lremove = os.stat(f'{remove}').st_size

                        if lremove == 0:
                            print(f'Nothing to remove {lremove}')
                        else:
                            print(f'{lremove} samples to remove')
                            if file2 not in filtered:
                                filtered.append(file2)

                            file2 = file2.split('.')[0]
                            print(f'file {file2} filtered')
                            shell('bcftools view -S ^{params.duplicated} {file2}.vcf.gz -O z -o {file2}.filtered.vcf.gz')
                            shell('bcftools index -f {file2}.filtered.vcf.gz')
                            shell('rm {params.duplicated}')

            for file in input:
                if file not in filtered:
                    file1 = file.split('.')[0]
                    print(f'file {file1} not filtered')
                    shell("cp {file1}.vcf.gz {file1}.filtered.vcf.gz && bcftools index -f {file1}.filtered.vcf.gz")
            shell('rm {params.duplicated}')

    rule merge:
        input:
            [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz', file=allFiles)]
        output:
            vcf=temp(f'{vcfdir}/{{chromosome}}_final.vcf.gz'),
            index=temp(f'{vcfdir}/{{chromosome}}_final.vcf.gz.csi')
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=60
        log: 'logs/Merge_{chromosome}.log'
        shell:
            """
            bcftools merge {input} -O z -o {output.vcf}
            bcftools index -f {output.vcf}
            """

else: # i.e. len(allFiles) == 1
    rule rename:
        input:
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixOne}.vcf.gz',
                    suffixOne = splitFiles)] if len(splitFiles) != 0 else [],
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixTwo}.vcf.gz',
                    suffixTwo = combinedFiles)] if len(combinedFiles) != 0 else []
        output:
            vcf = (f'{vcfdir}/{{chromosome}}_final.vcf.gz'),
            idx = (f'{vcfdir}/{{chromosome}}_final.vcf.gz.csi')
        conda: "bcftools"
        log: 'logs/Rename_{chromosome}.log'
        resources: cpus=1, mem_mb=32000, time_min=30
        shell:
            """
            if [ -h {input} ]; then
                ln -s $( realpath {input} ) {output.vcf}
                ln -s $( realpath {input} ).csi {output.vcf}
            else
                ln -s {input} {output.vcf}
                ln -s {input}.csi {output.idx}
            fi
            """

