import yaml


if len(allFiles) != 1:
    rule get_samples:
        input:
            f'{vcfOut}/{{chromosome}}/{{chromosome}}_{{file}}.vcf.gz'
        output: temp(f'{vcfOut}/{{chromosome}}/{{chromosome}}_{{file}}.txt')
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
            [f'{vcfOut}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{file}.txt',
                    file=allFiles)],
        output:
            temp([f'{vcfOut}' + x for x in
                expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.{ext}', file=allFiles, ext=['gz', 'gz.csi'])])
        params:
            duplicated = '{chromosome}.ids'
        resources: cpus=1, mem_mb=64000, time_min=60
        run:
            import os

            for file in input:
                shell("wc -l {file} >> files_len.txt")
            
            shell("sort -s -nr -k 1,1 files_len.txt | cut -d' ' -f 2 > files_reorderd.txt")

            with open('files_reorderd.txt', 'r') as file:
                files_reord = []
                for line in file:
                    line = line.strip()
                    files_reord.append(line)
            
            filtered = []

            for i in range(len(files_reord)):
                print(i)
                for j, file in enumerate(files_reord):
                    if file != files_reord[i] and j > i:
                        file1 = files_reord[i]
                        file2 = files_reord[j]
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

            for file in files_reord:
                if file not in filtered:
                    file1 = file.split('.')[0]
                    print(f'file {file1} not filtered')
                    shell("cp {file1}.vcf.gz {file1}.filtered.vcf.gz && bcftools index -f {file1}.filtered.vcf.gz")
            shell('rm {params.duplicated}')

    rule merge:
        input:
            [f'{vcfOut}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz', file=allFiles)]
        output:
            vcf=temp(f'{vcfOut}/{{chromosome}}_final.vcf.gz'),
            index=temp(f'{vcfOut}/{{chromosome}}_final.vcf.gz.csi')
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
            [f'{vcfOut}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixOne}.vcf.gz',
                    suffixOne = splitFiles)] if len(splitFiles) != 0 else [], 
                [f'{vcfOut}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixTwo}.vcf.gz',
                    suffixTwo = combinedFiles)] if len(combinedFiles) != 0 else []
        output:
            vcf = (f'{vcfOut}/{{chromosome}}_final.vcf.gz'),
            idx = (f'{vcfOut}/{{chromosome}}_final.vcf.gz.csi')
        conda: "bcftools"
        log: 'logs/Rename_{chromosome}.log'
        resources: cpus=1, mem_mb=32000, time_min=30
        shell:
            """
                ln -s $( realpath {input} ) {output.vcf}
                ln -s $( realpath {input} ).csi {output.idx}
            """

