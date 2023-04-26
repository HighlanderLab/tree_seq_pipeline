vcfdir = "vcfDir"
chromosomes = [f'Chr{n}' for n in range(1, 3)]

# FILE NAMES MUST BE:
#    SPLITFILES = CHR{N}_ORIGIN.VCF.GZ
#    COMBINEDFILES = COMBINE_ORIGIN.VCF.GZ
splitFiles = list(set([f.strip("vcf.gz").split("_")[1]
    for f in os.listdir(f'{vcfdir}/RawVCF')
        if (f.endswith("vcf.gz") and f.startswith("Chr"))]))
combinedFiles = list(set([f.strip("vcf.gz").split("_")[1]
    for f in os.listdir(f'{vcfdir}/RawVCF')
        if (f.endswith("vcf.gz") and f.startswith("Combine"))]))

rule all:
    input:
        expand(f'{vcfdir}/{{chromosome}}/MergeList.yaml',
            chromosome = chromosomes)

if len(splitFiles) > 0:
    rule move_vcf:
        input:
            f'{vcfdir}/RawVCF/{{chromosome}}_{{suffixOne}}.vcf.gz'
        output:
            f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{suffixOne}}.vcf.gz'
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

        shell:
            """
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})

            bcftools view -r ${{chr}} {input} -O z -o {output}
            bcftools index {output}
            """

# looks terrible but it is now running by chromosome
rule create_vcf_input:
    input:
        splitOutput=
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixOne}.vcf.gz',
                    suffixOne = splitFiles)] if len(splitFiles) != 0 else [],
        combinedOutput=
            [f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{suffixTwo}.vcf.gz',
                suffixTwo = combinedFiles)] if len(combinedFiles) != 0 else []
    output:
        [f'{vcfdir}' + x for x in expand('/{{chromosome}}/MergeList.yaml',
            chromosome = chromosomes)]
    run:
        import yaml

        name = str(wildcards.chromosome)
        filelst = [str(f) for f in input.splitOutput] + [str(f) for f in input.combinedOutput]

        chrDict = {name: filelst}
        outfile = str(output[0])
        f = open(outfile, 'w+')
        yaml.dump(chrDict, f, allow_unicode=False, default_flow_style=False)
        f.close()
