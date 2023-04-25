vcfdir = "/home/jana/Documents/1Projects/HoneybeeDemo/Data"
chromosomes = [f'Chr{n}' for n in range(1, 3)]

splitFiles = list(set([f.strip("vcf.gz").split("_")[1] for f in os.listdir(f'{vcfdir}/RawVCF') if (f.endswith("vcf.gz") and f.startswith("Chr"))]))
combinedFiles = list(set([f.strip("vcf.gz").split("_")[1] for f in os.listdir(f'{vcfdir}/RawVCF') if (f.endswith("vcf.gz") and f.startswith("Combine"))]))

rule all:
    input:
        expand(f'{vcfdir}/{{chromosome}}/MergeList.yaml', chromosome = chromosomes)

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
            [f'{vcfdir}' + x for x in expand('/{{chromosome}}/{{chromosome}}_{{suffixTwo}}.vcf.gz', chromosome = chromosomes, suffixTwo = combinedFiles)]

    #    conda:
    #        'envs/vcfEdit.yaml'
        # log:
        #     logs:'{chromosome}_split.log'
        shell:
            """
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})

            bcftools view -r ${{chr}} {input} -O z -o {output}
            bcftools index {output}
            """

rule create_vcf_input:
    input:
        splitOutput=expand(f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{suffixOne}}.vcf.gz', chromosome = chromosomes, suffixOne = splitFiles) if len(splitFiles) > 0 else [],
        combinedOutput=expand(f'{vcfdir}/{{chromosome}}/{{chromosome}}_{{suffixTwo}}.vcf.gz', chromosome = chromosomes, suffixTwo = combinedFiles) if len(combinedFiles) > 0 else []
    output:
        [f'{vcfdir}' + x for x in expand('/{{chromosome}}/MergeList.yaml', chromosome = chromosomes)]
    run:
        import yaml
        chrDict = {wildcards.chromosome: [input.splitOutput, input.combinedOutput]}
        chrYaml = yaml.dump(chrDict)
        outfile = str(output[0]        f = open(outfile, "w")
        f.write(chrYaml)
        f.close()
