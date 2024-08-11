# This script is used for managing VCF files in a directory previous to tree inference
# It identifies two types of VCF files in the 'RawVCF' directory: 
# 1. Files that start with either "Chr" or "chr" and end with "vcf.gz" (referred to as 'splitFiles')
# 2. Files that start with "Combined" and end with "vcf.gz" (referred to as 'combinedFiles')
# It then creates regex patterns for these file names for later use in rules.
# The script also combines the names of both types of files into a single list (allFiles).
# If there are any 'splitFiles', it defines a rule (move_vcf) to move these files 
# from the 'RawVCF' directory to a directory named after the chromosome in the file name.
# The rule also renames the files to include the chromosome name at the start of the file name.
import os
import re

print(chromosomes)
splitFiles = list(set([f.strip("vcf.gz").split("_")[1]
    for f in os.listdir(f'{vcfdir}/RawVCF')
        if (f.endswith("vcf.gz") and (f.startswith("Chr") | f.startswith("chr")) and "_" in f)]))

combinedFiles = list(set([f.strip("vcf.gz").split("_")[1]
    for f in os.listdir(f'{vcfdir}/RawVCF')
        if (f.endswith("vcf.gz") and f.startswith("Combined"))]))

# Prepare a regex pattern that matches any of the file names in splitFiles and
# combinedFiles so rules only look for them
splitFiles_pattern = "(" + "|".join(re.escape(name) for name in splitFiles) + ")"
combinedFiles_pattern = "(" + "|".join(re.escape(name) for name in combinedFiles) + ")"        

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
        benchmark: 'benchmarks/{chromosome}.{suffixOne}.move.benchmark.txt'
        wildcard_constraints:
            suffixOne = splitFiles_pattern
        shell:
            """
            ln -s {input.vcf} {output.vcf}
            ln -s {input.idx} {output.idx}
            """

if len(combinedFiles) > 0:
    rule split_and_move_vcfs:
        input:
            f'{vcfdir}/RawVCF/Combined_{{suffixTwo}}.vcf.gz'
        output:
            temp([f'{vcfdir}' + x
                for x in expand('/{{chromosome}}/{{chromosome}}_{{suffixTwo}}.vcf.gz',
                    chromosome = chromosomes, suffixTwo = combinedFiles)])
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=64000, time_min=60
        benchmark: 'benchmarks/{chromosome}.{suffixTwo}.split.benchmark.txt'
        wildcard_constraints:
            suffixTwo = combinedFiles_pattern
        shell:
            """
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})

            bcftools view -r ${{chr}} {input} -O z -o {output.vcf}
            bcftools index -f {output.vcf}
            """
