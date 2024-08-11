# Filter out INDELS and missing data (>20%) from the vcf file
# The input vcf file is the merged vcf file from the previous step
# or the output from split and move (if only one file is present)

rule remove_missing_indels:
    input: 
        vcf = input_vcf,
        index = input_index,
        summary = input_summary,
    output:    
        vcf = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}.recode.vcf.gz'),
        index = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}.recode.vcf.gz.csi')
    params:
        prefix = f'{vcfdir}/{{chromosome}}/{{chromosome}}',
        missing = config['missing_data'],
        indel = '--remove-indels' if config['remove_indel'] else ''
    conda: 'bcftools'
    benchmark: 'benchmarks/{chromosome}.missindel.benchmark.txt' 
    shell:
        r"""
        vcftools --gzvcf {input.vcf} \
            {params.indel} \
            --max-missing {params.missing} \
            --remove-filtered-all \
            --recode \
            --out {params.prefix}

        bgzip {params.prefix}.recode.vcf
        bcftools index -f {output.vcf}  

        # rm  {params.prefix}.recode.vcf 
        """ 
