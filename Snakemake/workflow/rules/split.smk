rule split_vcfs:
    input:
        expand(f'{vcfdir}/{{vcfcombined}}',
            vcfcombined=VCFs)
    output:
        [f'{vcfdir}' + x
            for x in expand('/{{chromosome}}/{{chromosome}}_split.vcf.gz',
                chromosome=chromosomes)]
    params:
        dir_dump=directory(f'{vcfdir}/combined/')
#    conda:
#        'envs/vcfEdit.yaml'
    log:
        logs:'{chromosome}_split.log'
    shell:
        """
        str='{wildcards.chromosome}'
        chr=$(echo ${{str:4}})

        bcftools view -r ${{chr}} {input} -O z -o {output}
        bcftools index {output}
        """
