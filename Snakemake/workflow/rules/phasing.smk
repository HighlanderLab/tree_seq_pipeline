if config['ploidy'] == 1:
    rule rename_phased:
        input: f'{vcfdir}/{{chromosome}}_final.vcf.gz'
        output: f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
        log: 'logs/rename_phased_{chromosome}.log'
        shell:
            """
            bcftools view {input} -O z -o {output}
            bcftools index {output}
            """

else:
    rule phase:
        input:
            vcf = f'{vcfdir}/{{chromosome}}_final.vcf.gz',
        output: f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
        params: 
            map = f'{mapdir}/{{chromosome}}.txt',
        log: 'logs/phase_{chromosome}.log'
        threads: 20
        resources: cpus=20, mem_mb=25000, time_min=5
        conda: 'shapeit4am'
        shell:
            """
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})
            shapeit4 --input {input.vcf} \
                             --map {params.map} \
                             --region ${{chr}} \
                             --output {output} \
                            # --sequencing \
                             --thread 20
            bcftools index {output}
            """
