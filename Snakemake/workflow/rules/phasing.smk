if config['ploidy'] == 1:
    rule rename_phased:
        input: rules.compress_vcf.output
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
            vcf = rules.compress_vcf.output,
            map = '../mapsDir/{chromosome}.gmap'

        output: f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
        log: 'logs/phase_{chromosome}.log'
        threads: 20
        resources: cpus=1, mem_mb=4000, time_min=5
        conda: 'shapeit4am'
        shell:
            """
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:4}})
            shapeit4 --input {input.vcf} \
                             --map {input.map} \
                             --region ${{chr}} \
                             --output {output} \
                             --sequencing \
                             --thread 20
            bcftools index {output}
            """
