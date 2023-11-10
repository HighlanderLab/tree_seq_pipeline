if int(config['ploidy']) == 1:
    rule rename_phased:
        input: f'{vcfdir}/{{chromosome}}_final.vcf.gz'
        output:
            file = f'{vcfdir}/{{chromosome}}_phased.vcf.gz',
            idx = f'{vcfdir}/{{chromosome}}_phased.vcf.gz.csi'
        log: 'logs/Rename_phased_{chromosome}.log'
        shell:
            """
            bcftools view {input} -O z -o {output.file}
            bcftools index {output.file}
            """

else:
    rule phase:
        input:
            vcf = f'{vcfdir}/{{chromosome}}_final.vcf.gz',
        output:
            file = f'{vcfdir}/{{chromosome}}_phased.vcf.gz',
            idx = f'{vcfdir}/{{chromosome}}_phased.vcf.gz.csi'
        params:
            map = f'{mapdir}/{{chromosome}}.txt',
        #log: 'logs/phase_{chromosome}.log'
        threads: 25
        #resources: cpus=20, mem_mb=25000, time_min=5
        conda: 'shapeit4am'
        log: 'logs/Phase_{chromosome}.log'
        shell:
            """
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})
            start=`date +%s`
            shapeit4 --input {input.vcf} \
                             --map {params.map} \
                             --region ${{chr}} \
                             --output {output} \
                             --thread {threads}
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds > shapeit_{wildcards.chromosome}.time
            bcftools index {output}
            """
# --sequencing
