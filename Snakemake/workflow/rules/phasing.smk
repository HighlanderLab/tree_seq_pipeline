# def getInputVcfFile_phasing(wildcards):
#     if os.path.isfile(vcfdir + "/" + wildcards.chromosome + '_final.vcf.gz'):
#         print("Final file found ")
#         vcf_file = vcfdir + "/" + wildcards.chromosome + '_final.vcf.gz'
#         print(vcf_file)
#     else:
#         print("No final file ")
#         file = open(vcfdir + '/Vcf_file_' + wildcards.chromosome + '.txt')
#         vcf_file = file.read().strip("\n")
#
#     return(vcf_file)


if config['ploidy'] == 1:
    rule rename_phased:
        input:
            vcf = f'{vcfOut}/{{chromosome}}_final.vcf.gz',
            idx = f'{vcfOut}/{{chromosome}}_final.vcf.gz.csi'
        output:
            vcf = f'{vcfOut}/{{chromosome}}_phased.vcf.gz',
            idx = f'{vcfOut}/{{chromosome}}_phased.vcf.gz.csi'
        log: 'logs/Rename_phased_{chromosome}.log'
        resources: cpus=1, mem_mb=32000, time_min=60
        shell:
            """
            if [ -h {input.vcf} ]; then
                ln -s $( realpath {input.vcf} ) {output.vcf}
                ln -s $( realpath {input.idx} ) {output.idx}
            else
                ln -s {input.vcf} {output.vcf}
                ln -s {input.idx} {output.idx}
            fi
            """

else:
    rule phase:
        input:
            vcf = f'{vcfOut}/{{chromosome}}_final.vcf.gz',
        output:
            file = f'{vcfOut}/{{chromosome}}_phased.vcf.gz',
            idx = f'{vcfOut}/{{chromosome}}_phased.vcf.gz.csi'
        params:
            map = f'{mapdir}/{{chromosome}}.txt',
        #log: 'logs/phase_{chromosome}.log'
        threads: 25
#        resources: cpus=20, mem_mb=25000, time_min=5
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
            bcftools index -f {output}
            """
