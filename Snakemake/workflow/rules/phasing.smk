# Takes the final vcf and the recombination map as inputs and run shapeit4.
# using a specific environment file so installs shapeit4 if needed.
# Outputs the phased vcf.gz with index for all chromosomes.
rule phase:
    input:
        vcf = rules.compress_vcf.output,
        map = '../mapsDir/{chromosome}.gmap'

    output: f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
#    conda:
#        'env/shapeit.yaml'
    log: 'logs/{chromosome}_phased.log'
    envmodules:
        config['bcftoolsModule']
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    shell:
        """
        str='{wildcards.chromosome}'
        chr=$(echo ${{str:4}})
        shapeit4 --input {input.vcf} \
                         --map {input.map} \
                         --region ${{chr}} \
                         --output {output} \
                         --sequencing \
                         --thread 10
        bcftools index {output}
        """
