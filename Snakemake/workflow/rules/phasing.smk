# Phase files if species is diploid using the output from the
# remove_missing_indels rule

file_to_phase = f'{vcfdir}/{{chromosome}}/{{chromosome}}.recode.vcf.gz'
file_to_phase_idx = f'{vcfdir}/{{chromosome}}/{{chromosome}}.recode.vcf.gz.csi'

phased_output = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz'
phased_output_idx = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz.csi'

checkpoint decide_start_point:
    output: "workflow_start.txt"
    shell:
        """
        if [ {config['phase_only']} == False ] && [ -f {phased_output} {phased_output_idx} ]; then
            echo "ancestral_info_to_vcf" > {output}
        else
            echo "start" > {output}
        fi
        """

if config['ploidy'] == 1:
    input_ancestral = file_to_phase
    input_ancestral_idx = file_to_phase_idx

else:
    rule phase:
        input:
            vcf = file_to_phase,
            idx = file_to_phase_idx,
        output:
            vcf = phased_output,
        params:
            map = f'{mapdir}/{{chromosome}}.txt',
        threads: 32
        resources: cpus = 32, mem_mb = 2048, time_min = 7200
        conda: 'shapeit4am'
        benchmark: 'benchmarks/{chromosome}.shapeit.benchmark.txt'
        shell:
            r"""
            str='{wildcards.chromosome}'
            chr=$(echo ${{str:3}})

            shapeit4 --input {input.vcf} \
                                --map {params.map} \
                                --region $chr \
                                --output {output.vcf} \
                                --thread {threads} \
                                --pbwt-depth 8 \
                                --sequencing \
                                --log {wildcards.chromosome}_phased.log                
            """         

    rule index_phased_cleanup: 
        input: rules.phase.output
        output: phased_output_idx
        params: 
            chr = lambda wildcards: wildcards.chromosome[3:],
            dire = f'{vcfdir}/{{chromosome}}',    
        conda: 'bcftools'
        benchmark: 'benchmarks/{chromosome}.bcftools_index.bechmark.txt'
        shell:
            r"""
            bcftools index {input}  
            """
    input_ancestral = phased_output
    input_ancestral_idx = phased_output_idx