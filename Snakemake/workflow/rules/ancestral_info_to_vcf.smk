
# extract all positions in vcf file using bcftools and save to regions file
# for rows in ancestral_alleles file where col 1 match wildcard.chromosome
# for each row in regions file, if value is in ancestral_alleles col 2, add
# ancestral_alleles col 3 to regions file col 2
# else add -1

ancestral_file = config['ancestral_allele']

if "global_cattle" in config['PROJECT']:
    # rule list_samples:
    #     input:
    #         vcf_file = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz',
    #         metadata = config['meta']
    #     output:
    #         samples = f'{vcfdir}/{{chromosome}}/samples_to_remove'
    #     conda: 'bcftools'
    #     shell:
    #         r"""
    #         bcftools query -l {input.vcf_file} > {wildcards.chromosome}.samples
    #         cut -d',' -f1 {input.metadata} | tail -n +2 > metadata.samples
    #         grep -v -f metadata.samples vcf.samples > unique.samples
    #         """
    
    rule remove_samples:
        input: 
            vcf_file = input_ancestral,
            idx_file = input_ancestral_idx,
        output: f'{vcfdir}/{{chromosome}}/{{chromosome}}_analysis.vcf',
        params:
            samples_to_remove = f'{vcfdir}/samples.remove'
        conda: 'bcftools'
        threads: 32
        resources: cpus=32, mem_mb=2048000, time_min=1200
        benchmark: 'benchmarks/{chromosome}.remove.benchmark.txt' 
        shell:
            r"""
            bcftools view -S ^{params.samples_to_remove} {input.vcf_file} --threads {threads} -O v -o {output}
            """

else:
    rule decompress:
        input: 
            vcf_file = input_ancestral,
            idx_file = input_ancestral_idx
        output: f'{vcfdir}/{{chromosome}}/{{chromosome}}_analysis.vcf'
        conda: 'bcftools'
        benchmark: 'benchmarks/{chromosome}.decompress.benchmark.txt' 
        shell:
            r"""
            bcftools view {input.vcf_file} -O v -o {output}
            """

input_ancestral = f'{vcfdir}/{{chromosome}}/{{chromosome}}_analysis.vcf'     

rule get_vcf_positions:
    input: 
        vcf_file = input_ancestral,
    output: 
        regions_file = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}.regions')
    conda: 'bcftools'
    benchmark: 'benchmarks/{chromosome}.positions.benchmark.txt' 
    shell:
        r"""
        bcftools query -f '%POS\n' {input.vcf_file} | awk '{{print $1, -1}}' > {output.regions_file}
        """

rule get_ancestral_alleles:
    input: 
        ancestral_alleles_file = ancestral_file
    output:
        ancestrals = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}.anc')
    params:
        chrnum = lambda wc: wc.get('chromosome')[3:]
    benchmark: 'benchmarks/{chromosome}.getaa.benchmark.txt'     
    shell:
       r"""
       # filter positions in ancestral that correspond to chr
        awk '$1 == {params.chrnum}' {input.ancestral_alleles_file} | cut -d' ' -f2,3 > {output.ancestrals}
       """       

rule update_regions_file:
    input:
        regions_file = rules.get_vcf_positions.output.regions_file,
        ancestrals = rules.get_ancestral_alleles.output.ancestrals,
    output:
        updated_regions_file = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}.updated.regions'),
        positions_exclude = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}.exclude'),
    benchmark: 'benchmarks/{chromosome}.regions.benchmark.txt'     
    shell:
        r"""
        awk 'NR==FNR{{a[$1]=$2;next}} ($1 in a){{$2=a[$1]}}1' {input.ancestrals} {input.regions_file} > {output.updated_regions_file}

        awk '$2 == -1' {output.updated_regions_file} | cut -d' ' -f1 > {output.positions_exclude}
        """

rule add_AA_to_regions_file:
    input:
        updated_regions_file = rules.update_regions_file.output.updated_regions_file,
    output:
        new_regions_file = f'{vcfdir}/{{chromosome}}/{{chromosome}}.final.regions'
    benchmark: 'benchmarks/{chromosome}.addaaregions.benchmark.txt'     
    shell:
        r"""
        echo "INFO" > {output.new_regions_file}
        awk '{{print "AA="$2}}' {input.updated_regions_file} >> {output.new_regions_file}
        """

rule add_info_to_vcf:
    input: 
        vcf_file = input_ancestral,
        # f'{vcfdir}/{{chromosome}}/{{chromosome}}_correctids.vcf',
        aa_info = rules.add_AA_to_regions_file.output.new_regions_file,
    output: 
        large_vcf = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf')
    conda: 'bcftools'
    benchmark: 'benchmarks/{chromosome}.aainfo.benchmark.txt' 
    shell:
        r"""
        awk 'NR==FNR{{a[NR] = $1; next}} FNR<15{{print}}; \
        FNR==6{{printf "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"}}; \
        FNR>14{{$8=a[FNR-14]; print}}' OFS="\t" {input.aa_info} {input.vcf_file} > {output.large_vcf}

        rm {input.vcf_file}
        """

rule compress_index_vcf:
    input: 
        large_vcf = rules.add_info_to_vcf.output.large_vcf
        # f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf'
    output: 
        vcf_file = f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf.gz',
        vcf_indx = f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf.gz.csi'
    conda: 'bcftools'
    benchmark: 'benchmarks/{chromosome}.compress.benchmark.txt' 
    shell:
        r"""
        bgzip -c {input} > {output.vcf_file}
        bcftools index {output.vcf_file}
        """
input_ts = rules.compress_index_vcf.output.vcf_file
