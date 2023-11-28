if config['ancestral_allele'] is not None:
    ancestral_file = Path(vcfdir, config['ancestral_allele'])
else:
    ancestral_file = rules.combine_pos_ancestral.output


rule get_af:
    input: f'{vcfOut}/{{chromosome}}_phased.vcf.gz'
    output:
        info = temp(f'{vcfOut}/Info{{chromosome}}.INFO')
    params:
        prefix=f'{vcfOut}/Info{{chromosome}}'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=32000, time_min=60
    log: 'logs/get_af_{chromosome}.log'
    shell:
        """
        bcftools +fill-tags {input} -Oz -o {input} -- -t AN,AC,AF
        vcftools --gzvcf {input} --out {params.prefix} --get-INFO AC --get-INFO AF

        LOGFILE={params.prefix}.log
        if test -f "$LOGFILE"; then
            rm $LOGFILE
        fi
        """

rule get_major:
    input: rules.get_af.output.info
    output: temp(f'{vcfOut}/Major{{chromosome}}.txt')
    threads: 1
    resources: cpus=1, mem_mb=32000, time_min=60
    log: 'logs/Get_major_{chromosome}.log'
#    wildcard_constraints:
#        chromosome="^chr | ^Chr"
    shell:
        """
        awk '{{if (NR!=1 && $5>=0.5) {{print $1"_"$2","$4}} else if (NR!=1 && $5<0.5) {{print $1"_"$2","$3}}}}' {input} > {output}
        """

rule extract_ancestral_chromosome:
    input:
        ancestral_file
    output:
        temp(f'{vcfOut}/Ancestral{{chromosome}}.txt')
    params:
        chrNum = lambda wc: wc.get('chromosome')[3:]
    log: 'logs/Extract_ancestral_chromosome_{chromosome}.log'
    resources: cpus=1, mem_mb=32000, time_min=30
    shell:
        """
        grep "{params.chrNum}_" {input} | grep -xv 'ambiguous' > {output}
        """

rule join_major_ancestral:
    input:
        ancestral = rules.extract_ancestral_chromosome.output,
        major = rules.get_major.output
    output:
        temp(f'{vcfOut}/MajAAnc_{{chromosome}}.txt')
    log: 'logs/Join_major_ancestral_{chromosome}.log'
    resources: cpus=1, mem_mb=32000, time_min=60
    shell:
        """
        awk -F"," 'NR==FNR{{A[$1]=$2;next}}{{print$0 FS (A[$1]?A[$1]:"0")}}' {input.ancestral} {input.major} > {output}
        """

rule determine_ancestral_major:
    input:
        rules.join_major_ancestral.output
    output:
        temp(f'{vcfOut}/MajOAnc_{{chromosome}}.txt')
    log: 'logs/Determine_major_ancestral_{chromosome}.log'
    shell:
        """
        awk -F","  '{{if ($3 == 0) {{print $1" "$2}} else {{print $1" "$3}}}}' {input} > {output}
        """

rule decompress:
    input:
        vcf=f'{vcfOut}/{{chromosome}}_phased.vcf.gz',
	majororaa=rules.determine_ancestral_major.output
    output: f'{vcfOut}/{{chromosome}}_phased.vcf' # this is removing both the .gz and the decompressed file
    threads: 1
    resources: cpus=1, mem_mb=32000, time_min=30
    log: 'logs/Decompress_{chromosome}.log'
    conda: "bcftools"
    shell:
        """
        bcftools view {input.vcf} -O v -o {output}
        """

rule change_infoAA_vcf:
    input:
        vcf=rules.decompress.output,
        ancestralAllele=rules.determine_ancestral_major.output
    output: f'{vcfOut}/{{chromosome}}_ancestral.vcf'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=64000, time_min=120
    log: 'logs/Change_infoAA_vcf_{chromosome}.log'
    shell:
        """
        HEADERNUM="$(( $(bcftools view -h {input.vcf} | wc -l) - 1 ))"
        INFOLINE=$(( $(bcftools view -h {input.vcf} | awk '/INFO/{{print NR}}' | head -n 1) ))
        awk -v OFS="\t" -v HEADER=$HEADERNUM -v INFO=$INFOLINE 'NR==FNR{{{{a[FNR] = $2; next}}}} FNR<=HEADER{{{{print}}}}; \
        FNR==INFO{{{{printf "##INFO=<ID=AA,Number=1,Type=String,Description=Ancestral Allele>\\n"}}}}; \
        FNR>HEADER{{{{$8="AA="a[FNR-HEADER]; print}}}}' OFS="\t" {input.ancestralAllele} {input.vcf} > {output}
        """

rule compress_vcf:
    input: rules.change_infoAA_vcf.output,
    output:
        file = f'{vcfOut}/{{chromosome}}_ancestral.vcf.gz',
        idx = f'{vcfOut}/{{chromosome}}_ancestral.vcf.gz.csi'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=64000, time_min=60
    log: 'logs/Compress_{chromosome}.log'
    shell:
        """
        bgzip {input}
        bcftools index -f {output.file}
        """
# keeping the phased files w/o ancestral alleles
