if config['ancestral_allele'] is not None:
    ancestral_file = Path(vcfdir, config['ancestral_allele'])
else:
    ancestral_file = rules.combine_pos_ancestral.output


rule get_af:
    input: f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz'
    output:
        vcf = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_info.vcf.gz'),
        info = temp(f'{vcfdir}/{{chromosome}}/Info{{chromosome}}.INFO')
    params:
        prefix=f'{vcfdir}/{{chromosome}}/Info{{chromosome}}'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=32000, time_min=60
    log: 'logs/get_af_{chromosome}.log'
    shell:
        """
        bcftools +fill-tags {input} -Oz -o {output.vcf} -- -t AN,AC,AF
        vcftools --gzvcf {output.vcf} --out {params.prefix} --get-INFO AC --get-INFO AF

        LOGFILE={params.prefix}.log
        if test -f "$LOGFILE"; then
            rm $LOGFILE
        fi
        """

rule get_major:
    input: rules.get_af.output.info
    output: temp(f'{vcfdir}/{{chromosome}}/Major{{chromosome}}.txt')
    threads: 1
    resources: cpus=1, mem_mb=32000, time_min=60
    log: 'logs/Get_major_{chromosome}.log'
    shell:
        """
        awk '{{if (NR!=1 && $5>=0.5) {{print $1"_"$2","$4}} else if (NR!=1 && $5<0.5) {{print $1"_"$2","$3}}}}' {input} > {output}
        """

rule extract_ancestral_chromosome:
    input:
        ancestral_file
    output:
        temp(f'{vcfdir}/{{chromosome}}/Ancestral{{chromosome}}.txt')
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
        temp(f'{vcfdir}/{{chromosome}}/MajAAnc_{{chromosome}}.txt')
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
        temp(f'{vcfdir}/{{chromosome}}/MajOAnc_{{chromosome}}.txt')
    log: 'logs/Determine_major_ancestral_{chromosome}.log'
    shell:
        """
        awk -F","  '{{if ($3 == 0) {{print $1" "$2}} else {{print $1" "$3}}}}' {input} > {output}
        """

rule decompress:
    input:
        vcf = rules.get_af.output.vcf,
	    majororaa = rules.determine_ancestral_major.output
    output: f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf' 
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
        vcf = rules.decompress.output,
        ancestralAllele = rules.determine_ancestral_major.output
    output: f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf'
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
    input: 
        vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf',
    output:
        vcf = f'{vcfdir}/{{chromosome}}_ancestral.vcf.gz',
        idx = f'{vcfdir}/{{chromosome}}_ancestral.vcf.gz.csi',
    params:
        chrdir = f'{vcfdir}/{{chromosome}}/',
        outdir = f'{vcfdir}/'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=64000, time_min=60
    log: 'logs/Compress_{chromosome}.log'
    shell:
       r"""
        bgzip {input.vcf}
        bcftools index {input.vcf}.gz

        mv {input}.gz {input}.gz.csi {params.outdir}
        rm -r {params.chrdir}
        """
