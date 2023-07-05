

rule get_af:
    input:
        f'{vcfdir}/{{chromosome}}_final.vcf.gz'
    output:
        info=temp(f'{vcfdir}/Info{{chromosome}}.INFO')
    params:
        prefix=f'{vcfdir}/Info{{chromosome}}'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
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
    input:
        rules.get_af.output.info
    output:
        temp('Major{chromosome}.txt')
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/get_major_{chromosome}.log'
    shell:
        """
        awk '{{if (NR!=1 && $5>=0.5) {{print $1"_"$2","$4}} else if (NR!=1 && $5<0.5) {{print $1"_"$2","$3}}}}' {input} > {output}
        """

rule decompress:
    input:
        vcf = f'{vcfdir}/{{chromosome}}_final.vcf.gz',
        major=rules.get_major.output
    output: temp(f'{vcfdir}/{{chromosome}}_final.vcf')
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/decompress_{chromosome}.log'
    shell:
        """
        gunzip {input.vcf}
        """

rule extract_vcf_pos:
    input:
        rules.decompress.output
    output:
        temp('VcfPos{chromosome}.txt')
    params:
        vcfDir=config['vcfDir']
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/extract_vcf_pos_{chromosome}.log'
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input} > tmp
        awk '{{print $1"_"$2}}' tmp > {output}
        rm tmp
        """

rule match_ancestral_vcf:
    input:
        vcfPos=rules.extract_vcf_pos.output,
        ancestral=config['ancestralAllele'],
        major=rules.get_major.output
    output: temp('AncestralVcfMatch{chromosome}.txt')
    params:
        chrNum=lambda wc: wc.get("chromosome")[3:]
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/match_ancestral_vcf_{chromosome}.log'
    shell:
        """
        grep "{params.chrNum}_" {input.ancestral} > aa_tmp

        for line in $(cat {input.vcfPos});
        do
            if grep -w "${{line}}" aa_tmp >> {output}; then
                continue
            else
                grep -w "${{line}}" {input.major} >> {output}
            fi
        done
	awk -F "," '{{print $1" "$2}}' {output} > tmp && mv tmp {output}
	rm aa_tmp
        """

rule change_infoAA_vcf:
    input:
        vcf=rules.decompress.output,
        ancestralAllele=rules.match_ancestral_vcf.output
    output: f'{vcfdir}/{{chromosome}}_ancestral.vcf'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/change_infoAA_vcf_{chromosome}.log'
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
        rules.change_infoAA_vcf.output
    output: f'{vcfdir}/{{chromosome}}_ancestral.vcf'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/compress_{chromosome}.log'
    shell:
        """
        bgzip {input}
        bcftools index {output}
        """
