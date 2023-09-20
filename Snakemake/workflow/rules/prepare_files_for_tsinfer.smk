

rule get_af:
    input: f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
    output:
        #new_file = temp(f'{vcfdir}/{{chromosome}}_phased_info_test.vcf.gz'), 
        info = temp(f'{vcfdir}/Info{{chromosome}}.INFO')
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
    input: rules.get_af.output.info
    output: temp('Major{chromosome}.txt')
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/get_major_{chromosome}.log'
    shell:
        """
        awk '{{if (NR!=1 && $5>=0.5) {{print $1"_"$2","$4}} else if (NR!=1 && $5<0.5) {{print $1"_"$2","$3}}}}' {input} > {output}
        """

rule decompress:
    input:
        vcf = f'{vcfdir}/{{chromosome}}_phased.vcf.gz',
        major = rules.get_major.output
    output: f'{vcfdir}/{{chromosome}}_phased.vcf' # this is removing both the .gz and the decompressed file 
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
        file = temp('VcfPos{chromosome}.txt'),
        sites = temp('{chromosome}_sites.list')
    params:
        vcfDir = config['vcfDir'],   
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/extract_vcf_pos_{chromosome}.log'
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input} > {output.sites}
        awk '{{print $1"_"$2}}' {output.sites} > {output.file}
        """
    # why not just take the first col of the table above??
    # cut -d',' -f 1 {rules.get_major.output}

rule match_ancestral_vcf:
    input:
        vcfPos = rules.extract_vcf_pos.output.file,
        ancestral = config['ancestralAllele'],
        major = rules.get_major.output
    output: 
        file = temp('AncestralVcfMatch{chromosome}.txt'),
    params:
        chrNum = lambda wc: wc.get("chromosome")[3:],
        ancestral_sites = '{chromosome}.aa'
    threads: 1
    resources: cpus=1, mem_mb=150000, time_min=5
    log: 'logs/match_ancestral_vcf_{chromosome}.log'
    shell:
        """
        grep "{params.chrNum}_" {input.ancestral} | grep -xv 'ambigous' > {params.ancestral_sites}
        echo done

        for line in $(cat {input.vcfPos});
        do
            if grep -w "${{line}}" {params.ancestral_sites} >> {output.file}; then
                echo line in aa
                continue
            else
                grep -w "${{line}}" {input.major} >> {output.file}
                echo line in major
            fi
        done
        rm {params.ancestral_sites}
        awk -F "," '{{print $1" "$2}}' {output} > {wildcards.chromosome}.tmp && mv {wildcards.chromosome}.tmp {output}
        """
#sed -i 's/,/ /g' {output.file}

rule change_infoAA_vcf:
    input:
        vcf=rules.decompress.output,
        ancestralAllele=rules.match_ancestral_vcf.output
    output: f'{vcfdir}/{{chromosome}}_ancestral.vcf'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=100000, time_min=5
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
        vcf_aa = rules.change_infoAA_vcf.output,
        vcf_phased = rules.decompress.output
    output: 
        vcf_aa = f'{vcfdir}/{{chromosome}}_ancestral.vcf.gz',
        vcf_phased = f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=64000, time_min=5
    log: 'logs/compress_{chromosome}.log'
    shell:
        """
        bgzip {input.vcf_aa}
        bcftools index {output.vcf_aa}

        bgzip {input.vcf_phased}
        bcftool index {output.vcf_phased}
        """
# keeping the phased files w/o ancestral alleles 