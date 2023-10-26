rule get_af:
    input: f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
    output:
        new_file = temp(f'{vcfdir}/{{chromosome}}_phased_info.vcf.gz'), 
        info = temp(f'{vcfdir}/Info{{chromosome}}.INFO')
    params:
        prefix=f'{vcfdir}/Info{{chromosome}}'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/get_af_{chromosome}.log'
    shell:
        """
        bcftools +fill-tags {input} -Oz -o {output.new_file} -- -t AN,AC,AF
        vcftools --gzvcf {output.new_file} --out {params.prefix} --get-INFO AC --get-INFO AF
	
	LOGFILE={params.prefix}.log
	if test -f "$LOGFILE"; then
	    rm $LOGFILE
	fi
        """

rule get_major:
    input: rules.get_af.output.info
    output: temp(f'{vcfdir}/Major{{chromosome}}.txt')
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/get_major_{chromosome}.log'
    shell:
        """
        awk '{{if (NR!=1 && $5>=0.5) {{print $1"_"$2","$4}} else if (NR!=1 && $5<0.5) {{print $1"_"$2","$3}}}}' {input} > {output}
        """

rule decompress:
    input:
        vcf = rules.get_af.output.new_file,
        #major = rules.get_major.output
    output: f'{vcfdir}/{{chromosome}}_phased_info.vcf' # this is removing both the .gz and the decompressed file 
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/decompress_{chromosome}.log'
    shell:
        """
        gunzip {input.vcf}
        """

rule extract_vcf_pos:
    input: rules.get_major.output
        #rules.decompress.output
    output:
        file = temp(f'{vcfdir}/VcfPos{{chromosome}}.txt'),
        #sites = temp('{chromosome}_sites.list')  
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/extract_vcf_pos_{chromosome}.log'
    shell:
        """
        cut -d',' -f 1 {input} > {output}
        """
    # why not just take the first col of the MajorAlle table??
    # # 
    # bcftools query -f '%CHROM %POS\n' {input} > {output.sites}
    #     awk '{{print $1"_"$2}}' {output.sites} > {output.file}

if config['ancestralAllele'] not null:
    ancestral_file = config['ancestralAllele']
else
   ancestral_file = "AncestralAllele/AncestralAllele_Vcf.txt"

rule match_ancestral_vcf:
    input:
        vcfPos = rules.extract_vcf_pos.output.file,
        ancestral = ancestral_file
        #ancestral = "AncestralAllele/AncestralAllele_Vcf.txt",
        major = rules.get_major.output
    output: 
        file = temp(f'{vcfdir}/AncestralVcfMatch{{chromosome}}.txt'),
    params:
        chrNum = lambda wc: wc.get('chromosome')[3:],
        ancestral_sites = f'{vcfdir}/{{chromosome}}.aa',
        tmp_file = f'{vcfdir}/{{chromosome}}.tmp',
    log: 'logs/match_ancestral_vcf_{chromosome}.log'
    shell:
        """
        grep "{params.chrNum}_" {input.ancestral} | grep -xv 'ambiguous' > {params.ancestral_sites}
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
        awk -F "," '{{print $1" "$2}}' {output} > {params.tmp_file} && mv {params.tmp_file} {output}
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
    input: rules.change_infoAA_vcf.output,
    output: 
        file = f'{vcfdir}/{{chromosome}}_ancestral.vcf.gz',
        idx = f'{vcfdir}/{{chromosome}}_ancestral.vcf.gz.csi'
    conda: "bcftools"
    threads: 1
    resources: cpus=1, mem_mb=64000, time_min=5
    log: 'logs/compress_{chromosome}.log'
    shell:
        """
        bgzip {input}
        bcftools index {output.file}
        """
# keeping the phased files w/o ancestral alleles 