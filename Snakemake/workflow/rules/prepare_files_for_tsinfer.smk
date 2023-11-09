if config['ancestralAllele'] is None:
    ancestral_file = "../Project/AncestralAllele/AncestralAllele_Vcf.txt"
else:
    ancestral_file = config['ancestralAllele']

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

    rule extract_ancestral_chromosome:
        input:
            ancestral_file
        output:
            temp(f'{vcfdir}/Ancestral{{chromosome}}.txt')
        params:
            chrNum = lambda wc: wc.get('chromosome')[3:]
        log: 'logs/extract_ancestral_chromosome_{chromosome}.log'
        shell:
            """
            grep "{params.chrNum}_" {input} | grep -xv 'ambiguous' > {output}
            """

    rule join_major_ancestral:
        input:
            ancestral = rules.extract_ancestral_chromosome.output,
            major = rules.get_major.output
        output:
            temp(f'{vcfdir}/Major_and_ancestral{{chromosome}}.txt')
        shell:
            "awk -F"," 'NR==FNR{A[$1]=$2;next}{print$0 FS (A[$1]?A[$1]:"0")}' {input.ancestral} {input.major} > {output}"

    rule determine_ancestral_major:
        input:
            rules.join_major_ancestral.output
        output:
            temp(f'{vcfdir}/Major_or_ancestral{{chromosome}}.txt')
        shell:
            "awk -F","  '{if ($3 == 0) {print $1" "$2} else {print $1" "$3}}' {input} > {output}"

else:
   ancestral_file = "../Project/AncestralAllele/AncestralAllele_Vcf.txt"

   #If this has been done by snakemake, this already includes eiher the ancestral or the Major_or_ancestral
   rule get_vcf_position:
       input:
            f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
        output:
            f'{vcfdir}/{{chromosome}}_position.txt'
        conda: "bcftools"
        shell:
            """
            bcftools query -f "%CHROM\_%POS\n" {input} {output}
            """

    rule sort_ancestral_to_vcf:
        input:
            aa=ancestral_file,
            vcfPos=rules.get_vcf_position.output
        output:
            temp(f'{vcfdir}/Major_or_ancestral{{chromosome}}.txt')
        shell:
            """
            awk 'FNR == NR {{ lineno[$1] = NR; next}} {{print lineno[$1], $0;}}' {input.vcfPos} {input.aa} | sort -k 1,1n | cut -d' ' -f2- > {output}
            """


rule decompress:
    input:
        f'{vcfdir}/{{chromosome}}_phased.vcf.gz'
    output: f'{vcfdir}/{{chromosome}}_phased_info.vcf' # this is removing both the .gz and the decompressed file
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    log: 'logs/decompress_{chromosome}.log'
    shell:
        """
        gunzip {input}
        """

rule change_infoAA_vcf:
    input:
        vcf=rules.decompress.output,
        ancestralAllele=f'{vcfdir}/Major_or_ancestral{{chromosome}}.txt'
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
