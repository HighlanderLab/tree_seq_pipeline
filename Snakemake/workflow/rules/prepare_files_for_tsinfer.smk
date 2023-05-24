

rule get_af:
    input:
        f'{vcfdir}/{{chromosome}}_final.vcf.gz'
    output:
        info=temp('Info{chromosome}.INFO'),
        log=temp('Info{chromosome}.log')
    params:
        prefix='Info{chromosome}'
    envmodules:
        config['bcftoolsModule'],
        config['vcftoolsModule']
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    shell:
        """
        bcftools +fill-tags {input} -Oz -o {input} -- -t AN,AC,AF
        vcftools --gzvcf {input} --out {params.prefix} --get-INFO AC --get-INFO AF
        """

rule get_major:
    input:
        rules.get_af.output.info
    output:
        temp('Major{chromosome}.txt')
    shell:
        """
        awk '{{if (NR!=1 && $5>=0.5) {{print $1"_"$2","$4}} else if (NR!=1 && $5<0.5) {{print $1"_"$2","$3}}}}' {input} > {output}
        """

rule decompress:
    input:
        vcf = f'{vcfdir}/{{chromosome}}_final.vcf.gz',
        major=rules.get_major.output
    output: temp(f'{vcfdir}/{{chromosome}}_final.vcf')
    shell:
        """
        gunzip {input.vcf}
        """

rule extract_vcf_pos:
    input:
        rules.decompress.output
    output:
        temp('VcfPos{chromosome}.txt')
    #conda:
    #    config['tsinferEnv']
    #envmodules:
    #    config['bcftoolsModule']
    params:
        vcfDir=config['vcfDir']
    envmodules:
        config['bcftoolsModule']
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input} > tmp
        awk '{{print $1"_"$2}}' tmp > {output}
        rm tmp
        """

# I changed here so the rule now takes as input
#   (1) positions in the vcf
#   (2) ancestral allele
#   (3) major allele
rule match_ancestral_vcf:
    input:
        vcfPos=rules.extract_vcf_pos.output, # This has more lines,
        ancestral=config['ancestralAllele'],
        major=rules.get_major.output
#        ancestralMajor=rules.combine_major_ancestral.output
        # The ancestral file has to have chr_pos and AA, split with a tab
    output: temp('AncestralVcfMatch{chromosome}.txt')
        # THis files needs to contain all the variants from the vcf (blank space)
    params:
        chrNum=lambda wc: wc.get("chromosome")[3:]
    shell:
    # take only the ancestrals referent to the current chr
    # read the position file (1) line by line;
    # test if the position existis in the ancestral file (2),
    # if TRUE print the line from file (2),
    # else print the line from file (3)
    # at the end remove POS column
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
# this INFOLINE=$(( $(grep -Fn "INFO" {input.vcf} | cut -d ":" -f 1  | head -n1) ))
# takes forever! quicker having bcftools to print the header and then grep INFO
# same here HEADERNUM=$(( $(grep "##" {input.vcf} | wc -l) + 1 ))
    input:
        vcf=rules.decompress.output,
        ancestralAllele=rules.match_ancestral_vcf.output
    output: f'{vcfdir}/{{chromosome}}_ancestral.vcf'
    envmodules:
        config['bcftoolsModule']
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
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
    envmodules:
        config['bcftoolsModule']
    conda: "HLab_tsinfer"
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=5
    shell:
        """
        bgzip {input}
        bcftools index {output}
        """
