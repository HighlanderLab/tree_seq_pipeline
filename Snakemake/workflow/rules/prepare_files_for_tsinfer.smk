rule get_af:
    input: f'{vcfdir}/{{chromosome}}.vcf.gz'
    output: temp('Info{chromosome}.INFO')
    params:
        prefix='Info{chromosome}'
    shell:
        """
        bcftools +fill-tags {input} -Oz -o {input} -- -t AN,AC,AF
        vcftools --gzvcf {input} --out {params.prefix} --get-INFO AC --get-INFO AF
        """

rule get_major:
    input:
        rules.get_af.output
    output: temp('Major{chromosome}.txt')
    shell:
        "awk '{{if (NR!=1 && $5>=0.5) {{print $1"_"$2","$4}} else"
        "if (NR!=1 && $5<0.5) {{print $1"_"$2","$3}}}}' {input} > {output}"

rule combine_major_ancestral:
    input:
        ancestral=config['ancestralAllele'],
        major=rules.get_major.output
    output: tmp('AncestralMajor{chromosome}.txt')
    shell:
        """
        join -a1 -t ","  -j 1 -o 1.1,1.2,2.2 <(sort -k1,1 {input.major}) <(sort -k1,1 {input.ancestral}) > tmpMA
        awk -F, '{{if ($3=="") {{print $1,$2}} else {{print $1,$3}}}}' tmpMA > {output}
        """

rule decompress:
    input:
        vcf = f'{vcfdir}/{{chromosome}}.vcf.gz',
        major=rules.combine_major_ancestral.output
    output: temp(f'{vcfdir}/{{chromosome}}.vcf')
    shell:
        """
        gunzip {input.vcf}
        """

rule extract_vcf_pos:
    input:
        rules.decompress.output
    output: temp('VcfPos{chromosome}.txt')
    #conda:
    #    config['tsinferEnv']
    #envmodules:
    #    config['bcftoolsModule']
    params:
        vcfDir=vcfDir
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input} > tmp
        awk '{{print $1"_"$2}}' tmp > {output}
        """

rule match_ancestral_vcf:
    input:
        vcfPos=rules.extract_vcf_pos.output, # This has more lines,
        ancestralMajor=rules.combine_major_ancestral.output
        # The ancestral file has to have chr_pos and AA, split with a tab
    output: temp('AncestralVcfMatch{chromosome}.txt')
        # THis files needs to contain all the variants from the vcf (blank space)
    shell:
        """
        for line in $(cat {input.vcfPos});
        do
          grep $line {input.ancestralMajor} || echo "";
        done > tmp
        awk -F"," '{{print $1"\t"$2}}' tmp > {output}
        """

rule change_infoAA_vcf:
    input:
        vcf=rules.decompress.output,
        ancestralAllele=rules.match_ancestral_vcf.output
    output: temp(f'{vcfdir}/{{chromosome}}_ancestral.vcf')
    shell:
        """
        HEADERNUM=$(( $(grep "##" {input.vcf} | wc -l) + 1 ))
        INFOLINE=$(( $(grep -Fn "INFO" {input.vcf} | cut -d ":" -f 1  | head -n1) ))
        awk -v HEADER=$HEADERNUM -v INFO=$INFOLINE 'NR==FNR{{{{a[FNR] = $2; next}}}} FNR<=HEADER{{{{print}}}}; \
        FNR==INFO{{{{printf "##INFO=<ID=AA,Number=1,Type=String,Description=Ancestral Allele>\\n"}}}}; \
        FNR>HEADER{{{{$8="AA="a[FNR-HEADER]; print}}}}' OFS="\t" {input.ancestralAllele} {input.vcf} > {output}
        """

rule compress_vcf:
    input:
        rules.change_infoAA_vcf.output
    output: 'Tsinfer/{chromosome}_ancestral.vcf.gz'
    shell:
        """
        bgzip {input}
        bcftools index {output}
        """
