

rule get_af:
    input:
        f'{vcfdir}/{{chromosome}}_final.vcf.gz'
    output:
        info=temp('Info{chromosome}.INFO'),
        log=temp('Info{chromosome}.log')
    params:
        prefix='Info{chromosome}'
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

# I think we should skip this - depending on the samples it can get very complicated
# rule combine_major_ancestral:
#     input:
#         ancestral=config['ancestralAllele'],
#         major=rules.get_major.output
#     output: temp('AncestralMajor{chromosome}.txt')
#     shell:
#         """
#         join -a1 -t ","  -j 1 -o 1.1,1.2,2.2 <(sort -t"," -k1,1 --version-sort {input.major}) <(sort -t"," -k1,1 --version-sort {input.ancestral}) > tmpMA
#         awk -F, '{{if ($3=="") {{print $1,$2}} else {{print $1,$3}}}}' tmpMA > {output}
#         rm tmpMA
#         """

rule decompress:
    input:
        vcf = f'{vcfdir}/{{chromosome}}_final.vcf.gz',
        major=rules.combine_major_ancestral.output
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
        major=rules.get_major.output,
#        ancestralMajor=rules.combine_major_ancestral.output
        # The ancestral file has to have chr_pos and AA, split with a tab
    output: temp('AncestralVcfMatch{chromosome}.txt')
        # THis files needs to contain all the variants from the vcf (blank space)
    shell:
    # take only the ancestrals referent to the current chr
    # read the position file (1) line by line;
    # test if the position existis in the ancestral file (2),
    # if TRUE print the line from file (2),
    # else print the line from file (3)
        """
        grep "{wildcards.chromosome}_" {input.ancestral} > aa_tmp

        for line in $(cat {input.vcfPos});
        do
            if grep -w "${line}" aa_tmp >> {output}; then
                continue
            else
                grep -w "${line}" {input.major} >> {output}
            fi
        done
#        awk -F"," '{{print $1"\t"$2}}' tmp > {output}
        """

rule change_infoAA_vcf:
# changed OFS for the pos/anc file to "," instead of "\t"
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
        FNR>HEADER{{{{$8="AA="a[FNR-HEADER]; print}}}}' OFS="," {input.ancestralAllele} {input.vcf} > {output}
        """

rule compress_vcf:
    input:
        rules.change_infoAA_vcf.output
    output: f'{vcfdir}/{{chromosome}}_ancestral.vcf'
    shell:
        """
        bgzip {input}
        bcftools index {output}
        """
