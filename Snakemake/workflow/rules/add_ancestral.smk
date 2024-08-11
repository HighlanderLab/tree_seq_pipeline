
ancestral_file = config['ancestral_allele']

# input and output
infile = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz'
infileidx = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz.csi'
outfile = f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf.gz'
outidx = f'{vcfdir}/{{chromosome}}_ancestral.vcf.gz.csi'

# temporary files
filtered_samples = f'{vcfdir}/{{chromosome}}/{{chromosome}}_sam.vcf.gz'
filtered_positions = f'{vcfdir}/{{chromosome}}/{{chromosome}}_sampos.vcf'
ancestral_not_compressed = f'{vcfdir}/{{chromosome}}/{{chromosome}}_ancestral.vcf'
regfile = f'{vcfdir}/{{chromosome}}_regions.tsv'
aaInfo = f'{vcfdir}/{{chromosome}}_aa.info'

if config['PROJECT'] == "global_cattle":
    rule remove_samples:
        input: infile
        output: filtered_samples
        params:
            samples_to_remove = f'{vcfdir}/{{chromosome}}/remove'
        conda: 'bcftools'
        threads: 32
        resources: cpus=32, mem_mb=2048000, time_min=1200
        shell:
            r"""
            bcftools view -S ^{params.samples_to_remove} {input} --threads {threads} -O z -o {output}

            bcftools index {output}
            """


# need to some helper files:
# 1) create a chr pos file for filtering
# 2) create an info file for adding ancestrals
# filter vcf positions so that it has the same as the aa file
# add ancestral allele to INFO

rule make_regions_file:
    input: 
        anc = ancestral_file,
        vcf = filtered_samples
    output: temp(regfile)
    params:
        chrnum = lambda wc: wc.get('chromosome')[3:]
    shell:
        r"""
        # filter positions in ancestral that correspond to chr
        cut -d' ' -f1,2 {input.ancestral_file} | awk '$1 == {params.chrnum}' | sed -e $'s/ /\t/g' > chr{params.chrnum}.anc

        cut -f2 chr{params.chrnum}.anc > chr{params.chrnum}.compare

        # get positions in vcf
        bcftools query -l '%POS\n' {input.vcf} > chr{params.chrnum}.file

        # match positions ancestral and vcf
        comm -12 <( sort chr{params.chrnum}.file ) <( chr{params.chrnum}.compare) > chr{params.chrnum}.keep

        # number of lines
        l=($(wc -l chr{params.chrnum}.keep))

        # create final regions file
        yes {params.chrnum} | head -n $l > chr{params.chrnum}
        paste chr{params.chrnum} chr{params.chrnum}.keep > {output}

        rm chr{params.chrnum} chr{params.chrnum}.keep chr{params.chrnum}.file chr{params.chrnum}.compare chr{params.chrnum}.anc
        """
#grep -vw 4579381 3565339

rule filter_positions:
    input: 
        vcf = filtered_samples,
        idx = f'{filtered_samples}.csi',
        regions = regfile,
    output: filtered_positions
    conda: 'bcftools'
    threads: 32
    resources: cpus = 32, mem_mb = 768000, time_min = 300
    shell:
        r"""
        bcftools view -R {input.regions} {input.vcf} \
        #            --threads {threads} \
                     -o {output}
        """

# take only the rows where col1 == current_chromosome
# add AA= in front of the aa col
# use only col pos & aa
# match with the left pos in filtered vcf

rule make_aaInfo_file:
    input: 
        aa_file = ancestral_file,
        vcf = filtered_positions
    output: temp(aaInfo)
    params:
        chrnum = lambda wc: wc.get('chromosome')[3:]
    shell:
        r"""
        # filter by chromosome and add AA= to the allele col
        awk '$1 == {params.chrnum}' {input} | awk '$3="AA="$3' | cut -d' ' -f2,3 > {wildcards.chromosome}.aa 

        # get the positions in the vcf
        bcftools query -f '%POS\n' {input.vcf} > {wildcards.chromosome}.pos 

        # add allele col to the pos in the vcf
        awk 'NR==FNR{{a[$1];next}}($1 in a){{print $0}}' {wildcards.chromosome}.pos {wildcards.chromosome}.aa > {wildcards.chromosome}.pos.aa

        cut -d' ' -f2  {wildcards.chromosome}.pos.aa | sed '1s/^/INFO\n/' > {output}

        rm {wildcards.chromosome}.aa {wildcards.chromosome}.pos {wildcards.chromosome}.pos.aa
        """

# need to run shapeit first to know line numbers
rule add_ancestral:
    input: 
        vcf = filtered_positions,
        aa_info = aaInfo,
    output: ancestral_not_compressed
    conda: 'bcftools'
    shell:
        r"""
        awk 'NR==FNR{{a[NR] = $1; next}} FNR<13{{print}}; \
        FNR==6{{printf "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"}}; \
        FNR>12{{$8=a[FNR-11]; print}}' OFS="\t" {input.aa_info} {input.vcf} > {output}
        """


# HEADERNUM="$(( $(bcftools view -h {input.vcf} | wc -l) - 1 ))"

#         INFOLINE=$(( $(bcftools view -h {input.vcf} | awk '/INFO/{{print NR}}' | head -n 1) ))

#         awk -v OFS="\t" -v HEADER=$HEADERNUM -v INFO=$INFOLINE 'NR==FNR{{{{a[FNR] = $2; next}}}} FNR<=HEADER{{{{print}}}}; \
#         FNR==INFO{{{{printf "##INFO=<ID=AA,Number=1,Type=String,Description=Ancestral Allele>\\n"}}}}; \
#         FNR>HEADER{{{{$8=a[FNR-HEADER]; print}}}}' OFS="\t" {input.aa_info} {input.vcf} > {output}
#  All_AAInfo.txt AllFiltered.recode.vcf > AllFiltered.recode.Ancestral.vcf

# compress and index final vcf
# rule match_samples:
#     input: ancestral_not_compressed
#     output: outfile
#     params: 
#         samples_to_keep = f'{vcfdir}/{{chromosome}}/keep.samples'
#     conda: 'bcftools'
#     shell:
#         r"""
#         bcftools view -S {params.samples_to_keep} {input} -O z -o {output}
#         """

rule compress_vcf:
    input: ancestral_not_compressed
    output: outfile,
    conda: 'bcftools'
    shell:
        r"""
        bgzip -c {input} > {output}
        """

rule index_output:
    input: outfile
    output: outidx
    conda: 'bcftools'
    shell:
        r"""
        bcftools index {input}
        """
