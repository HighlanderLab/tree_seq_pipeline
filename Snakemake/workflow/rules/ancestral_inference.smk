from __future__ import division
import sys

if config['ancestral_allele'] == None:
    rule extract_aligned_alleles_from_vcf:
        input:
            alignedFocal=Path(vcfdir, config['aligned_focal']),
            vcf=Path(vcfdir, config['raw_vcf'])
        output:
            info=f"{oDir}/AncestralAllele/Raw_vcf_info.INFO", 
            log=temp(f"{oDir}/AncestralAllele/Raw_vcf_info.log")
        params:
            prefix=f"{oDir}/AncestralAllele/Raw_vcf_info"
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=64000, time_min=300
        log: 'logs/Extract_aligned_alleles_from_VCF.log'
        shell:
            "scripts/Extract_vcf_pos.sh {input.vcf} {input.alignedFocal} {params.prefix}"


    rule extract_snps_from_info:
        input:
            f"{oDir}/AncestralAllele/Raw_vcf_info.INFO"
        output:
            f"{oDir}/AncestralAllele/Raw_vcf_info_SNPs.INFO"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=60
        log: 'logs/Extract_SNPs_from_info.log'
        shell:
            """
            head -n1 {input} > Header
            sed -i "s/\t/ /g" Header
            ./scripts/Extract_SNPs_from_info.awk {input} > tmp
            cat Header tmp > {output}
            rm Header tmp
            """

    rule extract_pos_aligned_alleles_from_vcf:
        input:
            rules.extract_snps_from_info.output
        output:
            f"{oDir}/AncestralAllele/Raw_vcf_full_pos.txt"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=60
        log: 'logs/Extract_pos_aligned_alleles_from_vcf.log'
        shell:
            """
            awk '{{print $1" "$2}}' {input} > {output}
            """

    rule extract_vcf_alleles_from_aligned:
        input:
            alignedFocal=Path(vcfdir, config['aligned_focal']),
            vcfPos=rules.extract_pos_aligned_alleles_from_vcf.output
        output:
            f"{oDir}/AncestralAllele/Aligned_snps_focal_sorted_raw_vcf.txt"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=120
        log: 'logs/Extract_vcf_aleles_from_aligned.log'
        shell:
            """
            grep -Fwf {input.vcfPos} {input.alignedFocal} > {output}
            """

    rule create_estsfs_dicts:
        input:
            alignedAlleles=rules.extract_vcf_alleles_from_aligned.output,
            vcfAlleles=rules.extract_snps_from_info.output
        conda: "HLab_tsinfer"
        output:
           [f"{oDir}/AncestralAllele/Estsfs/" + x for x in expand("EstSfs_Dict{{chunk}}.csv", chunk = range(config['no_estsfs_chunks']))]
        wildcard_constraints:
            chunk="\d+"
        params:
            outDir=f"{oDir}/AncestralAllele/Estsfs",
            noCycle=config['no_estsfs_chunks']
        threads: config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=64000, time_min=360
        log: 'logs/Create_estsfs_dicts{chunk}.log'
        shell:
            "python scripts/Create_estsfs_input_fromWGAbed_Cactus.py {wildcards.chunk} {params.noCycle} {input.alignedAlleles} {input.vcfAlleles} {params.outDir}"
            #CreateInputForEstsfs_Loop.sh This is qsub

    rule extract_estsfs_pos:
        input:
            rules.create_estsfs_dicts.output
        output:
            f"{oDir}/AncestralAllele/Estsfs/EstSfs_Pos{{chunk}}.csv"
        params:
            chunk=config['no_estsfs_chunks']
        wildcard_constraints:
            chunk="\d+"
        threads: config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/Extract_estsfs_dicts{chunk}_pos.log'
        shell:
            """
            cut -f1 {input} > {output}
            """

    rule edit_estsfs_dicts:
        input:
            rules.create_estsfs_dicts.output
        output:
            tmp1=temp(f"{oDir}/AncestralAllele/Estsfs/tmp1{{chunk}}.csv"),
            tmp2=temp(f"{oDir}/AncestralAllele/Estsfs/tmp2{{chunk}}.csv"),
            editedDict=f"{oDir}/AncestralAllele/Estsfs/EstSfs_Dict{{chunk}}E.csv"
        params:
            chunk=config['no_estsfs_chunks']
        wildcard_constraints:
            chunk="\d+"
        threads: config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/Edit_estsfs_dicts{chunk}.log'
        shell:
            """
            cut -f2,3,4,5 {input} | grep -v "()" > {output.tmp1}
            echo $PWD
            sed -i "s/ //g" {output.tmp1}
            # Set the correct separators for the file
            awk -F "\t" '{{print $1"\t"$2" "$3" "$4}}' {output.tmp1} > {output.tmp2}
            # Remove the parenthesis from the file
            sed -i "s/(//g" {output.tmp2}
            sed "s/)/ /g" {output.tmp2} > {output.editedDict}
            """

    rule combine_estsfs_pos:
        input:
            expand(f"{oDir}/AncestralAllele/Estsfs/EstSfs_Pos{{chunk}}.csv", chunk = range(config['no_estsfs_chunks']))
        output:
            f"{oDir}/AncestralAllele/Estsfs/EstSfs_Pos_Total.csv"
        wildcard_constraints:
            chunk="\d+"
        threads: 1
        resources: cpus=1, mem_mb=16000, time_min=30
        log: 'logs/Combine_estsfs_pos.log'
        shell:
            "cat {input} > {output}"

    rule combine_estsfs_dicts:
        input:
            expand(f"{oDir}/AncestralAllele/Estsfs/EstSfs_Dict{{chunk}}E.csv", chunk = range(config['no_estsfs_chunks']))
        output:
            f"{oDir}/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
        wildcard_constraints:
            chunk="\d+"
        threads: 1
        resources: cpus=1, mem_mb=16000, time_min=30
        log: 'logs/Combine_estsfs_dicts.log'
        shell:
            "cat {input} > {output}"

    rule run_estsfs:
        input:
            dict=f"{oDir}/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv",
            config=Path(vcfdir, config['estsfs_config']),
            seed=Path(vcfdir, config['estsfs_seed'])
        output:
            text=f"{oDir}/AncestralAllele/Estsfs/outputEtsfs.txt",
            pvalue=f"{oDir}/AncestralAllele/Estsfs/output-pvalues.txt"
        threads: 1 #config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=64000, time_min=2880
        log: 'logs/Run_estsfs.log'
        shell:
            """
            ./scripts/est-sfs {input.config} {input.dict} {input.seed} {output.text} {output.pvalue}
            """

    rule extract_major:
        input:
            f"{oDir}/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
        output:
            temp(f"{oDir}/AncestralAllele/MajorAllele.txt")
        threads: 1 #config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/Extract_major.log'
        shell:
            "./scripts/Extract_major_allele.awk {input} > {output}"

    rule extract_minor:
        input:
            f"{oDir}/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
        output:
            temp(f"{oDir}/AncestralAllele/MinorAllele.txt")
        threads: 1 #config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/Extract_minor.log'
        shell:
            "./scripts/Extract_minor_allele.awk {input} > {output}"

    rule extract_major_outgroup:
        input:
            f"{oDir}/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
        output:
            temp(f"{oDir}/AncestralAllele/MajorOutgroup.txt")
        threads: 1 #config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/Extract_major_outgroup.log'
        shell:
            "./scripts/Extract_major_outgroup.awk {input} > {output}"

    rule extract_estsfs_prob:
        input:
            f"{oDir}/AncestralAllele/Estsfs/output-pvalues.txt"
        output:
            temp(f"{oDir}/AncestralAllele/AncestralProb.txt")
        threads: 1 #config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/Extract_estsfs_prob.log'
        shell:
            "tail -n +9 {input} | cut -f3 -d' ' > {output}"


    rule determine_ancestral:
        input:
            ancProb=f"{oDir}/AncestralAllele/AncestralProb.txt",
            major=f"{oDir}/AncestralAllele/MajorAllele.txt",
            minor=f"{oDir}/AncestralAllele/MinorAllele.txt",
            majorOut=f"{oDir}/AncestralAllele/MajorOutgroup.txt"
        output:
            aa=f"{oDir}/AncestralAllele/AncestralAllele_Vcf_noPos.txt",
            prob=f"{oDir}/AncestralAllele/ProbMajorMinorOutgroup.txt"
        threads: 1 #config['no_estsfs_chunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/DetermineAncestral.log'
        shell:
            """
            paste {input.ancProb} {input.major} {input.minor} {input.majorOut} > {output.prob}
            awk '{{ if(($1 >= 0.5)) {{print $2}} else if (($1 < 0.5)) {{print $3}} else {{print $4}} }}' {output.prob} > {output.aa}
            #rm ProbMajorMinorOutgroup.txt
            """

    rule combine_pos_ancestral:
        input:
            aa=rules.determine_ancestral.output.aa,
            pos=rules.combine_estsfs_pos.output
        output:
            f"{oDir}/AncestralAllele/AncestralAllele_Vcf.txt"
        shell:
            """
	    paste -d "," {input.pos} {input.aa} > {output}
            """
# Comment this out because if this file does not exist, then the procedure in the prepare_files.smk is different
# else:
#     rule move_and_rename_aa:
#         input: config['ancestral_allele']
#         output: "../Project/AncestralAllele/AncestralAllele_Vcf.txt"
#         threads: 1 #config['no_estsfs_chunks']
#         resources: cpus=1, mem_mb=16000, time_min=60
#         log: 'logs/MoveAndRenameAA.log'
#         shell:
#             "mv {input} {output}"
