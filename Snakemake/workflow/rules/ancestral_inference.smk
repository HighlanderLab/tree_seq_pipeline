from __future__ import division
import sys

if config['ancestralAllele'] == None:
    rule extract_aligned_alleles_from_vcf:
        input:
            alignedFocal=config['alignedFocal'],
            vcf=config['rawVcf']
        output:
            info="{workdir}/AncestralAllele/RawVcfInfo.INFO",
            log=temp("{workdir}/AncestralAllele/RawVcfInfo.log")
        params:
            prefix="{workdir}/AncestralAllele/RawVcfInfo"
        conda: "bcftools"
        threads: 1
        resources: cpus=1, mem_mb=64000, time_min=300
        log: 'logs/ExtractAlignedAllelsFromVcf.log'
        shell:
            "scripts/ExtractPosVcf.sh {input.vcf} {input.alignedFocal} {params.prefix}"


    rule extract_snps_from_info:
        input:
            "../Project/AncestralAllele/RawVcfInfo.INFO"
        output:
            "../Project/AncestralAllele/RawVcfInfo_SNPs.INFO"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=60
        log: 'logs/ExtractSnpsFromInfo.log'
        shell:
            """
            head -n1 {input} > Header
            sed -i "s/\t/ /g" Header
            ./scripts/ExtractSNPsFromInfo.awk {input} > tmp
            cat Header tmp > {output}
            rm Header tmp
            """

    rule extract_pos_aligned_alleles_from_vcf:
        input:
            rules.extract_snps_from_info.output
        output:
            "../Project/AncestralAllele/RawVcfFullPos.txt"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=60
        log: 'logs/ExtractPosAlignedAllelesFromVcf.log'
        shell:
            """
            awk '{{print $1" "$2}}' {input} > {output}
            """

    rule extract_vcf_alleles_from_aligned:
        input:
            alignedFocal=config['alignedFocal'],
            vcfPos=rules.extract_pos_aligned_alleles_from_vcf.output
        output:
            "../Project/AncestralAllele/AlignedSnps_focal_SortedRawVcf.txt"
        threads: 1
        resources: cpus=1, mem_mb=32000, time_min=120
        log: 'logs/ExtractVcfAllelesFromAligned.log'
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
            expand("../Project/AncestralAllele/Estsfs/EstSfs_Dict{{chunk}}.csv")
        wildcard_constraints:
            chunk="\d+"
        params:
            outDir="../Project/AncestralAllele/Estsfs",
            noCycle=config['noEstSfsChunks']
        threads: config['noEstSfsChunks']
        resources: cpus=1, mem_mb=64000, time_min=360
        log: 'logs/CreateEstsfsDicts{chunk}.log'
        shell:
            "python scripts/CreateInputForEstsfs_fromWGAbed_Cactus.py {wildcards.chunk} {params.noCycle} {input.alignedAlleles} {input.vcfAlleles} {params.outDir}"
            #CreateInputForEstsfs_Loop.sh This is qsub

    rule edit_estsfs_dicts:
        input:
            rules.create_estsfs_dicts.output
        output:
            tmp1=temp("../Project/AncestralAllele/Estsfs/tmp1{chunk}.csv"),
            tmp2=temp("../Project/AncestralAllele/Estsfs/tmp2{chunk}.csv"),
            editedDict="../Project/AncestralAllele/Estsfs/EstSfs_Dict{chunk}E.csv"
        params:
            chunk=config['noEstSfsChunks']
        wildcard_constraints:
            chunk="\d+"
        threads: config['noEstSfsChunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/EditEstsfsDicts{chunk}.log'
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

    rule combine_estsfs_dicts:
        input:
            expand("../Project/AncestralAllele/Estsfs/EstSfs_Dict{chunk}E.csv", chunk = range(config['noEstSfsChunks']))
        output:
            "../Project/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
        wildcard_constraints:
            chunk="\d+"
        threads: 1
        resources: cpus=1, mem_mb=16000, time_min=30
        log: 'logs/CombineEstsfsDicts.log'
        shell:
            "cat {input} > {output}"

    rule run_estsfs:
        input:
            dict="../Project/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv",
            config=config['estsfsConfig'],
            seed=config['estsfsSeed']
        output:
            text="../Project/AncestralAllele/Estsfs/outputEtsfs.txt",
            pvalue="../Project/AncestralAllele/Estsfs/output-pvalues.txt"
        threads: 1 #config['noEstSfsChunks']
        resources: cpus=1, mem_mb=64000, time_min=2880
        log: 'logs/RunEstsfs.log'
        shell:
            """
            ./scripts/est-sfs {input.config} {input.dict} {input.seed} {output.text} {output.pvalue}
            """

    rule extract_major:
        input:
            "../Project/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
        output:
            temp("../Project/AncestralAllele/MajorAllele.txt")
        threads: 1 #config['noEstSfsChunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/ExtractMajor.log'
        shell:
            "./scripts/ExtractMajorAllele.awk {input} > {output}"

    rule extract_minor:
        input:
            "../Project/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
        output:
            temp("../Project/AncestralAllele/MinorAllele.txt")
        threads: 1 #config['noEstSfsChunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/ExtractMinor.log'
        shell:
            "./scripts/ExtractMinorAllele.awk {input} > {output}"

    rule extract_major_outgroup:
        input:
            "../Project/AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
        output:
            temp("../Project/AncestralAllele/MajorOutgroup.txt")
        threads: 1 #config['noEstSfsChunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/ExtractMajorOutgroup.log'
        shell:
            "./scripts/ExtractMajorOutgroup.awk {input} > {output}"

    rule extract_estsfs_prob:
        input:
            "../Project/AncestralAllele/Estsfs/output-pvalues.txt"
        output:
            temp("../Project/AncestralAllele/AncestralProb.txt")
        threads: 1 #config['noEstSfsChunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/ExtractEstsfsProb.log'
        shell:
            "tail -n +9 {input} | cut -f3 -d' ' > {output}"

    rule determine_ancestral:
        input:
            ancProb="../Project/AncestralAllele/AncestralProb.txt",
            major="../Project/AncestralAllele/MajorAllele.txt",
            minor="../Project/AncestralAllele/MinorAllele.txt",
            majorOut="../Project/AncestralAllele/MajorOutgroup.txt"
        output:
            aa="../Project/AncestralAllele/AncestralAllele_Vcf.txt",
	    prob="../Project/AncestralAllele/ProbMajorMinorOutgroup.txt"
        threads: 1 #config['noEstSfsChunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/DetermineAncestral.log'
        shell:
            """
            paste {input.ancProb} {input.major} {input.minor} {input.majorOut} > {output.prob}
            awk '{{ if(($1 >= 0.5)) {{print $2}} else if (($1 < 0.5)) {{print $3}} else {{print $4}} }}' {output.prob} > {output.aa}
            #rm ProbMajorMinorOutgroup.txt
            """
else:
    rule move_and_rename_aa:
        input: config['ancestralAllele']
        output: "../Project/AncestralAllele/AncestralAllele_Vcf.txt"
        threads: 1 #config['noEstSfsChunks']
        resources: cpus=1, mem_mb=16000, time_min=60
        log: 'logs/MoveAndRenameAA.log'
        shell:
            "mv {input} {output}"
