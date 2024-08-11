# if there is only one vcf files, use the output from split_and_move as input
# for the next step samples
# else, compare the files and filter out overlapping then merge the vcfs into
# one and index and use the output from the merge as input for the next step
if len(allFiles) == 1:
    input_vcf = rules.split_and_move_vcfs.output
    input_index = f'{rules.split_and_move_vcfs.output}.csi'
    input_summary = 'None'

else:
    input_to_filter = [f'{vcfdir}/{{chromosome}}/{{chromosome}}_{file}.vcf.gz' for file in allFiles]
    
    filtered_output = [f'{vcfdir}/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz' for file in allFiles]
    filtered_output_idx = [f'{vcfdir}/{{chromosome}}/{{chromosome}}_{file}.filtered.vcf.gz.csi' 
                            for file in allFiles]                    
    merged_output   = f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.vcf.gz'

    rule compare_and_filter_vcfs:
        input: input_to_filter 
        output: temp(filtered_output)
        params:
            file_names = [f'{{chromosome}}_{file}' for file in allFiles],
            path = f'{vcfdir}/{{chromosome}}'
        benchmark: 'benchmarks/{chromosome}.sample.benchmark.txt'    
        run:
            import os
            import shlex
            import itertools
            import subprocess
            from glob import glob

            path_names = [f'{params.path}/{name}' for name in params.file_names]
            pairs = list(itertools.combinations(path_names, 2))
            new_ext = '.filtered.vcf.gz'
            to_rename = path_names.copy()

            # to avoid shell injection, sanitize the file names
            # -- ensure speciall shell characters are escaped
            sanitized_path = shlex.quote(params.path)
            sanitized_new_ext = shlex.quote(new_ext)
            sanitized_idx = shlex.quote(new_ext + '.csi')
            
            for i, pair in enumerate(pairs):
                print(f'pair {i} {pair}')
                
                # Run vcftools to find if there are overlapping samples between pairs
                command = f"""
                        module load anaconda
                        source activate bcftools
                        vcftools --gzvcf {os.path.join(pair[0] + '.vcf.gz')} \
                                 --gzdiff {os.path.join(pair[1] + '.vcf.gz')} \
                                 --diff-indv
                        """
                subprocess.run(command, shell=True, check=True, cwd=sanitized_path)   
                
                # if vcftools find any overlapping samples, remove them from the second vcf file
                if "B" in subprocess.run("awk '{if ($2 == \"B\") print $2}' out.diff.indv_in_files", 
                                            shell=True, 
                                            capture_output=True, 
                                            text=True, 
                                            cwd=sanitized_path).stdout:
                    print(f'Overlapping samples found!! Removing samples from {pair[1]}')
                    command = f"""
                            grep B out.diff.indv_in_files | cut -f1 > samples.remove
                            rm out.diff.indv_in_files
                            """
                    subprocess.run(command, shell=True, check=True, cwd=sanitized_path)     
                    command = f"""
                              bcftools view -S ^samples.remove {os.path.join(pair[1] + '.vcf.gz')} \
                                            -Oz -o {os.path.join(pair[1] + sanitized_new_ext)}
                              
                              bcftools index -f {os.path.join(pair[1] + sanitized_new_ext)}      

                              rm samples.remove
                            """
                    subprocess.run(command, shell=True, check=True, cwd=sanitized_path)
                    # remove the file that has been sampled from renaming list
                    to_rename.remove(pair[1])
                    print(f'Files to rename are: {to_rename}')

                # if there are no overlapping samples, 
                # just rename the files that still have the original extension
                else:
                    print('No overlapping samples')

            print(f'Finished comparing and sampling vcfs')

            # rename the files that have not been sampled
            for name in to_rename:
                sanitized_name = shlex.quote(name)
                print(f'Renaming {sanitized_name}')
                command = f"""
                        mv {os.path.join(sanitized_name + '.vcf.gz')} {os.path.join(sanitized_name + sanitized_new_ext)}
                        mv {os.path.join(sanitized_name + '.vcf.gz.csi')} {os.path.join(sanitized_name + sanitized_idx)}
                        """
                subprocess.run(command, shell=True, check=True, cwd=sanitized_path)

            # Remove all files that do not contain the string new_ext
            command = f"""
                        find {sanitized_path} -type f ! -name '*{sanitized_new_ext}*' -exec rm -f {{}} \;
                        """
            subprocess.run(command, shell=True, check=True, cwd=sanitized_path)

    rule merge:
        input: rules.compare_and_filter_vcfs.output
        output: merged_output 
        params:
            files_path = f'{vcfdir}/{{chromosome}}/',
            indexes = filtered_output_idx
        conda: "bcftools"
        threads: 32
        resources: cpus=32, mem_mb=2048000, time_min=1200
        benchmark: 'benchmarks/{chromosome}.merge.benchmark.txt'
        shell:
            r"""
            printf '%s\n' {params.files_path}*filtered*.gz > {wildcards.chromosome}.merge
            
            bcftools merge \
                -m snps \
                -l {wildcards.chromosome}.merge \
                --threads {threads} -O z -o {output}
	    
            rm {wildcards.chromosome}.merge {params.indexes}
            """

    rule index_merged:
        input: rules.merge.output,
        output: f'{vcfdir}/{{chromosome}}/{{chromosome}}_merged.vcf.gz.csi'
        conda: 'bcftools'
        benchmark: 'benchmarks/{chromosome}.indxmerged.benchmark.txt'
        shell:
            r"""
            bcftools index -f {input}
            """                    
    
    rule vcf_summary:
        input: 
            vcf_file = rules.merge.output,
            index = rules.index_merged.output,
        output:
            summary = f'{vcfdir}/{{chromosome}}/{{chromosome}}_summary.vchk',
            # done = touch(f'{vcfdir}/{{chromosome}}/merged_summary.done')
        params:
            plots = f'{vcfdir}/{{chromosome}}/{{chromosome}}_summary/'
        conda: 'bcftools'
        benchmark: 'benchmarks/{chromosome}.summary.benchmark.txt'
        shell:
            r"""
            bcftools stats {input.vcf_file} > {output.summary}
            plot-vcfstats -p {params.plots} {output.summary}
            """

    input_vcf = rules.merge.output
    input_index = rules.index_merged.output
    input_summary = rules.vcf_summary.output.summary

# else, use the output from split_and_move as input for the next step
# else:
#     input_vcf = rules.split_and_move_vcfs.output