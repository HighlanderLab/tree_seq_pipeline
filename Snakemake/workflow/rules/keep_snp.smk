# rule gatk_index:
#     input:
#         vcf = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz',
#         idx = f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz.csi' 
#     output: temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_phased.vcf.gz.tbi'),
#     conda: 'gatk-env'
#     benchmark: 
#         'benchmark/{chromosome}.gatk_index.benchmark.txt'
#     shell:
#         r"""
#         gatk IndexFeatureFile -I {input.vcf}
#         """

# rule filter_snps:
#     input: 
#         vcf = rules.phase.output,
#         idx = rules.index_phased.output,
#     output: 
#         vcf_file = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_snps.vcf.gz'),
#         idx_file = temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_snps.vcf.gz.csi'),
#     params: config['reference_genome'],
#     conda: 'gatk-env'
#     benchmark: 'benchmarks/{chromosome}.filter_snps.benchmark.txt'
#     shell:
#         r"""
#         gatk SelectVariants \
#             -R {params} \
#             -V {input.vcf} \
#             --select-type-to-include SNP \
#             -O {output.vcf_file} 

#         bcftools index {output.vcf_file}           
#         """

# rule index:
#     input: f'{vcfdir}/{{chromosome}}/{{chromosome}}_snps.vcf.gz',
#     output: temp(f'{vcfdir}/{{chromosome}}/{{chromosome}}_snps.vcf.gz.csi')
#     conda: 'bcftools'
#     benchmark: 'benchmarks/{chromosome}.index_filtered_snps.benchmark.txt'
#     shell:
#         r"""
#         bcftools index {input}
#         """
