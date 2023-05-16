1) Enter an interactive mode on Eddie (you don't need that much memory)
2) Navigate to HighlanderLab/share/Snakemake. Here, you'll find
	a) Data: in here, there is the meta data, the ancestral allele and RawVCF with two VCF files with different samples, but both with all chromosomes
	b) tree_seq_pipeline: this it the github directory but with an added RunSnakemake.sh file
3) Move to "tree_seq_pipeline/Snakemake/workflow". Here, you'll find the RunSnakemake.sh. This file:
	a) loads the miniconda module
	b) activates the correct environment
	c) loads bcftools and vcftools
	d) runs the Snakefile
4) The Snakefile creates:
	a) tree_seq_pipeline/Snakemake/Project directory that stores the output
	b) in the Data directory (specified in the config), it creates a folder for each chromosome (e.g. Chr1, Chr2) and the final VCFs
	c) You need to remove these files/folders if you want to run it again, otherwise it will not run!
	

NOTE!!!
1) There might be problems with using the correct python (I had that). Before running the Snakemake, just activate the conda env and type "which python". If the path is not "/exports/cmvm/eddie/eb/groups/HighlanderLab/anaconda/envs/HLab_tsinfer/bin/python", then run 
"export PATH=/exports/cmvm/eddie/eb/groups/HighlanderLab/anaconda/envs/HLab_tsinfer/bin/:$PATH" and check again "which python". It has to point to the env python


