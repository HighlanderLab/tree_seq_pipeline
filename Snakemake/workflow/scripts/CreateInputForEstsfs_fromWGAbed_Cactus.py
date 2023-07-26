#!/usr/bin/env python3
import sys
import pandas as pd
import os
from collections import defaultdict
from math import ceil

# This is a script to prepare the input file for est-sfs for the inferrence of the ancestral allele
# The target specie is Apis mellifera and outgroups are Apis cerana, Apis dorsata and Apis florea
# The input files are the output from the multiple whole genome alignment with Mauve and a vcf file of Apis mellifera samples
# A auxiliary input file is the chromosome-name conversion file holding different names for Apis mellifera chromosomes
# 1) to convert .xmfa to .snps file with Mauve and its class SnpExporter
# 2) extract the INFO (AF and AC) for the SNPs in the .snps file from the vcf file
# 3) create a dictionary to hold the est-sfs coded SNPs for Apis mellifera samples from the vcf and outfroups from .snps
args = sys.argv
cycle = int(args[1]) #snakemake.wildcards['chunk']
noCycle = int(args[2])
print("NoCycle is " + str(noCycle))
alignedAlleles = args[3]#snakemake.input[0]
infoFile = args[4]#snakemake.input[1]
outDir = args[5]


######### --- 2 ---###########
# First extract the first three columns from .snps file (the snp pattern, Apis mellifera config and Apis mellifera position)
# Also remove lines with NULL for Apis mellifera from the .snps file
#print("Extracting columns from the WGA bed file and removing the null lines")
#os.system("cut -f1,2,4 " + wgabed + " > AlignedSnps.txt")
# Extract only SNPs aligned in your focal specie (the first character in the chromosome column not ?)
#os.system("""awk '{if ($3 ~ /^\?/) next; else print $0 }' AlignedSnps.txt > AlignedSnps_focal.txt""")
#os.system("""awk -F"\t" '{print $0 FS $1"_"$2}' AlignedSnps_focal.txt > tmp && mv tmp AlignedSnps_focal.txt""")

#Reaplce the chromosome names in outspecies aligned file and sort it
#os.system("qsub Sort_qsub.sh") # This was run on Eddie


# Add FullPos to INFO
#os.system("""awk -F"\t" '{print $1"_"$2}' """ + infoFile + " > VcfFullPos.txt")
#os.system("""awk -F"\t" '{print $0 FS $1"_"$2}' """ + infoFile + " > VCF.INFO")
# # Extract only the VCF SNps from the Alignment file
#os.system("grep -Fwf VcfFullPos.txt AlignedSnps_focal_CarnicaChrSorted.txt > AlignedSnps_focal_CarnicaChrSortedVcf.txt")

# Set the name of the modified .snps file
# Read in the snps file
snps = pd.read_csv(alignedAlleles, sep=" ", header=None)
print("Length of SNPS is " + str(len(snps)))
snps = snps[[0,1,2,4]]
snps.columns =["Chromosome", "Position", "SnpPattern", "FullPos"]
nRowCycle = ceil(len(snps) / noCycle)
print("Nrow cycle is " + str(nRowCycle))
snps = snps[nRowCycle * cycle: nRowCycle * (cycle+1)]
print(len(snps))
# Focal chromosomes
# Combine the columns into a sigle column position
# Positions in WGAbed are 0 based (hence +1)!


# Prepare a list of SNPs to extract from the vcf file
#print("Extracting SNP info from the vcf file")
# Use vcftools to extract the alternate allele count (AC) and allele frequency (AF) for the given set of SNPs from the vcf file
#os.system("vcftools --vcf " + vcfFile + " --positions " + focalPos + "  --get-INFO AC --get-INFO AF --out Test_Focal")
# Read in the vcf info
snpsInfo = pd.read_csv(infoFile, sep=" ")
snpsInfo.loc[:, "Pos"] = snpsInfo['CHROM'].apply(str).str.cat(snpsInfo['POS'].apply(str),sep="_")



######### --- 3 ---###########
# Create representation of the four SNPs for the est-sfs
snpDict = {"A": (1,0,0,0), "C": (0,1,0,0), "G":(0,0,1,0), "T":(0,0,0,1), "?":(0,0,0,0), "N":(0,0,0,0)}

# Create a list of species following the order on the .snps file
print("Creating a dictionary with est-sfs coded snps")
species = {"Amel":0, "Acer":1, "Ador":2, "Aflo":3}

# Create a dictionary holding est-sfs coded info for each snp
specieValues = defaultdict()


# # Create a loop to write the snp info to dictionary
print("Creating an est-sfs dictionary of alelle counts for each specie.")
for snpPos, snpPattern in iter(zip(snps.FullPos, snps.SnpPattern)):
    if len(snpPattern.split(",")[0]) == 1:
        # Check whether the SNP is also in the vcf file
        if snpPos in list(snpsInfo.Pos):
            snpSpecDict = {}
            snpLine = snpsInfo.query("Pos == @snpPos")
            # Extract the alleles of the focal species and the outgroups
            specieAlleles = snpPattern.split(",") #[snpPattern.split(",")[x] for x in species.values()]
            focalAllele = snpPattern.split(",")[0].upper()
            # Extract count and frequency for the alternate allele and compute the total allele count
            # This is now an iterator!!! Iterate only once!!!! Even with a sum!
            altFreq = {altAl:float(altValue) for altAl, altValue in zip(list(snpLine.ALT)[0].split(","), list(snpLine.AF.astype(str))[0].split(","))}
            altCount = {altAl:int(altValue) for altAl, altValue in zip(list(snpLine.ALT)[0].split(","), list(snpLine.AC.astype(str))[0].split(","))}
            totalCount = max([round(ac / af) if af != 0 else 0 for ac, af in zip(altCount.values(), altFreq.values())])
            altSum = sum(altCount.values())
            refSum = totalCount - altSum
            refAllele = list(snpLine.REF)[0]

            # Create an est-sfs dict for the reference allele
            vcfDict = []
            refCountDict = tuple((tuple(x*refSum for x in snpDict.get(refAllel)) for refAllel in  refAllele if refAllele in snpDict.keys()))
            if refCountDict:
                vcfDict.append(refCountDict[0])
            # Create an est-sfs dict for the alternate alleles
            altCountDict = list((tuple(x*altAlleleCount for x in snpDict.get(altAllele))
                                 for altAllele, altAlleleCount in zip(altCount.keys(), altCount.values())
                                 if altAllele in snpDict.keys()))
            vcfDict = vcfDict + altCountDict

            # Zip and sum all the tuples within the vcfDict list
            FocalCount = tuple([sum(x) for x in zip(*vcfDict)])

            # Add est-sfs code for the SNP into the dictionary
            # Create a dict for the three outgroups from the mauve output
            snpCountDict = [snpDict[x.upper()] for x in specieAlleles[1:]]
            # Create a dict for the Apis mellifera from the mauve output
            FocalSnpCountDict = snpDict[focalAllele.upper()]
            # Sum the mauve Apis mellifera and vcf Apis mellifera alleles
            FocalCount = tuple([sum(x) for x in zip(FocalCount, FocalSnpCountDict)])
            # Add the focal summed counts as the first element of the list
            snpCountDict.insert(0, FocalCount)
            # Assign the list of tuples for each Apis species to the dictionary under the key = SNPName
            specieValues[snpPos] = snpCountDict
            snpSpecDict[snpPos] = snpCountDict

pd.DataFrame.from_dict(specieValues, orient='index').to_csv(outDir + "/EstSfs_Dict" + str(cycle) + ".csv", header=None, sep="\t")

# Write the dictionary of tupples into a dataframe
#pd.DataFrame.from_dict(specieValues, orient='index').to_csv("EstSfs_Dict.csv", index=None, header=None, sep="\t")
# Remove the spaces in the output file
#os.system("""sed -i "s/ //g" EstSfs_Dict.csv""")
# Set the correct separators for the file
#os.system("""awk -F "\t" '{print $1"\t"$2" "$3" "$4}' EstSfs_Dict.csv > tmp && mv tmp EstSfs_Dict.csv""")
# Remove the parenthesis from the file
#os.system("""sed -i "s/(//g" EstSfs_Dict.csv""")
#os.system("""sed -i "s/)//g" EstSfs_Dict.csv""")
