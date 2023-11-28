#!/usr/bin/env python
# coding: utf-8
import os
import pip
import sys
import json
import subprocess
import numpy as np
import pandas as pd
import pkg_resources

print(sys.executable)

required = {'cyvcf2', 'tsinfer', 'tqdm'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
# try:
#     import tqdm as tqdm
#     print("module 'tqdm' is installed")
# except ModuleNotFoundError:
    print(f"modules {missing} not installed")
    pip.main(['install', *missing])

import cyvcf2
import tsinfer
import tqdm as tqdm

# chromosome = snakemake.wildcards['chromosome']
# vcfFile = snakemake.input[0]
# meta = Path(vcfdir, snakemake.config['meta'])
# sampleFile = snakemake.output[0]
# ploidy = snakemake.config['ploidy']
# chrLength = snakemake.config['chrLength']

args = sys.argv
vcfFile = args[1]
meta = args[2]
sampleFile = args[3]
ploidy = int(args[4])
chrLength = args[5]

#######################################################################
# Define the functions to read in the vcf
#######################################################################
def add_individuals(vcf, samples, populations, metaData, ploidy):
    """
    The function to store the information about which nodes pertain to which sample
    Could be modified in case of a different ploidy
    :param vcf: Vcf file of your samples
    :param samples: tsinfer samples file
    :param populations: tsinfer population object
    :param metaData: meta data file with ID, Pop, Subpop, and time
    """
    for name, population in zip(vcf.samples, populations):
        samples.add_individual(ploidy=ploidy, 
                               metadata={"name": name}, 
                               population=population,
                               #time=list(metaData.Time[metaData.ID == name])[0]
                               )

def add_sites(vcf, samples, ploidy):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
    pos = 0
    progressbar = tqdm.tqdm(total=samples.sequence_length, desc="Read VCF", unit='bp')
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        progressbar.update(variant.POS - pos)
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        if ploidy > 1:
            if any([not phased for _, _, phased in variant.genotypes]):
                raise ValueError("Unphased genotypes for variant at position", pos)

        alleles = [variant.REF]
        #ALT = []
        for string in variant.ALT:
             alleles = alleles + [letter for letter in string]
             #ALT = ALT + [letter for letter in string]

        print(alleles)
        # ignores non-bialllelic sites
        if len(alleles) > 2:
            print('Ignoring non-biallelic site')
            continue

        # if len([*ancestral]) > 1:
        #     print('Ambigous ancestral - assigning as missing')
        #     ancestral_allele = tsinfer.MISSING_DATA
        #     #continue

        ancestral = variant.INFO.get("AA", variant.REF)
        print(f'alleles: {alleles} | ancestal: {[*ancestral]}')
        
        try:
           ancestral_allele = alleles.index(ancestral[0])
        except:
            print('Ancestral not in alleles list, assigning as missing')
            ancestral_allele = tsinfer.MISSING_DATA
            continue

        genotypes = [g for row in variant.genotypes for g in row[0:ploidy]]
        samples.add_site(position=pos,
                             genotypes=genotypes,
                             alleles=alleles,
                             ancestral_allele=alleles.index(ancestral[0]))


def add_populations(vcf, samples, metaData):
    """
    This is a modified add populations function for drones from David Wragg data
    Add tsinfer Population objects for drone subspecies, stored in metaPop object, and return a list od IDs correposning to the VCF samples
    Input: vcf = cyvcf2 VCF object, samples is tsinfer SampleData and metaData is a pandas dataframe with id in ID column and populations in Type column
    pops = dictinary holding name and additional information about populations
    Return: A list of population indexes
    """
    pop_lookup = {}
    metaPop = metaData.iloc[:, 1:].drop_duplicates().dropna(axis = 0, how = 'all')
    sample_subpop = [list(metaData.SubPop[metaData.ID == x]) for x in vcf.samples]
    for subpop, pop in zip(metaPop.SubPop, metaPop.Pop):
        pop_lookup[subpop] = samples.add_population(metadata={"pop": pop, "subpop": subpop})
    return [pop_lookup[subpop[0]] for subpop in sample_subpop]

#######################################################################
# Read in the information and prepare the data
#######################################################################
# Read in the meta file
metaFile = pd.read_csv(meta)

# Create a population (subspecie) list for the samples in the VCF
vcfD = cyvcf2.VCF(vcfFile, strict_gt=True)

# Create samples for haploid data
with tsinfer.SampleData(path=sampleFile,
                        sequence_length=chrLength,
                        num_flush_threads=20, max_file_size=2**30) as samples:
   populations = add_populations(vcf = vcfD, samples = samples, metaData = metaFile)
   print("populations determined")
   add_individuals(vcf = vcfD, samples = samples, populations = populations, metaData = metaFile, ploidy = ploidy)
   print("Inds added")
   add_sites(vcf = vcfD, samples = samples, ploidy = ploidy)

print(
   "Sample file created for {} samples ".format(samples.num_samples)
   + "({} individuals) ".format(samples.num_individuals)
   + "with {} variable sites.".format(samples.num_sites)
)
