#!/usr/bin/env python
# coding: utf-8

import os
import sys
import json
import cyvcf2
import tsinfer
import numcodecs
import numpy as np
import pandas as pd
from tqdm import tqdm

# create a GZip compressor instance with a specific compression level
gzip_compressor = numcodecs.GZip(level=1)

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
    
    # to account for same subpop in different countries - creating a 'fake subpop' to be used as key
    # if subpop is repeated new code is subpop+origin, else subpop
    d = [a if not (sum(j == a for j in metaPop.SubPop[:i])) else f'{a}{c}' for (i,a),c in zip(enumerate(metaPop.SubPop),metaPop.Origin)]
    metaPop['code'] = d
    # add code col to the full meatafile
    metaData = pd.merge(metaData, metaPop[['Pop', 'SubPop', 'Origin', 'code']], on = ['Pop', 'SubPop', 'Origin'], how = 'inner')
    
    # list population per sample based on recoded subpop
    sample_subpop = [list(metaData.code[metaData.ID == x]) for x in vcf.samples]
    
    # the dictionary key will be based on the recoding - so keeps the correct amount of population
    # but the metadata takes in the real subpop value
    for subpop, pop, origin, lat, lon, code in zip(metaPop.SubPop, metaPop.Pop, metaPop.Origin, metaPop.Lat, metaPop.Lon, metaPop.code):
        pop_lookup[code] = samples.add_population(metadata={"pop": pop, "subpop": subpop, "origin": origin, "coord1": lat, "coord2":lon})
    return [pop_lookup[code[0]] for code in sample_subpop]

def add_diploid_individuals(vcf, samples, populations):
    for iid, population in zip(vcf.samples, populations):
        samples.add_individual(ploidy=2, metadata={'ID': iid}, population=population)


def add_diploid_sites(vcf, samples):
    """
    Read the sites in the vcf and add them to the samples object.
    """
    # You may want to change the following line, e.g. here we allow
    # "*" (a spanning deletion) to be a valid allele state
    allele_chars = set("ATGCatgc*")
    
    pos = 0
    
    progressbar = tqdm(total=samples.sequence_length, desc="Read VCF", unit='bp')
    
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        progressbar.update(variant.POS - pos)
        
        if pos == variant.POS:
            print(f"Duplicate entries at position {pos}, ignoring all but the first")
            continue
        else:
            pos = variant.POS
        
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        
        alleles = [variant.REF.upper()] + [v.upper() for v in variant.ALT]
        
        ancestral = variant.INFO.get("AA", ".")  # "." means unknown
        
        if(ancestral == 'ambiguous' or ancestral not in alleles): # set ancestral == reference
            ancestral_allele = 0
        elif ancestral == '-1':
            ancestral_allele = -1    
        else:
            ancestral_allele = alleles.index(ancestral)
        
        # Check we have ATCG alleles
        for a in alleles:
            if len(set(a) - allele_chars) > 0:
                print(f"Ignoring site at pos {pos}: allele {a} not in {allele_chars}")
                continue
        
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [g for row in variant.genotypes for g in row[0:2]]
        samples.add_site(pos, genotypes, alleles, ancestral_allele=ancestral_allele)

def process_samples(file, metadata, sample_path, chr_length, compressor, threads):
    with tsinfer.SampleData(path=sample_path, 
                            sequence_length=chr_length,
                            compressor=compressor, 
                            num_flush_threads=threads, max_file_size=2**30) as samples:
        populations = add_populations(vcf=file, samples=samples, metaData=metadata)
        print("populations determined")
        add_diploid_individuals(vcf=file, samples=samples, populations=populations)
        print("individuals added")
        add_diploid_sites(vcf=file, samples=samples)
        print("sites added")
        
        print(
            "Sample file created for {} samples ".format(samples.num_samples)
            + "({} individuals) ".format(samples.num_individuals)
            + "with {} variable sites.".format(samples.num_sites),
            flush=True,
        )

# -----------------------------------------------------------------------------------------

args = sys.argv
vcfFile = args[1]
sampleFile = args[2]
meta = args[3]
chrLength = args[4]
nthreads = int(args[5])

metaFile = pd.read_csv(meta)

# Create a population (subspecie) list for the samples in the VCF
vcfD = cyvcf2.VCF(vcfFile, strict_gt=True)

# Create samples for haploid data
try:
    process_samples(vcfD, metaFile, sampleFile, chrLength, gzip_compressor, nthreads)
except Exception as e:
    print(f"Error encountered: {e}. Retrying with num_flush_threads set to 16.")
    process_samples(vcfD, metaFile, sampleFile, chrLength, gzip_compressor, 16)
        
# with tsinfer.SampleData(path=sampleFile,
#                         sequence_length=chrLength,
#                         num_flush_threads=nthreads, max_file_size=2**30) as samples:
#    populations = add_populations(vcf = vcfD, samples = samples, metaData = metaFile)
#    print("populations determined")
#    add_diploid_individuals(vcf = vcfD, samples = samples, populations = populations)
#    print("individuals added")
#    add_diploid_sites(vcf = vcfD, samples = samples)
#    print("sites added")

# -----------------------------------------------------------------------------------------
# def execute_with_threads(nthreads):
#     try:
#         with tsinfer.SampleData(path=sampleFile, 
#                                 sequence_length=chrLength,
#                                 compressor=gzip_compressor, 
#                                 num_flush_threads=nthreads, max_file_size=2**30) as samples:
#             populations = add_populations(vcf=vcfD, samples=samples, metaData=metaFile)
#             print("populations determined")
#             add_diploid_individuals(vcf=vcfD, samples=samples, populations=populations)
#             print("individuals added")
#             add_diploid_sites(vcf=vcfD, samples=samples)
#             print("sites added")
#     except Exception as e:
#         print(f"Error encountered: {e}. Retrying with num_flush_threads set to 16.")
#         with tsinfer.SampleData(path=sampleFile,
#                                 sequence_length=chrLength,
#                                 compressor=gzip_compressor,
#                                 num_flush_threads=16, max_file_size=2**30) as samples:
#             populations = add_populations(vcf=vcfD, samples=samples, metaData=metaFile)
#             print("populations determined")
#             add_diploid_individuals(vcf=vcfD, samples=samples, populations=populations)
#             print("individuals added")
#             add_diploid_sites(vcf=vcfD, samples=samples)
#             print("sites added")
# execute_with_threads(nthreads)
