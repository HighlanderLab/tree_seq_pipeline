#!/usr/bin/env python
# coding: utf-8

import cyvcf2
import tsinfer
import pandas as pd
import json
import numpy as np
import sys
import os

sampleFile = snakemake.input[0]
outputFile = snakemake.output[0]

#######################################################################
# Do the inference and write the outputs
#######################################################################

# Do the inference on the 10 SNPs
samples = tsinfer.load(sampleFile)

ts = tsinfer.infer(samples)
print(
    "Inferred tree sequence `{}`: {} trees over {} Mb".format(
        "drone_ts", ts.num_trees, ts.sequence_length / 1e6
    )
)

# Check the metadata
for sample_node_id in ts.samples():
    individual_id = ts.node(sample_node_id).individual
    population_id = ts.node(sample_node_id).population
    print(
        "Node",
        sample_node_id,
        "labels genome sampled from",
        json.loads(ts.individual(individual_id).metadata),
        "in",
        json.loads(ts.population(population_id).metadata)["subspecies"],
    )

ts.dump(outputFile)
