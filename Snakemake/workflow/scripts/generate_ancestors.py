import sys
import tsinfer
import numcodecs
import numpy as np

args = sys.argv

samples = tsinfer.load(args[1])
sites_to_exclude = np.loadtxt(args[2])
output = args[3]
threads = int(args[4])

gzip_compressor = numcodecs.GZip(level=1)

ancestors = tsinfer.generate_ancestors(
    sample_data=samples,
    path=output,
    exclude_positions=sites_to_exclude,
    progress_monitor=True,
    num_threads=threads,
    compressor=gzip_compressor,
)