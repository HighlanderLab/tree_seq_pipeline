import sys
import tsinfer

args = sys.argv
ancestors = tsinfer.load(args[1])
uprtime = float(args[2])
lwertime = float(args[3])
output = args[4]

ancestors.truncate_ancestors(
    path=output, 
    max_file_size=2**38,
    lower_time_bound=lwertime,
    upper_time_bound=uprtime,
    length_multiplier=2,
)

