import sys
import tskit
import subprocess
# upgrade to the latest version of tsdate
subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "tsdate"])
import tsdate

args = sys.argv
ts = tskit.load(args[1])
generation_interval = float(args[2])
mu = float(args[3])
output = args[4]

# Ne = ts.diversity()/4/mu  

processed_ts = tsdate.preprocess_ts(
    ts,
    split_disjoint=True,
    keep_unary=False,
    )

dated_ts = tsdate.variational_gamma(tree_sequence = processed_ts, 
                                  mutation_rate = mu, 
                                  eps = generation_interval,
                                  max_iterations = 100,
                                  time_units = 'generations',
                                  progress = True,)

dated_ts.dump(output)