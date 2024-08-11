# Script for (iteractively) dating the tree sequence 
# trees must be dated, re-inferred and dated again
# to get the correct dates (test number of iterations)

# if this is the first iteration, start from the existing tree sequence;
# otherwise, re-infer the trees using the dated sampleFile and date the trees again

input_ts = 

print(f'Starting dating of tree sequence {input_ts}...')

rule date_trees:
    input: rules.simplify_trees.output
    output: f'{_input:n}.dated'
    conda: ''


print('A rough estimate of the effective population size is', ts.diversity() / (4 * 1e-6))

