from itertools import combinations
import os

args = sys.argv
samples = sys.argv[1]
vcflist_input = args[2]
vcflist_output = args[3]
samplelist=[]
track=[]

#            os.system("module load igmm/apps/bcftools/1.9")
for ifile in samples:
    with open(ifile, 'r') as f:
        samplelist.append(f.read().splitlines())

# Take sublist in pairs (return index [0] and sublist [1]) and compare
for a, b in combinations(enumerate(samplelist), 2):
        # if there are overlapping entries, remove duplicates
        if not set(a[1]).isdisjoint(b[1]) == True:
            [b[1].remove(element) for element in a[1] if element in b[1]]
            # keeps track of the sublists that changed using index
            track.append(b[0])

# filter files if sample list has changed, otherwise only rename
for i in range(len(samplelist)):
    vcf=vcflist_input[i]

    ovcf=vcflist_output[i]
    samples=samplelist[i]
    if (i in set(track)) and (len(samples) != 0):
        shell('bcftools view -S {samples} --force-samples {vcf} -O z -o {ovcf}')
        shell('bcftools index {ovcf}')

    elif len(samples) == 0:
        shell('touch {ovcf}.ignore')
    else:
        shell('bcftools index {ovcf}')
