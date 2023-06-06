from itertools import combinations
import os
import sys

args = sys.argv
argsList = args[1:]
print(len(argsList))
part = int(len(argsList) / 3)
inputSamples = argsList[0:part]
vcflist_input = argsList[part:(part*2)]
vcflist_output = argsList[(part*2):(part*3)]
print(vcflist_input)
samplelist=[]
track=[]

for ifile in inputSamples:
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
        os.system('bcftools view -S {samples} --force-samples ' + vcf + ' -O z -o ' + ovcf)
        os.system('bcftools index ' + ovcf)

    elif len(samples) == 0:
        os.system('touch ' + ovcf + '.ignore')
    else:
        os.system('bcftools view ' + vcf + ' -O z -o ' + ovcf)
        os.system('bcftools index ' + ovcf)
