import sys

inputVCF=open(sys.argv[1])
outputVCF=open(sys.argv[2],'w')
toolsVCF=sys.argv[3:]

ids2GT={}
for t in toolsVCF:
    for l in open(t):
        if l[0] != "#":
            l=l.split("\t")
            ids2GT[l[2]]=(l[8],l[9])
            
for l in inputVCF:
    if l[0] == "#":
        outputVCF.write(l)
    else:
        l=l.strip().split("\t")
        l[8]=ids2GT[l[2]][0]
        l[9]=ids2GT[l[2]][1]
        outputVCF.write("\t".join(l))
