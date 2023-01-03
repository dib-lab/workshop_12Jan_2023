import sys
from itertools import groupby
from pysam import VariantFile

def fasta_iter(fasta_name):
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)



vcfFile =VariantFile(sys.argv[1])
refIter=fasta_iter(sys.argv[2])
out=VariantFile(sys.argv[3],'w',header=vcfFile.header)
reference={} 

for header,seq in refIter:
    reference[header.split(" ")[0]]=seq


numCorrect=0
numFailed=0
numNotInReference=0
numNegative=0
total=0
for rec in vcfFile.fetch():
    if "SVTYPE" not in rec.info:
        allele=rec.ref
        check=reference[rec.chrom][rec.start:rec.start+len(allele)]
        if allele != check:
            print("Mismatch %s != %s"%(allele,check))
            print(rec)
            numFailed+=1
        else:
            out.write(rec)
            numCorrect+=1
        continue

    if rec.info["SVTYPE"] in ["TRA","INV"] or ("SVLEN" not in rec.info):
        continue
    total+=1
    svlen=rec.info["SVLEN"]
    if isinstance(svlen,tuple):
        svlen=svlen[0]

    if rec.pos + svlen <0 or rec.stop <0 or rec.start + svlen <0 or rec.start -svlen < 0 or rec.pos -svlen < 0:
        print("Negative end position")
        print(rec)
        numNegative+=1
        continue
    if rec.chrom not in reference:
        numNotInReference+=1
        print("Not in Reference")
        print(rec)
        continue
    allele=rec.ref
    check=reference[rec.chrom][rec.start:rec.start+len(allele)]
    if allele != check:
        print("Mismatch %s != %s"%(allele,check))
        print(reference[rec.chrom][rec.start-5:rec.start+len(allele)+5])
        print(rec)
        numFailed+=1
    else:
        out.write(rec)
        numCorrect+=1


print("Correct = %d"%numCorrect)
print("Failed = %d" % numFailed)
print("Negative END= %d" % numNegative)
print("Not in ref = %d" %numNotInReference)
print("Total = %d" % total)


# Popsize=float(len(sys.argv)-2)*2

# variants={}
# breedName=sys.argv[1]

# for vcfs in sys.argv[2:]:
#     vcf_in = VariantFile(vcfs)
#     for rec in vcf_in.fetch():
#         if rec.id not in variants:
#             variants[rec.id]=0
#         sampleName = rec.samples.keys()[0]     
#         if rec.samples[sampleName]['GT'] == (0,1) :
#             variants[rec.id]+=1
#         if rec.samples[sampleName]['GT'] == (1,1) :
#             variants[rec.id]+=2


# print("variantId,%s"%(breedName))            
# for vid,found in variants.items():
#     popPerc=(float(found) /popSize)*100.0 
#     print("%s,%.2f"%(vid,popPerc,))
    
