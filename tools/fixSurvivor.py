import sys
from itertools import groupby
from pysam import VariantFile


vcfFile =VariantFile(sys.argv[1])
out=VariantFile(sys.argv[2],'w',header=vcfFile.header)

for rec in vcfFile.fetch():
    if "SVTYPE" in rec.info and rec.info["SVTYPE"] =="DEL":
        print(rec.alts)
        rec.alts=rec.ref[0]
    out.write(rec)
    
