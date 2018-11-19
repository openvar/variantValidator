import sys
import vvTestFunctions as fn

argv=sys.argv
#print(argv)
if len(argv)!=4:
    print("Syntax: python mergeInputVariants.py path1 path2 pathOut")
    print("Pass two paths (path1 and path2) to this script to merge the variant strings within")
    print("and save the output to the file pathOut, removing duplicates in the process.")
else:
    v1=fn.loadVariantFile(argv[1])
    v2=fn.loadVariantFile(argv[2])
    vOut=fn.mergeVariantList(v1,v2)
    fn.saveVariantFile(argv[3],vOut)
    print("Merged "+str(len(v1))+" + "+str(len(v2))+" into "+str(len(vOut))+" unique variants."   )
