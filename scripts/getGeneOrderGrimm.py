import sys

# this script reads in all the gene order files
# limits the analysis to those clusters that are present in
# all the supplied genomes
# and renumbers the result

# first, it is necessary to get the list of what is present

noFiles = len(sys.argv)

presentGenes = {}
geneLists = []
genomes = {}

# read the input files
# note that sys.argv has the script itself in the 0th position 
for i in range(1, noFiles):
    # open each file and read lines 
    fileHandle = open(sys.argv[i], "r")
    lines = fileHandle.readlines()
    
    # assume the input is a single line
    # of space-separated numbers
    genes = lines[0].split()
    
    # make a dict with the file name and the genes from that file
    # this will later be used for assigning each gene to the right genome
    genomes[(sys.argv[i].split("_")[1])] = genes

    if len(geneLists)<1:
        geneLists = genes # save this list
    else:
        # if the list is defined then we select the overlapping genes from the previous files and the new file
        geneLists = list(set(geneLists)&set(genes))

    # close the file so not weird things happens to the data
    fileHandle.close()


# assign new indexing from 1 to end
id = 1
newMap = {}
for gene in geneLists:
    newMap[gene] = str(id)
    id = id + 1

# print out the revised gene lists
for genome in genomes:
    presentGenes = list()
    for gene in genomes[genome]:
        if gene in newMap and newMap[gene] not in presentGenes:
            presentGenes.append(newMap[gene])

    print (">Genome_"+genome+"\n"+" ".join(presentGenes))
