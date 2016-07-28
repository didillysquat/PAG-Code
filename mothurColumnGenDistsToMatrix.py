def readDefinedFileToList(filename):
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList

def createMatrixFromMothurColOutput():
    # Read in the Mothur output distances in column format
    rootDir = r'SomeDir'
    distColFile = readDefinedFileToList(rootDir + '\distances.dist')
    distDict = {}
    seqList = []
    # Create a dictionary that has the pair of sequences as a key
    # (frozen set so that order doesn't matter)
    # and genetic distance as value
    # Then to make a list of the sequences, add any sequence names
    # that aren't already in the list, in
    for row in distColFile:
        # First put the entry into the dict
        distDict[frozenset(row.split(' ')[:-1])] = float(row.split(' ')[2])
        # Then add any new seq names to list
        for seq in row.split(' ')[:-1]:
            if seq not in seqList:
                seqList.append(seq)

    # At this point we have a list of seqs from the Mothur output
    # and a dict with seq pairs as key and gen distance as value
    # Create empty matrix that is of dimensions
    # [rows(n seqs + 1)][columns(number of seqs + 1)]
    distMatrix = [list(seqList)]
    distMatrix[0].insert(0, 'Samples')
    # Populate the matrix with 'N/A' for the time being
    for sample in seqList:
        blankRow = ['N/A']*len(seqList)
        blankRow.insert(0, sample)
        distMatrix.append(blankRow)

    # Matrix is FstMatrix[row][column]
    # Fill matrix from the column genetic distances in FstColDist
    row = 1
    while row < len(distMatrix): # For each row
        col = 1
        while col < len(distMatrix[row]): # For each col
            sampleOne = distMatrix[0][col]
            sampleTwo = distMatrix[row][0]
            if sampleOne == sampleTwo:
                distMatrix[row][col] = 0.00
            else:
                    distMatrix[row][col] = distMatrix[col][row] = \
                        distDict[frozenset([sampleOne, sampleTwo])]
            col += 1
        row += 1

    # Get rid of the column and row headers to leave just pure
    # numerical matrix
    del distMatrix[0]
    distMatrix = [a[1:] for a in distMatrix]

    # Write out the matrix and the ordered list of seqs
    write2DListToDestination(r'someFileDirectory', distMatrix)
    writeListToDestination(r'someFileDirectory', seqList)

createMatrixFromMothurColOutput()