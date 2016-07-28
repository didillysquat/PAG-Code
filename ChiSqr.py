# This script will take in the new plot details csv which contains information for each sequence and which
# cluster it fell within according to being clustered by the kmeans
# that were conducted in R.
# The distanceMatrixScript.py created the data that was required for the Kmeans clustering.

def readDefinedFileToList(filename):
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList


def main():
    rootDir = r'C:\Users\User\Google Drive\BS work\Figures\MDS rDNA seqs'
    dataFile = readDefinedFileToList(rootDir + r'\newPlotDetails.txt')

    # We want to produce a table that looks like

    #   cluster 1   cluster 2   cluster 3
    #POB
    #PLT
    #POH

    # So we will count each seq and add it to one of the 9 cell counts
    # Numbered from top left to bottom right
    cellCounts = [0]*9

    for line in dataFile:
        if line.split(',')[3] == 'POB':
            if int(line.split(',')[4]) == 1:
                cellCounts[0] += 1
            elif int(line.split(',')[4]) == 2:
                cellCounts[1] += 1
            elif int(line.split(',')[4]) == 3:
                cellCounts[2] += 1
        elif line.split(',')[3] == 'PLT':
            if int(line.split(',')[4]) == 1:
                cellCounts[3] += 1
            elif int(line.split(',')[4]) == 2:
                cellCounts[4] += 1
            elif int(line.split(',')[4]) == 3:
                cellCounts[5] += 1
        elif line.split(',')[3] == 'POH':
            if int(line.split(',')[4]) == 1:
                cellCounts[6] += 1
            elif int(line.split(',')[4]) == 2:
                cellCounts[7] += 1
            elif int(line.split(',')[4]) == 3:
                cellCounts[8] += 1
    # This prints out the observed table
    print('\tCluster1\tCluster2\tCluster3')
    print('POB\t' + '\t'.join([str(a) for a in cellCounts[0:3]]) + '\t' + str(sum([a for a in cellCounts[0:3]])))
    print('PLT\t' + '\t'.join([str(a) for a in cellCounts[3:6]]) + '\t' + str(sum([a for a in cellCounts[3:6]])))
    print('POH\t' + '\t'.join([str(a) for a in cellCounts[6:]])+ '\t' + str(sum([a for a in cellCounts[6:]])))
    print('\t' + '\t'.join([str(a) for a in [sum([cellCounts[0], cellCounts[3], cellCounts[6]]), sum([cellCounts[1], cellCounts[4], cellCounts[7]]), sum([cellCounts[2], cellCounts[5], cellCounts[8]])]]))
    print('Total seqs = ' + str(sum(cellCounts)))
    # Now we need to calculate the expected table
    #http://www.stat.yale.edu/Courses/1997-98/101/chisq.htm
    # Each expected cell is calculated as (row total*column total)/total number of observations
    cellCountsExpected = [0]*9
    i = 0
    while i < len(cellCountsExpected):
        if (i+1)/3 <= 1: # Then we are in the first row
            if (i+1)%3 == 1: # Then we are in the first col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[0:3]])*sum([cellCounts[0], cellCounts[3], cellCounts[6]]))/sum(cellCounts))
            elif (i+1)%3 == 2: # Then we are in the second col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[0:3]]) * sum([cellCounts[1], cellCounts[4], cellCounts[7]])) / sum(cellCounts))
            elif (i+1)%3 == 0: # Then we are in the third col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[0:3]]) * sum([cellCounts[2], cellCounts[5], cellCounts[8]])) / sum(cellCounts))
        elif (i+1)/3 <= 2: # Then we are in the second row
            if (i + 1) % 3 == 1:  # Then we are in the first col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[0], cellCounts[3], cellCounts[6]])) / sum(cellCounts))
            elif (i + 1) % 3 == 2:  # Then we are in the second col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[1], cellCounts[4], cellCounts[7]])) / sum(cellCounts))
            elif (i + 1) % 3 == 0:  # Then we are in the third col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[2], cellCounts[5], cellCounts[8]])) / sum(cellCounts))
        elif (i+1)/3 <= 3: # Then we are in the third row
            if (i + 1) % 3 == 1:  # Then we are in the first col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[0], cellCounts[3], cellCounts[6]])) / sum(cellCounts))
            elif (i + 1) % 3 == 2:  # Then we are in the second col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[1], cellCounts[4], cellCounts[7]])) / sum(cellCounts))
            elif (i + 1) % 3 == 0:  # Then we are in the third col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[2], cellCounts[5], cellCounts[8]])) / sum(cellCounts))
        i += 1
# Now we print the expected table
    print('EXPECTED')
    print('\tCluster1\tCluster2\tCluster3')
    print('POB\t' + '\t'.join([str(a) for a in cellCountsExpected[0:3]]) )
    print('PLT\t' + '\t'.join([str(a) for a in cellCountsExpected[3:6]]) )
    print('POH\t' + '\t'.join([str(a) for a in cellCountsExpected[6:]]))
    print('Total seqs = ' + str(sum(cellCountsExpected)))

    # Now we compute the chi squared value which is the sum of ((o-e)^2)/e
    chiSqrResult = []
    i = 0
    while i < len(cellCounts):
        chiSqrResult.append(((cellCounts[i]-cellCountsExpected[i])**2)/cellCountsExpected[i])
        i += 1
    print('differences')
    print('\t'.join([str(a) for a in chiSqrResult[:3]]))
    print('\t'.join([str(a) for a in chiSqrResult[3:6]]))
    print('\t'.join([str(a) for a in chiSqrResult[6:]]))
    print('The Chi Sqr value is: ' + str(sum(chiSqrResult)))

    # Critical chi2 value for 4 df at alpha = 0.05 is 9.488 ours came to 5.94 so we are good! Cannot reject null hypothesis. So no significant effect of speceies on clustering exists
# main()

def mainpsba():
    rootDir = r'C:\Users\User\Google Drive\BS work\Figures\MDS psba'
    dataFile = readDefinedFileToList(rootDir + r'\newPlotDetails.txt')

    # We want to produce a table that looks like

    #   cluster 1   cluster 2   cluster 3
    #POB
    #PLT
    #POH

    # So we will count each seq and add it to one of the 9 cell counts
    # Numbered from top left to bottom right
    cellCounts = [0]*9

    for line in dataFile:
        if line.split(',')[3] == 'P. lobata':
            if int(line.split(',')[4]) == 1:
                cellCounts[0] += 1
            elif int(line.split(',')[4]) == 2:
                cellCounts[1] += 1
            elif int(line.split(',')[4]) == 3:
                cellCounts[2] += 1
        elif line.split(',')[3] == 'P. lutea':
            if int(line.split(',')[4]) == 1:
                cellCounts[3] += 1
            elif int(line.split(',')[4]) == 2:
                cellCounts[4] += 1
            elif int(line.split(',')[4]) == 3:
                cellCounts[5] += 1
        elif line.split(',')[3] == 'P. harrisoni':
            if int(line.split(',')[4]) == 1:
                cellCounts[6] += 1
            elif int(line.split(',')[4]) == 2:
                cellCounts[7] += 1
            elif int(line.split(',')[4]) == 3:
                cellCounts[8] += 1
    # This prints out the observed table
    print('\tCluster1\tCluster2\tCluster3')
    print('POB\t' + '\t'.join([str(a) for a in cellCounts[0:3]]) + '\t' + str(sum([a for a in cellCounts[0:3]])))
    print('PLT\t' + '\t'.join([str(a) for a in cellCounts[3:6]]) + '\t' + str(sum([a for a in cellCounts[3:6]])))
    print('POH\t' + '\t'.join([str(a) for a in cellCounts[6:]])+ '\t' + str(sum([a for a in cellCounts[6:]])))
    print('\t' + '\t'.join([str(a) for a in [sum([cellCounts[0], cellCounts[3], cellCounts[6]]), sum([cellCounts[1], cellCounts[4], cellCounts[7]]), sum([cellCounts[2], cellCounts[5], cellCounts[8]])]]))
    print('Total seqs = ' + str(sum(cellCounts)))
    # Now we need to calculate the expected table
    #http://www.stat.yale.edu/Courses/1997-98/101/chisq.htm
    # Each expected cell is calculated as (row total*column total)/total number of observations
    cellCountsExpected = [0]*9
    i = 0
    while i < len(cellCountsExpected):
        if (i+1)/3 <= 1: # Then we are in the first row
            if (i+1)%3 == 1: # Then we are in the first col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[0:3]])*sum([cellCounts[0], cellCounts[3], cellCounts[6]]))/sum(cellCounts))
            elif (i+1)%3 == 2: # Then we are in the second col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[0:3]]) * sum([cellCounts[1], cellCounts[4], cellCounts[7]])) / sum(cellCounts))
            elif (i+1)%3 == 0: # Then we are in the third col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[0:3]]) * sum([cellCounts[2], cellCounts[5], cellCounts[8]])) / sum(cellCounts))
        elif (i+1)/3 <= 2: # Then we are in the second row
            if (i + 1) % 3 == 1:  # Then we are in the first col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[0], cellCounts[3], cellCounts[6]])) / sum(cellCounts))
            elif (i + 1) % 3 == 2:  # Then we are in the second col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[1], cellCounts[4], cellCounts[7]])) / sum(cellCounts))
            elif (i + 1) % 3 == 0:  # Then we are in the third col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[2], cellCounts[5], cellCounts[8]])) / sum(cellCounts))
        elif (i+1)/3 <= 3: # Then we are in the third row
            if (i + 1) % 3 == 1:  # Then we are in the first col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[0], cellCounts[3], cellCounts[6]])) / sum(cellCounts))
            elif (i + 1) % 3 == 2:  # Then we are in the second col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[1], cellCounts[4], cellCounts[7]])) / sum(cellCounts))
            elif (i + 1) % 3 == 0:  # Then we are in the third col
                cellCountsExpected[i] = int((sum([a for a in cellCounts[3:6]]) * sum([cellCounts[2], cellCounts[5], cellCounts[8]])) / sum(cellCounts))
        i += 1
# Now we print the expected table
    print('EXPECTED')
    print('\tCluster1\tCluster2\tCluster3')
    print('POB\t' + '\t'.join([str(a) for a in cellCountsExpected[0:3]]) )
    print('PLT\t' + '\t'.join([str(a) for a in cellCountsExpected[3:6]]) )
    print('POH\t' + '\t'.join([str(a) for a in cellCountsExpected[6:]]))
    print('Total seqs = ' + str(sum(cellCountsExpected)))

    # Now we compute the chi squared value which is the sum of ((o-e)^2)/e
    chiSqrResult = 0
    differences = []
    i = 0
    while i < len(cellCounts):
        chiSqrResult += (((cellCounts[i]-cellCountsExpected[i])**2)/cellCountsExpected[i])
        differences.append(str((((cellCounts[i]-cellCountsExpected[i])**2)/cellCountsExpected[i])))
        i += 1
    print(chiSqrResult)
    print('\t'.join(differences[0:3]))
    print('\t'.join(differences[3:6]))
    print('\t'.join(differences[6:]))

main()