import math

def readDefinedFileToList(filename):
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList

def calculateGlobalErrorsOfTheMean():
    # Read in a file that has the information on which cluster each of the sequences fell into according
    # to the kmeans clustering. The newPlotDetails file is produced by the MDS and kmeans .r file
    rootDir = r'yourRootDir'
    dataFile = readDefinedFileToList(rootDir + r'\newPlotDetails.txt')
    # Go through each seq (line) and store the x and y values dependent on which cluster the points fall into within either the species of K-means strategy for clustering.
    # With this we will work out the centroids (the average x and average y coordinates) for each of the clusters for each of the clustering strategies (i.e. kMeans, vs. species)
    speciesCluster = [[],[],[]] # list order POB, PLT POH
    kMeansCluster = [[],[],[]] # list order cluster 1, 2 and 3
    for seq in dataFile:
        # Do species cluster first
        if seq.split(',')[3] == 'POB':
            speciesCluster[0].append([float(seq.split(',')[0]), float(seq.split(',')[1])])
        elif seq.split(',')[3] == 'PLT':
            speciesCluster[1].append([float(seq.split(',')[0]), float(seq.split(',')[1])])
        elif seq.split(',')[3] == 'POH':
            speciesCluster[2].append([float(seq.split(',')[0]), float(seq.split(',')[1])])
        # Do kmeans cluster second
        if seq.split(',')[4] == '1':
            kMeansCluster[0].append([float(seq.split(',')[0]), float(seq.split(',')[1])])
        elif seq.split(',')[4] == '2':
            kMeansCluster[1].append([float(seq.split(',')[0]), float(seq.split(',')[1])])
        elif seq.split(',')[4] == '3':
            kMeansCluster[2].append([float(seq.split(',')[0]), float(seq.split(',')[1])])

    # At this point we should have the x and y coords for all of the seqs put in the appropriate lists
    # Run through and work out the x and y means for each of the clusters for each of the clustering strategies

    # Species cluster first
    i = 0
    speciesClusterMeansList = []
    while i < len(speciesCluster):
        xtotal = 0
        ytotal = 0
        count = 0
        for seq in speciesCluster[i]:
            count += 1
            xtotal = xtotal + seq[0]
            ytotal = ytotal + seq[1]
        speciesClusterMeansList.append([float(xtotal/count), float(ytotal/count)])
        i += 1

    # Kmeans cluster second
    i = 0
    kMeansClusterMeansList = []
    while i < len(kMeansCluster):
        xtotal = 0
        ytotal = 0
        count = 0
        for seq in kMeansCluster[i]:
            count += 1
            xtotal = xtotal + seq[0]
            ytotal = ytotal + seq[1]
        kMeansClusterMeansList.append([float(xtotal/count), float(ytotal/count)])
        i += 1

    # Out put the centroids for species clustering
    print('Cluster means for the species method of clustering are:')
    print(speciesClusterMeansList[0])
    print(speciesClusterMeansList[1])
    print(speciesClusterMeansList[2])
    # Output the centroids for kmeans clustering - these are given in the r output when the kmeans
    # clustering is done (cluster$centers), so you can verify that the values agree.
    print('Cluster means for the kmeans method of clustering are:')
    print(kMeansClusterMeansList[0])
    print(kMeansClusterMeansList[1])
    print(kMeansClusterMeansList[2])

    # At this point we know each of the means for each of the clusters for each of the cluster strategies
    # Now for each cluster of each strategy, calculate the standard error of the mean

    # Species clustering first
    i = 0
    speciesClusterStandardErrors = []
    while i < len(speciesCluster): # For each of the clusters e.g. POB, PLT, POH
        sumOfSquaredDeviations = 0
        count = 0
        for seq in speciesCluster[i]:
            count += 1
            hypotenuse = math.sqrt(((seq[0] - speciesClusterMeansList[i][0])**2) + ((seq[1] - speciesClusterMeansList[i][1])**2))
            sqrHypotenuse = hypotenuse**2
            sumOfSquaredDeviations += sqrHypotenuse

        # At this point we should have the sum of squared deviations
        SD = math.sqrt(sumOfSquaredDeviations/(count - 1))
        SE = SD/(math.sqrt(count))
        speciesClusterStandardErrors.append(SE)
        i += 1

    # At this point we have the three SE for the species clustering
    globalSEForSpeciesClustering = sum(speciesClusterStandardErrors)/3

    # Kmeans clustering second
    i = 0
    kmeansClusterStandardErrors = []
    while i < len(kMeansCluster): # For each of the clusters e.g. POB, PLT, POH
        sumOfSquaredDeviations = 0
        count = 0
        for seq in kMeansCluster[i]:
            count += 1
            hypotenuse = math.sqrt(((seq[0] - kMeansClusterMeansList[i][0])**2) + ((seq[1] - kMeansClusterMeansList[i][1])**2))
            sqrHypotenuse = hypotenuse**2
            sumOfSquaredDeviations += sqrHypotenuse

        # At this point we should have the sum of squared deviations
        SD = math.sqrt(sumOfSquaredDeviations/(count - 1))
        SE = SD/(math.sqrt(count))
        kmeansClusterStandardErrors.append(SE)
        i += 1

    # At this point we have the three SEs for the kmeans clustering
    globalSEForKMeansClustering = sum(kmeansClusterStandardErrors)/3

    # At this point we have both of the global standard errors for both clustering types
    # Print them
    print('The global standard error of the mean for the species clustering approach is: ' + str(globalSEForSpeciesClustering))
    print('The global standard error of the mean for the species clustering approach is: ' + str(globalSEForKMeansClustering))