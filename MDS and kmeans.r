install.packages('bios2mds')
library('bios2mds')
# Your working directory
setwd('foo')

# File containing a list of sequence names in same order of the matrix
seqNames = read.table(file='seqNames', sep = ',')
# csv file that is the matrix of genetic distances
myTable = read.table(file='matrixCSV', sep = ',')
colnames(myTable) = seqNames$V1
rownames(myTable) = seqNames$V1
myMatrix = as.matrix(myTable)
# File containing a csv of the species of each of the seqs, again in same order as the name file and the matrix
speciesInfo = read.table(file='seqSpeciesInfo', sep=',')
colnames(speciesInfo) = c('Seq', 'Species')

# Perform the MDS
myMDS = mmds(myMatrix, pc=2)
MDSCOORDS = myMDS$coord

# Produce the data frame from which to make the plots
plotDetails = cbind(MDSCOORDS, speciesInfo)

# Perform the kmeans clustering for the plot data
set.seed = 20
myCluster = kmeans(plotDetails[,1:2], 3, nstart=50)                                                                                 

#If you want to produce a shoulder plot
wss = (nrow(plotDetails[,1:2])-1)*sum(apply(plotDetails[,1:2], 2, var))
for (i in 2:15) wss[i] <- kmeans(plotDetails[,1:2], i, nstart=50)$tot.withinss
plot(1:15, wss, type='b', xlab = 'number of clusters', ylab = 'within groups sum of squares')
#...

newPlotDetails = cbind(plotDetails, as.factor(myCluster$cluster))

# Write out this data as csv for use in chi squared Python script part of 
# which will calculate the centroids for plotting and calculating the SEMs
write.table(newPlotDetails, file = "newPlotDetails.txt", quote=FALSE, row.names=FALSE, col.names = FALSE, sep=',')


# THese are the centroids for the speceis clusters
# Calculated as part of the Chi squared Python script
xvals = c(0.0199411, -0.02711428, 0.03194736)
yvals = c(-0.11282, -0.0748857, 0.2387894)
speciesCenters = cbind(xvals, yvals)

# Purely for visualising the clustering
plot(plotDetails$PC1, plotDetails$PC2, col=c('red', 'green3', 'darkblue')[as.numeric(newPlotDetails$`as.factor(myCluster$cluster)`)], pch=c(0, 1, 2)[as.numeric(newPlotDetails$`as.factor(myCluster$cluster)`)])

#This is for the figure Plot saved as SVG
svg("specDistrFig.svg")
plot(plotDetails$PC1, plotDetails$PC2, col=c('darkblue', 'green3', 'red')[as.numeric(newPlotDetails$Species)], pch=c(2, 1, 0)[as.numeric(newPlotDetails$Species)])
points(speciesCenters, col=c('green3', 'red', 'darkblue'), bg='white', pch=c(21, 22, 24), cex=2)
points(speciesCenters, col=c('green3', 'red', 'darkblue'),  pch=3, cex=8)
points(myCluster$centers, col='black', bg='grey', pch=21, cex=1)
points(myCluster$centers, col='black', pch=3, cex=8)
dev.off()
