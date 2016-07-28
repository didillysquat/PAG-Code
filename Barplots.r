# This code is to make a 'grouped bar chart'
# There are a few pitfalls worth knowing about
# Namely, you don't have to use the 'table()' function unless you have your data uncounted e.g. one row for one seq would be POB, SAD, rDNA and you would have lots of rows like that
# The better way is to already have done your counts so your table is a data frame and looks like:
# lobata  lutea harri site
# 0       0     4     DAL
# 4       9     4     SAD
#etc.
# THen to make the bar plot simply set the height to one of the variables. e.g. data$lobata, and the names.arg to the sites.
# However, if you want the sites to be in a particular order you have to set them as a factor. The easiest way to do this is to add a new variable that is a factor to the dataset
# e.g. in the above example you would do data$orderedsites= factor(data$sites, levels=c('DAL', 'SAD', 'RAG', 'UMM', 'RAK', 'MW', 'ME', 'FUJ', 'MUC'))
# THen, and only then when you make your barplots you can choose data$orderedsites as your names.arg.
# Also when you want to plot several variables, instead of just putting in one e.g. data$loabata
# You put in barplot(rbind(data$lobata, data$lutea, etc.), names... etc.)
# You will see this in action below. You are welcome.
# http://onertipaday.blogspot.co.uk/2007/05/barplots-of-two-sets.html
# http://stackoverflow.com/questions/24247460/r-ordering-a-barplot-histogram-generated-from-a-data-frame

setwd('foo')
# A file in the format
# P. lobata rDNA,P. lutea rDNA,P. harrisoni rDNA, P. lobata psbA, etc..., Sites(this is the final column)
# XXX,XXX,XXX,XXX, etc.
#Where each row is a different location in the order seen below e.g. DAL, SAD, etc.
# The final column is sites and is what you see as counts2$sites below
counts2 = read.csv('counts2.csv')
counts2$orderedsites = factor(counts2$sites, levels = c('DAL', 'SAD', 'RAG', 'UMM', 'RAK', 'MW', 'ME', 'FUJ', 'MUC'))
# Produce the plots
svg('SpDistFiginset.svg')
par(mfrow=c(1,2))
# Create the rDNA plots
# Plot either one of the two below
barplot(rbind(counts2$P..lobata.rDNA, counts2$P..lutea.rDNA, counts2$P..harrisoni.RDNA, counts2$Porites.spp...unidentified..rDNA), names.arg=counts2$orderedsites, las=3, beside=TRUE, col = c('green3', 'red', 'darkblue', 'grey'), ylim=c(0, 20))
barplot(rbind(counts2$P..lobata.rDNAseqs, counts2$P..lutea.rDNAseqs, counts2$P..harrisoni.rDNAseqs, counts2$POR.rDNAseqs), names.arg=counts2$orderedsites, las=3, beside=TRUE, col = c('green3', 'red', 'darkblue', 'grey'), ylim=c(0, 80))
# Create the psbA plots
barplot(rbind(counts2$P..lobata.psba, counts2$P..lutea.psba, counts2$P..harrisoni.psba, counts2$Porites.spp...unidentified..psba), names.arg=counts2$orderedsites, las=3, beside=TRUE, col = c('green3', 'red', 'darkblue', 'grey'), ylim=c(0, 20))

dev.off()
