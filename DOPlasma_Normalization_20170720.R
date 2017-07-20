########## DO Mice: Plasma Normalization ##########
################################################################################
# DO Plasma GC-MS Metabolite Data.
# Samples ran.
# Anji Trujillo
# etrujillo2@wisc.edu
# July 20, 2017
################################################################################
########## Loading packages ############
library(pheatmap)
#library(preprocessCore)
#library(sva)
library(ggplot2)
#library(gplots)
library(ggfortify)
#library(caret)
library(nortest)
library(qqtest)
library(reshape2)
#library(csvread)


#####################
# Working Directory #
#####################

getwd() #[1] "C:/Users/etrujillo/Documents"

setwd("C:/Users/etrujillo/Desktop/DOProjectFolder")

#####################
# Load in the data. #
#####################

DOPlasma <- read.csv("19July2017DOPlasmaMetabolites_RAW_NoTransformation_Tier5_ForR.csv", 
                     header = TRUE, sep = ",", stringsAsFactors = FALSE) # Read in csv and skip first line with batch info
DOPlasma <- DOPlasma[,-c(3)] # Remove the third columns
DOPlasmaNoRT <- DOPlasma[,-c(1:2)] #No Retention Time Info
DOPlasma_Header <- read.csv("19July2017DOPlasmaMetabolites_RAW_NoTransformation_Tier5_ForR.csv", 
                            header = FALSE, sep = ",", stringsAsFactors = FALSE) #read in a csv that includes the batch information

#####################
# Explore the data. #
#####################

str(DOPlasma)
summary(DOPlasma)

###########################################
# Restructure pivot table into long form. #
###########################################

molten.DOPlasma <- melt(DOPlasma, variable.name = "Mouse.ID", 
                        value.names = "Intensity", id.vars = c("RetentionTime", "Feature.ID") ) #melt the data into long form, renaming the columns appropriately
BatchDF <- t(DOPlasma_Header[1:2, 4:ncol(DOPlasma_Header)]) # create a two lists with Batch and MouseID for later use in loop
colnames(BatchDF)[1] <- "Batch.ID" #rename column names in BatchDF 
colnames(BatchDF)[2] <- "Mouse.ID"
Batch_vector <- NULL # create an empty batch vector

#For loop for parsing through the long form DO Plasma data and appending the batch information to the correct mouse sample
for(i in 1:nrow(molten.DOPlasma)){ 
  ID = molten.DOPlasma[i, 3]
  batch <- BatchDF[grep(ID, BatchDF[,2])[1],1]
  #ifelse(is.na(molten.DOPlasma[i,2]), "", molten.DOPlasma[i,2] <- batch)
  Batch_vector <- append(Batch_vector,batch)
  }

#bind vector and long form Plasma data  
Long_DOPlasma <- cbind(molten.DOPlasma, Batch_vector)

##################
# Plotting data. #
##################

Transposed_DOPlasma <- t(DOPlasma[,-c(2)]) # Transpose plasma data frame so that the rows are samples. Correct format for prcomp PCA function
colnames(Transposed_DOPlasma) <- Transposed_DOPlasma[1,] # rename column names
autoplot(prcomp(Transposed_DOPlasma[-c(1),])) #Plot PCA
#autoplot(prcomp(Long_DOPlasma[,c(2,4)]), data = Long_DOPlasma, color = 'Batch_vector', label = TRUE, shape = FALSE, label.size = 3)

plot(Transposed_DOPlasma, type = "l")
pca_plasma <- prcomp(Transposed_DOPlasma[-c(1),])
  

plot(density(t(DOPlasmaNoRT))) # Plot density

metaboliteID <- NULL
  
  metaboliteID <- for(i in ID){print(rep(1:245))}
  
molten.DOPlasma <- cbind(molten.DOPlasma, metaboliteID)  
  


autoplot(prcomp(DOPlasma), label = TRUE, loadings = TRUE)

autoplot(princomp(mydf2))

head(mydf2[1:3])

mydf2[2] <- lapply(mydf2[2], as.numeric)


ggplot(mydf, aes(x=Intensity, y=metabolite)) +
  geom_point(colour="lightblue", alpha=0.1, position="jitter") +
  geom_boxplot(alpha=0.2)

ggplot(mydf, aes(x=Intensity, fill=metabolite)) +
  geom_histogram(position="identity", binwidth=0.25, alpha=0.3)

barplot(t(all.dataMatrix2), beside = TRUE)
barplot(t(cecumGroup), beside = TRUE)

ggplot(mydf, aes(x=metabolite, y=Intensity)) +
  geom_point(position="jitter", alpha=0.3) +
  geom_boxplot(outlier.size=1, alpha=0.2)

ggplot(mydf2, aes(x=metabolite, y=Intensity)) +
  geom_point(position="jitter", alpha=0.3) +
  geom_boxplot(outlier.size=1, alpha=0.2)

boxplot(t(numericMetaboliteCecum))

boxplot(t(all.dataMatrix2))

Transposed.all.dataMatrix2 <- t(all.dataMatrix2)

write.table(Transposed.all.dataMatrix2, "~Transposed.all.dataMatrix2.txt", sep="\t")

barchart(t(cecumGroup))

barchart(all.dataMatrix2)

dotplot(cecumGroup)

dotplot(all.dataMatrix2)

boxplot(t(cecumGroup))

boxplot(t(MetabolitesPlasma), notch = FALSE)
boxplot(t(MetabolitesCecum))


barplot(t(cecumGroup), beside = T)

metabolite.df = data.frame(
  all.dataMatrix2,
  Biofluid = factor(c(rep(Cecum, 10), rep(Plasma, 9))))

# function provides a summary of centeral tendencies and distributions
transposedMetabolitesCecum <- t(MetabolitesCecum)
df1=stack(MetabolitesCecum)
mydf=data.frame(Intensity=df1$values,Mouse=df1$ind,metabolite=sampleNamesCecum)
mydf$Intensity <- as.numeric(as.character(mydf$Intensity))

df2=stack(MetabolitesPlasma)
mydf2=data.frame(Intensity=df2$values,Mouse=df2$ind,metabolite=sampleNames)
mydf2$Intensity <- as.numeric(as.character(mydf2$Intensity))

mydf$Intensity <- as.numeric(as.character(mydf$Intensity))


dotplot(Mouse ~ Intensity, data = mydf, groups = metabolite,
       main = "Irrigation Area by Region")

boxplot(Mouse ~ Intensity, data = mydf,
        xlab = "Month", ylab = "Maximum Temperature",
        main = "Temperature at Southampton Weather Station (1950-1999)"
)

######################## Use the sapply function to apply a function to every element in your dataframe

# Create a variable that saves the summary for each column
PeakHeightSummary <- sapply(PeakHeight, summary.list )

firstBatchPeakHeightSummary <- sapply(firstBatchPeakHeight[-c(30:36),], summary.list ) # on entire list of founderStrain controls
secondBatchPeakHeightSummary <- sapply(secondBatchPeakHeight[-c(25:30),], summary.list ) # on entire list of founderStrain controls

######################## Look at the distribution summaries of the Controls

# Use the grep function to separate the all the controls from my data
PeakHeight.Control <- PeakHeight[grep("_Control", rownames(PeakHeight)),]

# index number for the samples with the following pattern
grep("_Control", rownames(firstBatchPeakHeight)) # Controls in batch 1 are index 34 35 36
grep("_Blank", rownames(firstBatchPeakHeight)) # Blanks in batch 1 are index 30 31 32 33






