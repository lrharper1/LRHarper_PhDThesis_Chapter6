#' ---
#' Title: "Mammal eDNA metabarcoding data analysis"
#' Author: "Lynsey Rebecca Harper"
#' Date: "20th August 2018"
#' ---
#' 
#' 
#' eDNA metabarcoding has shown potential to detect semi-aquatic and terrestrial
#' wildlife as well as aquatic species. We conducted two experiments to validate
#' this molecular approach as a monitoring tool for mammals that are conservation,
#' management, and reintroduction priority species in the UK.
#' 
#' The first experiment was based at two wildlife parks that house these mammal
#' species and designed to understand how mammal behaviour influences eDNA 
#' detection under controlled conditions.
#' 
#' The second experiment was designed to validate the eDNA metabarcoding approach
#' under natural conditions. eDNA metabarcoding was compared to camera trapping 
#' at two ponds each from 3 locations across the UK where our target species are
#' confirmed as present.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. Then,
#' load functions for calculating kappa coefficient, plotting model 
#' residuals and testing model fit.
#' 

## Clear memory
rm(list=ls())

## Direct R to folder containing packages on workstation
#.libPaths(c("C:\\R\\rlib", .libPaths("rlib")))

## set working directory to the location of the script
# install.packages("rstudioapi") # first time for each computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Check working directory
getwd()

## Load required packages
p <- c("ggplot2","ggpubr","lemon","munsell","lazyeval","grid","gridExtra",
       "forcats", "coda","MASS","car","scales","AICcmodavg","xtable",
       "gtools","xlsx","reshape2","plyr","dplyr","tidyr","arm","RVAideMemoire",
       "permute","ResourceSelection","bbmle","RColorBrewer", "MuMIn",
       "ggmap","mapproj","geosphere","jpeg","proto","rjson","RgoogleMaps",
       "maps","labeling","ggsn","png","coin","modeltools","mvtnorm","pROC",
       "vegan","multcomp","betapart","adespatial","adegraphics","ade4","gdata")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/",
                                          dependencies=TRUE)
lapply(p, require, character.only = TRUE)
devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
library(glmmTMB)

## Load custom functions
f <- c("CheckResidsFunction.R", "OverdispersalFunction.R", 
       "CheckConvergenceFunction.R","glmmTMB_posthoc_test_function.R")
lapply(f, source)

#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()


#' ---
#' 
#' ## 1) Raw data processing
#' 
#' Examine the raw data output from metaBEAT.
#' 

## Original BLAST 
ass.raw <- read.csv("../Data/Mammal_metabarcoding_assigned_raw.csv", header=TRUE)
summary(ass.raw[,c(1:5,length(ass.raw))])
head(ass.raw)
names(ass.raw)
str(ass.raw)

## Unassigned BLAST
unass.raw <- read.csv("../Data/Mammal_metabarcoding_unassigned_raw.csv", header=TRUE)
summary(unass.raw[,c(1:5,length(unass.raw))])
head(unass.raw)
names(unass.raw)
str(unass.raw)

## Make first row header names
names(ass.raw) <- lapply(ass.raw[1,], as.character)
names(unass.raw) <- lapply(unass.raw[1,], as.character)
ass.raw <- ass.raw[-1,]
unass.raw <- unass.raw[-1,]

## Reset row names
rownames(ass.raw) <- NULL
rownames(unass.raw) <- NULL

## Remove '-nc.blast'
colnames(ass.raw) <- gsub('-nc.blast', '', colnames(ass.raw))
colnames(unass.raw) <- gsub('-nc.blast.blast', '', colnames(unass.raw))

## Rename first column
colnames(ass.raw)[1] <- "Assignment"
colnames(unass.raw)[1] <- "Assignment"

## Create new dataframe for calculating the proportional read counts 
## for each mammal species downstream.
## Duplicate raw read data
raw1 <- ass.raw
raw2 <- unass.raw

## Remove last column containing taxonomy
raw1 <- raw1[,-327]
raw2 <- raw2[,-327]

## Bind data frames
raw.merged <- rbind(raw1, raw2)

## Merge read counts by taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
raw.merged[,2:326] <- lapply(raw.merged[,2:326], function(x) as.numeric(as.character(x)))
raw.merged <- ddply(raw.merged, .(Assignment), numcolwise(sum))

## Create copy of dataframe without positive or negative controls
raw.df <- raw.merged[,which(!grepl("CB|FB|EB|Negative|Positive", 
                                   colnames(raw.merged)))]

## Make Assignment column row names
rownames(raw.df) <- raw.df$Assignment
raw.df <- raw.df[,-1]

## Calculate total number of reads in samples
raw.df <- rbind(raw.df, colSums(raw.df))

## Make new dataframe containing sample ID and total number of reads 
## per sample
raw.total <- raw.df[204,]
raw.total$Assignment <- "Total"
raw.total <- raw.total[,c(221,1:220)]
rownames(raw.total) <- NULL

## Remove any taxonomic assignments that aren't vertebrate from each dataframe
## using the taxonomy column created during processing with metaBEAT
ass <- ass.raw[which(grepl("Chordata|unassigned", ass.raw$taxomomy)),]
unass <- unass.raw[which(grepl("Chordata|unassigned", unass.raw$taxomomy)),]

## Remove last column containing taxonomy
ass <- ass[-327]
unass <- unass[-327]

## Bind data frames
merged.df <- rbind(ass, unass)

## Merge read counts by taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
merged.df[,2:326] <- lapply(merged.df[,2:326], function(x) as.numeric(as.character(x)))
merged.df <- ddply(merged.df, .(Assignment), numcolwise(sum))

## Export as .csv file
write.csv(merged.df, "../Data/Mammal_metabarcoding_merged.csv", row.names=FALSE)


#' --- 
#' 
#' ## 2) Refine dataset
#' 
#' Now, we need to further refine the metabarcoding dataset.
#' 
#' 1. Any spurious species must be removed. Use NBN atlas to check species
#'    occurrence records and ensure they match with sampling locations. This
#'    is also a good source for checking current taxonomy.
#' 2. Any Genus or Family assignments containing only one species in 
#'    the UK must be changed to that species.
#' 3. Likewise, for any species assignment which is the only species in 
#'    the UK and also has genus/family reads, read counts from all
#'    assignments must be merged.
#'  
#' A record of these changes will be kept as these changes are made.
#'      

## Inspect new dataframe
summary(merged.df[,1:6])
names(merged.df)
str(merged.df)

#' 
#' First, remove spurious assignments.
#' 
#' - *Pinicola enucleator* pine grosbeak, row 77
#' - *Schoeniclus rusticus* rustic bunting, row 90
#' - *Atherina boyeri* big-scale sand smelt, row 106
#' - *Bos mutus* domestic yak, row 110
#' - Chordata_environmental_sample, row 114
#' - *Pan*, row 127
#' - *Sardina pilchardus*, row 131
#' 

true.assign <- merged.df[-c(77,90,106,110,114,127,131),]

## Reset row names of data frame for further indexing
rownames(true.assign) <- NULL

## Remove underscore from species names
true.assign$Assignment <- gsub("_", " ", true.assign$Assignment)


#'
#' Now, correct species names:
#' 
#' - *Anas carolinensis* = *Anas*, row 5
#' - *Buteo* = *Buteo buteo*, row 18
#' - *Canis lupus* = *Canis lupus familiaris*, row 20
#' - *Castor* = *Castor fiber*, row 22
#' - *Emberiza chrysophrys* = *Emberiza*, row 36
#' - *Larus glaucoides* = *Larus*, row 49
#' - *Pelophylax* = *Pelophylax ridibundus*, row 70
#' - *Sus scrofa* = *Sus scrofa domesticus*, row 94
#' - *Triturus* = *Triturus cristatus*, row 96
#' - *Bison* = *Bison bonasus*, row 104
#' - *Bos* = *Bos taurus*, row 106
#' - Cichlidae = *Maylandia zebra*, row 110
#' - *Haplochromis burtoni* = *Maylandia zebra*, row 115
#' - *Lynx pardinus* = *Lynx lynx*, row 118
#' - *Meleagris* = *Meleagris gallopavo*, row 119
#' - *Oreochromis niloticus* = *Maylandia zebra*, row 121
#' - *Pundamilia nyererei* = *Maylandia zebra*, row 123
#' - *Sprattus* = *Sprattus sprattus*, row 126
#' - *Strix* = *Strix aluco*, row 127

true.assign$Assignment <- as.character(true.assign$Assignment)
true.assign[5, "Assignment"] <- "Anas"
true.assign[18, "Assignment"] <- "Buteo buteo"
true.assign[20, "Assignment"] <- "Canis lupus familiaris"
true.assign[22, "Assignment"] <- "Castor fiber"
true.assign[36, "Assignment"] <- "Emberiza"
true.assign[49, "Assignment"] <- "Larus"
true.assign[70, "Assignment"] <- "Pelophylax ridibundus"
true.assign[94, "Assignment"] <- "Sus scrofa domesticus"
true.assign[96, "Assignment"] <- "Triturus cristatus"
true.assign[104, "Assignment"] <- "Bison bonasus"
true.assign[106, "Assignment"] <- "Bos taurus"
true.assign[c(110,115,121,123), "Assignment"] <- "Maylandia zebra"
true.assign[118, "Assignment"] <- "Lynx lynx"
true.assign[119, "Assignment"] <- "Meleagris gallopavo"
true.assign[126, "Assignment"] <- "Sprattus sprattus"
true.assign[127, "Assignment"] <- "Strix aluco"

## Now merge read counts again by taxonomic assignment
## Anything with the same name in the data frame will be merged
true.assign <- ddply(true.assign, .(Assignment), numcolwise(sum))



#' ---
#' 
#' ## 3) Clean up dataset
#' 
#' Now the data is in a form that can be manipulated easily, it must be 
#' filtered to remove potential contaminants and false positives. 
#' 
#' There are several ways of doing this:
#' 1. Identify highest level of cichlid DNA contamination across all eDNA
#'    samples.
#' 3. Identify the highest level of contamination in positive (non-cichlid
#'    DNA) controls.
#' 4. Identify taxon-specific thresholds using positive controls, i.e. 
#'    the frequency required to remove a given taxa from the positive 
#'    control as only cichlid DNA should be present.
#' 
#' Arguably, taxon-specific thresholds based on positive controls are 
#' more effective as negative controls have no template DNA for 
#' contaminant DNA to compete with for amplification, thus contaminant
#' DNA amplifies exponentially. However, only 16 positive controls were
#' included on this MiSeq run which may render this approach ineffective.
#'

####################################################
# OPTION 1: highest level of cichlid contamination #
####################################################

## Create copy of dataframe without positive or negative controls for 
## threshold determination
cichlid <- true.assign[,which(!grepl("CB|FB|EB|Negative|Positive", 
                                     colnames(true.assign)))]

## Make Assignment column row names
rownames(cichlid) <- cichlid$Assignment
cichlid <- cichlid[,-1]

## Calculate total number of reads in samples
cichlid <- rbind(cichlid, colSums(cichlid))

## Create new dataframe containing the frequency of reads in each sample
cichlid.freq  <- cichlid/c(cichlid[115,])
cichlid.freq[is.na(cichlid.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency of
## DNA from each taxon across samples
cichlid.freq$Threshold <- apply(cichlid.freq, 1, max)

## Check Threshold column has been created properly
head(cichlid.freq[,220:221])

## Combine read frequencies with taxonomic assignment
species <- data.frame(true.assign$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
cichlid.freq <- cbind(species, cichlid.freq)

## Print contamination threshold based on max. level of cichlid
## DNA contamination
max(cichlid.freq[64,222])   # 0.3077%

## In this scenario, any assignments <0.3077% total reads in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
cichlid.test <- cichlid.freq
cichlid.test[cichlid.test <= 0.3077] <- 0

## Now convert back into read counts.
## Remove last row containing frequencies and add the total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
cichlid.test <- cichlid.test[-115,-222]
rownames(cichlid.test) <- NULL

total.counts <- data.frame(colSums(cichlid[-115,]))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
cichlid.conversion <- smartbind(cichlid.test, total.counts)
rownames(cichlid.conversion) <- cichlid.conversion$Assignment
cichlid.conversion <- cichlid.conversion[,-1]
cichlid.FP <- cichlid.conversion*c(cichlid.conversion[115,])

## Remove total row, reset row names, recreate Assignment column
cichlid.FP <- cichlid.FP[-115,]
cichlid.FP$Assignment <- rownames(cichlid.FP)
cichlid.FP <- cichlid.FP[,c(221,1:220)]
rownames(cichlid.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- cichlid[-115,]
temp$Assignment <- cichlid.FP$Assignment
temp <- temp[,c(221,1:220)]
temp[,2:221][temp[,2:221] > 0] <- 1
cichlid.FP[,2:221][cichlid.FP[,2:221] > 0] <- 1
cichlid1 <- data.frame(colSums(temp[,2:221]))
cichlid2 <- data.frame(colSums(cichlid.FP[,2:221]))
cichlid.compare <- cbind(cichlid1, cichlid2)

## Calculate proportion of information lost
cichlid.compare$proportion <- cichlid.compare[,2]/cichlid.compare[,1]*100
cichlid.compare[is.na(cichlid.compare)] <- 0

## Examine mean and range of proportion of taxa retained
range(cichlid.compare$proportion)
mean(cichlid.compare$proportion[cichlid.compare$proportion!=0])

## This would result in up to 100% taxa being removed.
## On average, 16.62% species detections are retained.


########################################################
# OPTION 2: highest level of contamination in controls #
########################################################

## Store mock community data, extraction blanks and controls in new 
## dataframes
neg.controls <- true.assign %>% select(Assignment, matches("CB|FB|EB|Negative"))
pos.controls <- true.assign %>% select(Assignment, contains("Positive"))

## Positive and negative controls were subset into seperate dataframes
## earlier in script
head(pos.controls)
head(neg.controls)

## Check positive controls for highest level of contamination (any DNA
## except cichlid). Remove column with taxonomic assignment.
pos <- pos.controls[,-1]

## Calculate total number of reads in samples
pos <- rbind(pos, colSums(pos))

## Create new dataframe containing the frequency of reads in each sample
pos.freq <- pos/c(pos[115,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
pos.freq$Threshold <- apply(pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(pos.freq)

## Combine read frequencies with taxonomic assignment
species <- data.frame(pos.controls$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
pos.freq <- cbind(species, pos.freq)

## Print contamination threshold based on max. level of non-cichlid
## DNA contamination
## Exclude unassigned from threshold determination due to large number of
## reads, which cannot be distinguished as poorly amplified cichlid DNA 
## or other taxa
rownames(pos.freq) <- pos.freq$Assignment
pos.freq <- pos.freq[,-1]
max(pos.freq[-c(64,113,115),17])   # 0.0644%

## In this scenario, any assignments <0.0644% total reads in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
pos.test <- cichlid.freq[,-222]
pos.test[pos.test <= 0.0644] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
pos.test <- pos.test[-115,]
rownames(pos.test) <- NULL

total.counts <- data.frame(colSums(cichlid[-115,]))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
pos.conversion <- smartbind(pos.test, total.counts)
rownames(pos.conversion) <- pos.conversion$Assignment
pos.conversion <- pos.conversion[,-1]
pos.FP <- pos.conversion*c(pos.conversion[115,])

## Remove total row, reset row names and recreate Assignment column
pos.FP <- pos.FP[-115,]
pos.FP$Assignment <- rownames(pos.FP)
pos.FP <- pos.FP[,c(221,1:220)]
rownames(pos.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- cichlid[-115,]
temp$Assignment <- pos.FP$Assignment
temp <- temp[,c(221,1:220)]
temp[,2:221][temp[,2:221] > 0] <- 1
pos.FP[,2:221][pos.FP[,2:221] > 0] <- 1
pos1 <- data.frame(colSums(temp[,2:221]))
pos2 <- data.frame(colSums(pos.FP[,2:221]))
pos.compare <- cbind(pos1, pos2)

## Calculate proportion of information lost
pos.compare$proportion <- pos.compare[,2]/pos.compare[,1]*100
pos.compare[is.na(pos.compare)] <- 0

## Examine mean and range of proportion of taxa retained
range(pos.compare$proportion)
mean(pos.compare$proportion)

## This would result in up to 100% taxa being removed.
## On average, 30.38% species detections are retained.


#######################################
# OPTION 3: Taxon-specific thresholds #
#######################################

## Check positive controls for highest level of contamination (any DNA
## except cichlid). Remove column with taxonomic assignment.
pos <- pos.controls[,-1]

## Calculate total number of reads in samples
pos <- rbind(pos, colSums(pos))

## Create new dataframe containing the frequency of reads in each sample
pos.freq <- pos/c(pos[115,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
pos.freq$Threshold <- apply(pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(pos.freq)

## Combine read frequencies with taxonomic assignment
species <- data.frame(pos.controls$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
pos.freq <- cbind(species, pos.freq)

## Manually change the species threshold value for cichlid to 0 as 
## this will be a true contaminant in the dataset and no false positive 
## threshold required.
pos.freq$Threshold[64] <- 0

## Check this manual edit has occurred
head(pos.freq[64,18])

## In this scenario, any assignments less than threshold in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Apply thresholds
SS.test <- cichlid.freq[,-222]
SS.test[SS.test <= pos.freq$Threshold] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
SS.test <- SS.test[-115,]
rownames(SS.test) <- NULL

total.counts <- data.frame(colSums(cichlid[-115,]))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
SS.conversion <- smartbind(SS.test, total.counts)
rownames(SS.conversion) <- SS.conversion$Assignment
SS.conversion <- SS.conversion[,-1]
SS.FP <- SS.conversion*c(SS.conversion[115,])

## Remove total row, reset row names and recreate Assignment column
SS.FP <- SS.FP[-115,]
SS.FP$Assignment <- rownames(SS.FP)
SS.FP <- SS.FP[,c(221,1:220)]
rownames(SS.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- cichlid[-115,]
temp$Assignment <- SS.FP$Assignment
temp <- temp[,c(221,1:220)]
temp[,2:221][temp[,2:221] > 0] <- 1
SS.FP[,2:221][SS.FP[,2:221] > 0] <- 1
SS1 <- data.frame(colSums(temp[,2:221]))
SS2 <- data.frame(colSums(SS.FP[,2:221]))
SS.compare <- cbind(SS1,SS2)

## Calculate proportion of information lost
SS.compare$proportion <- SS.compare[,2]/SS.compare[,1]*100
SS.compare[is.na(SS.compare)] <- 0

## Examine mean and range of proportion of taxa retained
range(SS.compare$proportion)
mean(SS.compare$proportion)

## This would result in up to 43% taxa being removed.
## On average, 93.72% species detections are retained.


###########
# SUMMARY # 
###########

## Tidy dataframes
cichlid.compare$Sample <- rownames(cichlid.compare)
cichlid.compare <- cichlid.compare[,c(4,1:3)]
rownames(cichlid.compare) <- NULL
colnames(cichlid.compare)[2:4] <- c("sp_richness_NT",
                                    "sp_richness_TA",
                                    "prop_TA")
cichlid.compare$type <- "Cichlid"

pos.compare$Sample <- rownames(pos.compare)
pos.compare <- pos.compare[,c(4,1:3)]
rownames(pos.compare) <- NULL
colnames(pos.compare)[2:4] <- c("sp_richness_NT",
                                "sp_richness_TA",
                                "prop_TA")
pos.compare$type <- "Positive"

SS.compare$Sample <- rownames(SS.compare)
SS.compare <- SS.compare[,c(4,1:3)]
rownames(SS.compare) <- NULL
colnames(SS.compare)[2:4] <- c("sp_richness_NT",
                               "sp_richness_TA",
                               "prop_TA")
SS.compare$type <- "Taxon-specific"

## Combine dataframes
threshold.df <- rbind(cichlid.compare, pos.compare, SS.compare)

## Plot number of species retained after thresholds applied
p1 <- ggplot(dat=threshold.df, aes(x=Sample, y=sp_richness_TA, fill=type))
p1 <- p1 + geom_bar(stat="identity", position=position_dodge())
p1 <- p1 + scale_y_continuous(limits=c(0,30))
p1 <- p1 + labs(x="Sample", y="Detections remaining after threshold application")
p1 <- p1 + theme_bw()
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.y = element_text(colour="black"),
                 legend.position = "none",
                 text = element_text(size=20))
p1 <- p1 + facet_grid(type ~ .)
p1

## Plot proportion of species retained after thresholds applied
p1 <- ggplot(dat=threshold.df, aes(x=Sample, y=prop_TA, fill=type))
p1 <- p1 + geom_bar(stat="identity", position=position_dodge())
p1 <- p1 + scale_y_continuous(limits=c(0,100))
p1 <- p1 + labs(x="Sample", 
                y="Detections remaining after threshold application (%)")
p1 <- p1 + theme_bw()
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 #axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.y = element_text(colour="black"),
                 legend.position = "none",
                 text = element_text(size=20))
p1 <- p1 + facet_grid(type ~ .)
p1

## Going forward the taxon-specific thresholds will be used as these
## retain the majority of biological information.
## Recreate dataframe with threshold applied.
SS.FP <- SS.conversion*c(SS.conversion[115,])
SS.FP <- SS.FP[-115,]
SS.FP$Assignment <- rownames(SS.FP)
SS.FP <- SS.FP[,c(221,1:220)]
rownames(SS.FP) <- NULL



#################
# CONTAMINATION #
#################

## Examine how much contamination occured in negative controls (field blanks,
## filtration blanks, extraction blanks, and PCR negative controls)
contamination <- neg.controls

## Make Assignment column row names
rownames(contamination) <- contamination$Assignment
contamination <- contamination[,-1]

## Calculate total number of reads in each control
contamination <- rbind(contamination, colSums(contamination))

## Calculate frequency of reads in each control
contamination <- contamination/c(contamination[115,])

## Replace NA values with 0
contamination[is.na(contamination)] <- 0

## Remove total frequency row from dataframe
contamination <- contamination[-115,]

## Make row names first column in dataframe
contamination$Assignment <- rownames(contamination)
contamination <- contamination[,c(90,1:89)]
rownames(contamination) <- NULL

## Transpose dataframe for negative controls
contaminants <- setNames(data.frame(t(contamination[,-1])), contamination[,1])

## Make row names first column in dataframe
contaminants$ID <- rownames(contaminants)
contaminants <- contaminants[,c(115,1:114)]
rownames(contaminants) <- NULL

## Create column specifying type of negative control
contaminants$Type <- ifelse(grepl("CB", contaminants$ID), "Field",
                            ifelse(grepl("FB", contaminants$ID), "Filtration",
                                   ifelse(grepl("EB", contaminants$ID), "Extraction",
                                          ifelse(grepl("Negative", contaminants$ID), "PCR",
                                                 "Other"))))

## Move to start of dataframe
contaminants <- contaminants[,c(1,116,2:115)]

## Melt dataframe for plotting
contaminants <- melt(contaminants, id=c("ID","Type"))

## Rename columns
colnames(contaminants)[3:4] <- c("Assignment", "Frequency")

## Create factor to order heatmap by
contaminants$ftype <- factor(contaminants$Type, 
                             levels=c("Field","Filtration","Extraction","PCR"))

## Plot contamination found in negative controls
hm0 <- ggplot(contaminants, aes(x=ID, 
                                y=fct_rev(as_factor(Assignment)), 
                                fill=Frequency))
hm0 <- hm0 + geom_tile(colour="grey80")
hm0 <- hm0 + scale_fill_gradientn(name="Proportional\nread counts", 
                                  limits=c(0,1),
                                  breaks=c(0,0.25,0.50,0.75,1),
                                  labels=scales::number_format(accuracy=0.01,
                                                               decimal.mark="."),
                                  colours=c("white","red","black"), 
                                  values=c(0,0.1,1))
hm0 <- hm0 + labs(x="Process controls", y="Taxonomic assignment")
hm0 <- hm0 + theme_bw()
hm0 <- hm0 + theme(panel.grid.major = element_line(colour="white"),
                   panel.grid.minor = element_line(colour="white"), 
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(colour="black"),
                   text = element_text(size=20),
                   legend.key.size = unit(1, 'lines'))
hm0 <- hm0 + facet_grid(. ~ ftype, scale="free_x",space="free_x")
hm0



#' ---
#'
#' ## 4) Spurious assigments
#'
#' Remove any non-vertebrate assignments and vertebrate assignments higher 
#' than species level as these are too coarse. Also remove positive control
#' and domestic species:
#' 

## Remove any assignments above species level, excluding the genera Anas
## Emberiza, and Larus
spp.df <- SS.FP[which(grepl(" |Anas|Emberiza|Larus", SS.FP$Assignment)),]

## Reset row names of data frame for further indexing
rownames(spp.df) <- NULL

## Remove positive control, true contaminants (i.e. penguin, reindeer), 
## prey items given to captive animals (i.e. chicks)
spp.df <- spp.df[-c(24,27,34,51,53),]


##########################
# EXPERIMENT 1 DATAFRAME #
##########################

## Create dataframe for Experiment 1: 
## eDNA detection and signal strength in artificial systems
cap.df <- spp.df[,which(grepl("Assignment|HWP|WWT", colnames(spp.df)))]
cap.df <- cap.df[which(grepl("Lutra|Arvicola|Castor|Lynx|Erinaceus|Meles|Cervus|Sciurus vulgaris|Martes", 
                             cap.df$Assignment)),]
rownames(cap.df) <- NULL

## Extract samples belonging to Experiment 1 from raw totals dataframe
exp1.total <- raw.total[,which(grepl("Assignment|HWP|WWT", colnames(raw.total)))]

## Bind to previously created dataframe
cap.df <- rbind(cap.df, exp1.total)

## Make Assignment column row names
rownames(cap.df) <- cap.df$Assignment
cap.df <- cap.df[,-1]

## Calculate proportional read counts for each mammal species and normalise by
## multiplying by the average read count per sample
cap.df <- cap.df/c(cap.df[10,])

## Remove total row
cap.df <- cap.df[-10,]

## Make row names first column of dataframe
cap.df$Assignment <- rownames(cap.df)
cap.df <- cap.df[,c(81,1:80)]
rownames(cap.df) <- NULL


##########################
# EXPERIMENT 2 DATAFRAME #
##########################

## Create dataframe for Experiment 2:
## eDNA detection and signal strength in natural systems
wild.df <- spp.df[,which(!grepl("HWP|WWT", colnames(spp.df)))]
rownames(wild.df) <- wild.df$Assignment
wild.df <- wild.df[,-1]

## Extract samples belonging to Experiment 2 from raw totals dataframe
exp2.total <- raw.total[,which(!grepl("Assignment|HWP|WWT", colnames(raw.total)))]

## Bind to previously created dataframe
wild.df <- rbind(wild.df, exp2.total)

## Calculate proportional read counts for each mammal species and normalise by
## multiplying by the average read count per sample
wild.df <- wild.df/c(wild.df[63,])

## Remove total row
wild.df <- wild.df[-63,]

## Make row names first column of dataframe
wild.df$Assignment <- rownames(wild.df)
wild.df <- wild.df[,c(141,1:140)]
rownames(wild.df) <- NULL



#' ---
#' 
#' ## 5) Experiment 1: data analysis
#' 
#' The primary research questions for this experiment are:
#' 
#' - Is eDNA detection with the passive sampling strategy comparable to 
#'   the targeted sampling strategy based on behavioural observation? 
#'   More simply, is eDNA from each mammal species locally or evenly 
#'   distributed in pond water? 
#' - Is eDNA distribution and strength related to type of mammal species 
#'   (semi-aquatic, ground-dwelling, arboreal)?
#' - Do read counts for target mammal spp correlate with mammal type 
#'   and density? 
#' - Do read counts correlate with behaviour frequency/duration?
#' 

## Import metadata for each sample (i.e. volume filtered, filters used)
cap.sample.metadata <- read.csv("../Data/Captive_sample_metadata.csv", header=TRUE)
head(cap.sample.metadata)
str(cap.sample.metadata)
names(cap.sample.metadata)
dim(cap.sample.metadata)

## Remove column containing notes
cap.sample.metadata <- cap.sample.metadata[,-6]

## Import metadata for each waterbody in enclosures
cap.pond.metadata <- read.csv("../Data/Captive_pond_metadata.csv", header=TRUE)
head(cap.pond.metadata)
str(cap.pond.metadata)
names(cap.pond.metadata)
dim(cap.pond.metadata)

## Transpose dataframe for captive mammal samples
cap.sample.df <- setNames(data.frame(t(cap.df[,-1])), cap.df[,1])

## Make row names a column in data frame
cap.sample.df$Sample <- rownames(cap.sample.df)

## Move to start of dataframe
cap.sample.df <- cap.sample.df[,c(10,1:9)]

## Reset row names
rownames(cap.sample.df) <- NULL

## Bind dataframe with metadata for samples
cap.sample.df <- merge(cap.sample.df, cap.sample.metadata, by="Sample")

## Move metadata to beginning of dataframe
cap.sample.df <- cap.sample.df[,c(1,11:14,2:10)]

## Set non-target species detections to 0 in each sample
cap.sample.df[c(1:4,45:46,75:76),c(6:10,12:14)][cap.sample.df[c(1:4,45:46,75:76),c(6:10,12:14)] > 0] <- 0
cap.sample.df[c(6:12,47:50,69:74),c(6,8:14)][cap.sample.df[c(6:12,47:50,69:74),c(6,8:14)] > 0] <- 0
cap.sample.df[c(13:22,51:56,77:80),c(6:7,9:14)][cap.sample.df[c(13:22,51:56,77:80),c(6:7,9:14)] > 0] <- 0
cap.sample.df[c(23:27,57:59),c(6:11,13:14)][cap.sample.df[c(23:27,57:59),c(6:11,13:14)] > 0] <- 0
cap.sample.df[28:31,6:13][cap.sample.df[28:31,6:13] > 0] <- 0
cap.sample.df[32:33,7:14][cap.sample.df[32:33,7:14] > 0] <- 0
cap.sample.df[34:35,c(6:8,10:14)][cap.sample.df[34:35,c(6:8,10:14)] > 0] <- 0
cap.sample.df[c(36:40,60:65),c(6:9,11:14)][cap.sample.df[c(36:40,60:65),c(6:9,11:14)] > 0] <- 0
cap.sample.df[c(41:44,66:68),c(6:12,14)][cap.sample.df[c(41:44,66:68),c(6:12,14)] > 0] <- 0

## Melt dataframe
cap.sample.df <- melt(cap.sample.df, id=c("Sample","Location","Individuals",
                                          "Volume","No_filters"))

## Rename columns 
colnames(cap.sample.df) <- c("ID","Location","No_ind","Volume",
                             "No_filters","Species","Prop_reads")

## Replace POOL for Lynx at Wildwood Trust with T01 in ID and sample number columns
cap.sample.df$ID <- gsub("POOL", "T01", cap.sample.df$ID)

## Create column specifying whether samples are passive, targeted or other
## (i.e. designated drinking source, quarantine enclosure)
cap.sample.df$Sample_type <- ifelse(grepl("P0",cap.sample.df$ID), "Stratified",
                                    ifelse(grepl("T0",cap.sample.df$ID), "Directed", "Other"))

## Duplicate Sample column
cap.sample.df$Enclosure <- cap.sample.df$ID

## Remove abbreviated location, filtration day and round from sample ID
cap.sample.df$Enclosure  <- gsub("HWP|WWT", "", cap.sample.df$Enclosure)
cap.sample.df$Enclosure  <- gsub("-D1|-D2|-D3|-D01|-D02|-D03", "", cap.sample.df$Enclosure)
cap.sample.df$Enclosure  <- gsub("-Round01-|-Round02-|-Round03-|-Round04-", "", cap.sample.df$Enclosure)

## Split Enclosure column into multiple columns
cap.sample.df <- cap.sample.df %>% separate(Enclosure, into=paste("id", 1:4, sep = ""))

## Remove redundant columns
cap.sample.df <- cap.sample.df[,-c(11:12)]

## Rename new columns
names(cap.sample.df)[9:10] <- c("Enclosure","Sample_no")

## Rename sample numbers to reflect directed or stratified sampling strategy
cap.sample.df$Sample_no <- gsub("P01","STR01",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("P02","STR02",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("P03","STR03",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("P04","STR04",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("P05","STR05",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("P06","STR06",cap.sample.df$Sample_no)

cap.sample.df$Sample_no <- gsub("T01","DIR01",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("T02","DIR02",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("T03","DIR03",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("T04","DIR04",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("T05","DIR05",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("T06","DIR06",cap.sample.df$Sample_no)

## Change names of water vole samples
cap.sample.df$Sample_no <- gsub("1IND","QUAR1",cap.sample.df$Sample_no)
cap.sample.df$Sample_no <- gsub("4IND","QUAR2",cap.sample.df$Sample_no)

## Swap beaver and zoo in enclosure and sample number columns
cap.sample.df$Enclosure <- gsub("ZOO", "BEAV", cap.sample.df$Enclosure)
cap.sample.df$Sample_no <- gsub("BEAV", "ZOO", cap.sample.df$Sample_no)

## Create column specifying the type of mammal species
cap.sample.df$Species_type <- ifelse(cap.sample.df$Enclosure %in% c("OTT","WV","BEAV"), "Semi-aquatic",
                                     ifelse(cap.sample.df$Enclosure %in% c("BAD","DEER","HH","LYNX"), "Ground-dwelling",
                                            ifelse(cap.sample.df$Enclosure %in% c("PM","SQ"), "Arboreal", "Other")))

## Reorder dataframe
cap.sample.df <- cap.sample.df[,c(1:2,9,3,11,10,8,4:7)]

## Now, we need to remove all species assignments that shouldn't be in a 
## given enclosure. Create dataframe for each species and then concatenate.
otter <- subset(cap.sample.df , Enclosure == "OTT" & Species == "Lutra lutra")
wv <- subset(cap.sample.df , Enclosure == "WV" & Species == "Arvicola amphibius")
beav <- subset(cap.sample.df , Enclosure == "BEAV" & Species == "Castor fiber")
hh <- subset(cap.sample.df , Enclosure == "HH" & Species == "Erinaceus europaeus")
bad <- subset(cap.sample.df , Enclosure == "BAD" & Species == "Meles meles")
deer <- subset(cap.sample.df , Enclosure == "DEER" & Species == "Cervus elaphus")
lynx <- subset(cap.sample.df , Enclosure == "LYNX" & Species == "Lynx lynx")
pm <- subset(cap.sample.df , Enclosure == "PM" & Species == "Martes martes")
sq <- subset(cap.sample.df , Enclosure == "SQ" & Species %in% c("Sciurus vulgaris","Meles meles"))
sq <- sq[-c(2:5),]
sq[1,3] <- "BAD"
sq[1,5] <- "Ground-dwelling"

## Concatenate dataframes
cap.spp <- rbind(otter,wv,beav,hh,bad,deer,lynx,sq,pm)

## Check structure of dataframe before proceeding with downstream analyses
str(cap.spp)

## Make several variables in dataframe factors instead of character/integer
cap.spp$ID <- as.factor(cap.spp$ID)
cap.spp$Enclosure <- as.factor(cap.spp$Enclosure)
cap.spp$Species_type <- as.factor(cap.spp$Species_type)
cap.spp$Sample_no <- as.factor(cap.spp$Sample_no)
cap.spp$Sample_type <- as.factor(cap.spp$Sample_type)
cap.spp$Species <- as.factor(cap.spp$Species)

## Finally, reset row names in case of further indexing
rownames(cap.spp) <- NULL


############
# RAW DATA #
############

## Only hedgehog (drinking bowl) and red deer (random samples) were not detected
## in every sample collected from species enclosures.
## Make ordered sample type variable in dataframe for plotting
cap.spp$fSample_type <- factor(cap.spp$Sample_type, levels=c("Directed","Stratified","Other"))

## Disable scientific notation
options(scipen = 999)

## Plot with common species names and tiles coloured by species lifestyle
cap.spp$fSpecies <- factor(cap.spp$Species, 
                           levels=c("Arvicola amphibius",
                                    "Castor fiber",
                                    "Lutra lutra",
                                    "Erinaceus europaeus",
                                    "Meles meles",
                                    "Lynx lynx",
                                    "Cervus elaphus",
                                    "Sciurus vulgaris",
                                    "Martes martes"))

cap.spp$fLifestyle <- factor(cap.spp$Species_type,
                             levels=c("Semi-aquatic",
                                      "Ground-dwelling",
                                      "Arboreal"))

hm1 <- ggplot(cap.spp, aes(x=Sample_no, y=fct_rev(as_factor(fSpecies))))
hm1 <- hm1 + geom_tile(aes(fill=Prop_reads, colour=fLifestyle), size=1)
hm1 <- hm1 + geom_text(aes(label=round(Prop_reads, 4)), cex=2.5)
hm1 <- hm1 + scale_colour_manual(name="Species lifestyle",
                                 values=c("dodgerblue","limegreen","goldenrod1"))
hm1 <- hm1 + scale_fill_gradient(name="Proportional\nread counts", 
                                 limits=c(0,1), 
                                 breaks=c(0,0.25,0.5,0.75,1),
                                 low="white", high="gray50",
                                 guide=guide_colourbar(frame.colour="black",
                                                       ticks.colour="black"))
hm1 <- hm1 + labs(x="Sample", y="Species")
hm1 <- hm1 + theme_bw()
hm1 <- hm1 + theme(panel.grid.major = element_line(colour="white"),
                   panel.grid.minor = element_line(colour="white"), 
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black", angle = 45, vjust=1, hjust=1),
                   axis.text.y = element_text(face="italic", colour="black"),
                   text = element_text(size=20),
                   legend.key.size = unit(1, 'lines'))
hm1 <- hm1 + facet_grid(Location ~ fSample_type, scale="free_x", space="free_x")
hm1


## Data exploration
summary(cap.spp)         # summary of all variables in dataframe
max(cap.spp$Prop_reads)  # max. proportional read count
min(cap.spp$Prop_reads)  # min. proportional read count
which(cap.spp$Prop_reads == (max(cap.spp$Prop_reads)))  # observation with highest proportional read count (pine marten)
which(cap.spp$Prop_reads == (min(cap.spp$Prop_reads)))  # observation with lowest proportional read count (hedgehog, deer)
range(cap.spp$Prop_reads)     # range of proportional read counts
quantile(cap.spp$Prop_reads)  # quantiles
IQR(cap.spp$Prop_reads)       # Interquartile range
quantile(cap.spp$Prop_reads, c(0.05, 0.95))  # 5th and 95th percentiles
var(cap.spp$Prop_reads, na.rm=TRUE)  # variance
sd(cap.spp$Prop_reads)               # standard deviation

## Exploratory plots
boxplot(cap.spp$Prop_reads)           
plot(cap.spp$Prop_reads ~ cap.spp$ID)           
plot(cap.spp$Prop_reads ~ cap.spp$Location)     # Model as random
plot(cap.spp$Prop_reads ~ cap.spp$Enclosure)    # Model as fixed
plot(cap.spp$Prop_reads ~ cap.spp$No_ind)       # Omit, not enough data points
plot(cap.spp$Prop_reads ~ cap.spp$Species_type) # Model as fixed
plot(cap.spp$Prop_reads ~ cap.spp$Sample_no)    # Only used for summarising detection
plot(cap.spp$Prop_reads ~ cap.spp$Sample_type)  # Model as fixed
plot(cap.spp$Prop_reads ~ cap.spp$Volume)       # Model as fixed
plot(cap.spp$Prop_reads ~ cap.spp$No_filters)   # Model as fixed

## Check distribution of data using histogram
hist(cap.spp$Prop_reads)  # somewhat skewed to right

## Plot quantile-quantile plot to compare observations to a normal distribution
qqnorm(cap.spp$Prop_reads)
qqline(cap.spp$Prop_reads)  # data tails away from normal distribution

## Run quantitative tests for the normal distribution
shapiro.test(cap.spp$Prop_reads)   # data not normal, P = 0.00005276
ks.test(cap.spp$Prop_reads, pnorm) # data not normal, P < 2.2e-16

## The proportional data is not normally distributed. However, this is not
## unexpected. Strictly bounded data (i.e. proportions) are continuous
## but often strongly non-normal. If the data cannot be transformed to make
## it more normal (as is the case with our data), then Generalised Linear
## Models are most appropriate for performing regression analyses. These
## are equipped to deal with dependent variables where the variance is not 
## uniform across its range, the residuals from a General Linear Model are 
## not normally distributed, and there may be non-linear relationships 
## between dependent and independent variables. Proportional data falls
## under a binomial distribution, as the data is bounded by 0 and 1.
## As we also want to account for the influence of sampling location and
## occasionally species but not model the effect of these variables directly, 
## we need to use Generalised Linear Mixed Models (GLMMs) to model these 
## variables as random effects. Therefore, we will use binomial GLMMs with 
## the logit link function to model the proportional data.


#################################################
# DIFFERENCE BETWEEN SAMPLE TYPE ACROSS SPECIES #
#################################################

## Make vectors for stratified and directed samples
stratified <- cap.spp[which(!grepl("Stratified", cap.spp$Sample_type)),]
stratified <- stratified$Prop_reads
directed <- cap.spp[which(!grepl("Directed", cap.spp$Sample_type)),]
directed <- directed$Prop_reads

## Make dataframe containing only directed and stratified samples
sampling.df <- cap.spp[which(!grepl("Other", cap.spp$Sample_type)),]
sampling.df$Sample_type <- droplevels(sampling.df$Sample_type)

## Test for difference between mean read counts of targeted and passive samples 
## across all species.
## Check variance around each group:
leveneTest(sampling.df$Prop_reads, sampling.df$Sample_type, location = "mean")  

## Medians for each group
median(sampling.df$Prop_reads[sampling.df$Sample_type == "Stratified"])
median(sampling.df$Prop_reads[sampling.df$Sample_type == "Directed"])

## Mann-Whitney test (used for non-normal distributions)
wilcox.test(stratified, directed, exact=TRUE, conf.int=TRUE)

## Plot results
p1 <- ggplot(sampling.df, aes(x=Sample_type, y=Prop_reads))
p1 <- p1 + geom_jitter(aes(colour=Sample_type), pch=16, cex=2, width=0.2)
p1 <- p1 + geom_boxplot(alpha=0.7, outlier.colour="black", outlier.size=2)
p1 <- p1 + scale_colour_manual(values=c("dodgerblue","goldenrod1"))
p1 <- p1 + labs(x="Sample type", y="Proportional read counts")
p1 <- p1 + scale_y_continuous(limits=c(0,1))
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0, colour="black"),
                 plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
                 text = element_text(size=20),
                 legend.position="none",
                 legend.key=element_blank(),
                 legend.key.size = unit(2, 'lines'))
p1


##############################################
# DIFFERENCES BETWEEN SAMPLE TYPE BY SPECIES #
##############################################

## First, we need to remove species that don't have observations of both
## sample types. Only targeted samples were taken for water vole, hedgehog, 
## and red squirrel thus these must be removed
m1.df <- sampling.df[which(!grepl("WV|HH|SQ", sampling.df$Enclosure)),]
m1.df$Enclosure <- droplevels(m1.df$Enclosure)
m1.df$Species <- droplevels(m1.df$Species)

## Fit Poisson GLMM to examine the difference between targeted and passive 
## samples for each species. Use GLMM to account for the spatial variation 
## in samples by modeling wildlife park as a random factor.
m1 <- glmmTMB(Prop_reads ~ (1|Location) + Species/Sample_type,
              data=m1.df,
              family=binomial(link="logit"))
summary(m1)
drop1(m1, test = "Chi")

## Calculate R-squared of model
r.squaredGLMM(m1)
## marginal R-squared = 39.21% (proportion of variance explained by fixed effects)
## conditional R-squared = 39.21% (proportion of variance explained by fixed + random effects)

## Examine difference between species and between sample types
summary(glht(m1, linfct=mcp(Species="Tukey")))
summary(glht(m1, linfct=mcp(Sample_type="Tukey")))

## Test for overdispersion
library(sjstats)
overdisp(m1)

## Plot the fitted data against the observed data
plot(m1.df$Prop_reads ~ fitted(m1))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(m1.df$Prop_reads, fitted(m1))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(m1, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.3713

## No deviation from normality as residuals normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ m1.df$Sample_type)   
plot(sresid ~ m1.df$Species)   

## Assumption 3: no collinearity
## All variables are factors.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## PLOT MODEL FIT
## Make a table of prediction data, containing the values of fixed factors
## for which to make predictions of the dependent variable
pdat <- expand.grid(Location = levels(m1.df$Location),
                    Species = levels(m1.df$Species),
                    Sample_type = levels(m1.df$Sample_type))
pdat

## Make a dataframe containing the predicted data
## re.form=NA means predictions are made for random effect as a whole
## se.fit=TRUE obtains the standard error for each of the predictions
predictions <- data.frame(predict(m1, newdata=pdat, 
                                  se.fit=TRUE, 
                                  na.action=na.exclude, 
                                  type="response"))

## Combine predictions with the prediction data into a dataframe for plotting
preds.for.plot <- data.frame(pdat, predictions)

## Remove predicted data for Highland Wildlife Park
preds.for.plot <- preds.for.plot[which(preds.for.plot$Location=="Wildwood Trust"),]

## Calculate and add upper/lower 95% CIs from the SEs of the predictions
preds.for.plot$ciu <- (preds.for.plot$fit+1.96*preds.for.plot$se.fit) 
preds.for.plot$cil <- (preds.for.plot$fit-1.96*preds.for.plot$se.fit)

## Check dataframe for plotting
head(preds.for.plot) 

## Plot predicted data against observed data
m1.df$fSpecies <- factor(m1.df$fSpecies)

p2 <- ggplot()
p2 <- p2 + geom_point(data=m1.df,
                      aes(x=fSpecies, y=Prop_reads, 
                          colour=fLifestyle,
                          shape=Sample_type,
                          group=Sample_type), 
                      position=position_jitterdodge(dodge.width=0.7), 
                      cex=5)
p2 <- p2 + geom_point(aes(x=preds.for.plot$Species, 
                          y=preds.for.plot$fit, 
                          shape=preds.for.plot$Sample_type), 
                      position=position_dodge(width=0.5), cex=5)
p2 <- p2 + geom_errorbar(aes(x=preds.for.plot$Species, 
                             group=preds.for.plot$Sample_type,
                             ymin=preds.for.plot$cil, ymax=preds.for.plot$ciu),
                         colour="black", size=0.8, width=0.3,
                         position=position_dodge(width=0.5))
p2 <- p2 + labs(title="(a)",
                x="Species", 
                y="Proportional read counts")
p2 <- p2 + scale_y_continuous(limits=c(0,1),
                              breaks=c(0,0.25,0.50,0.75,1.0),
                              oob=squish)
p2 <- p2 + scale_shape_manual(name="Sample type",
                              values=c(16,21))
p2 <- p2 + scale_colour_manual(name="Species lifestyle",
                               values=c("dodgerblue","limegreen","goldenrod2"))
p2 <- p2 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(face="italic", colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0, colour="black"),
                 plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
                 text = element_text(size=20),
                 legend.position="bottom",
                 legend.direction="horizontal",
                 legend.box="vertical",
                 legend.key=element_blank(),
                 legend.key.size = unit(2, 'lines'))
p2


#######################################
# 3. DIFFERENCE BETWEEN SPECIES TYPES #
#######################################

## Replicate dataframe used for first model
m2.df <- m1.df

## Examine influence of species type on read counts
m2 <- glmmTMB(Prop_reads ~ (1|Location/Species) + Species_type,
              data=m2.df,
              family=binomial(link="logit"))
summary(m2)
drop1(m2, test = "Chi")

## Calculate R-squared of model
r.squaredGLMM(m2)
## marginal R-squared = 6.86% (proportion of variance explained by fixed effects)
## conditional R-squared = 18.07% (proportion of variance explained by fixed + random effects)

## Examine difference between species and between sample types
summary(glht(m2, linfct=mcp(Species_type="Tukey")))

## Test for overdispersion
overdisp(m2)

## Plot the fitted data against the observed data
plot(m2.df$Prop_reads ~ fitted(m2))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(m2.df$Prop_reads, fitted(m2))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(m2, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.006619

## Some deviation from normality as residuals are not normally distributed
## therefore model may not be that reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ m2.df$Species_type)   

## Assumption 3: no collinearity
## All variables are factors.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag


## PLOT MODEL FIT
## Obtain predicted values for full data set: cap.spp
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
pdat2 <- expand.grid(Location = levels(m2.df$Location),
                     Species = levels(m2.df$Species),
                     Species_type = levels(m2.df$Species_type))
pdat2

## Make a dataframe containing the predicted data
## re.form=NA means predictions are made for random effect as a whole
## se.fit=TRUE obtains the standard error for each of the predictions
predictions2 <- data.frame(predict(m2, newdata=pdat2, se.fit=TRUE, 
                                   na.action=na.exclude, type="response"))

## Combine predictions with the prediction data into a dataframe for plotting
preds.for.plot2 <- data.frame(pdat2, predictions2)

## Keep only predicted data for one species and one location
preds.for.plot2 <- preds.for.plot2[which(preds.for.plot2$Location=="Wildwood Trust"),]
preds.for.plot2 <- preds.for.plot2[which(preds.for.plot2$Species=="Castor fiber"),]

## Calculate and add upper/lower 95% CIs from the SEs of the predictions
preds.for.plot2$ciu <- (preds.for.plot2$fit+1.96*preds.for.plot2$se.fit) 
preds.for.plot2$cil <- (preds.for.plot2$fit-1.96*preds.for.plot2$se.fit)

## Check dataframe for plotting
head(preds.for.plot2) 

## Plot predicted data against observed data
p3 <- ggplot()
p3 <- p3 + geom_point(data=m2.df,
                      aes(x=fLifestyle, y=Prop_reads, 
                          colour=fSpecies, 
                          shape=Location,
                          group=Location), 
                      position=position_jitterdodge(dodge.width=0.3),
                      cex=5)
p3 <- p3 + geom_point(aes(x=preds.for.plot2$Species_type, 
                          y=preds.for.plot2$fit), 
                      position=position_dodge(width=0.5), cex=5)
p3 <- p3 + geom_errorbar(aes(x=preds.for.plot2$Species_type, 
                             ymin=preds.for.plot2$cil, 
                             ymax=preds.for.plot2$ciu),
                         colour="black", size=0.8, width=0.3)
p3 <- p3 + labs(title="(b)",
                x="Species lifestyle", 
                y="Proportional read counts")
p3 <- p3 + scale_y_continuous(limits=c(0,1),
                              breaks=c(0,0.25,0.50,0.75,1.0),
                              oob=squish)
p3 <- p3 + scale_colour_manual(guide=guide_legend(nrow=1),
                               name="Species",
                               breaks=c("Castor fiber","Lutra lutra","Meles meles",
                                        "Lynx lynx","Cervus elaphus","Martes martes"),
                               labels=c(expression(italic("Castor fiber")),
                                        expression(italic("Lutra lutra")),
                                        expression(italic("Meles meles")),
                                        expression(italic("Lynx lynx")),
                                        expression(italic("Cervus elaphus")),
                                        expression(italic("Martes martes"))),
                               values=c("dodgerblue","dodgerblue4","greenyellow",
                                        "limegreen","green4","goldenrod2"))
p3 <- p3 + scale_shape_manual(name="Location",
                              values=c(15,17))
p3 <- p3 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0, colour="black"),
                 plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
                 text = element_text(size=20),
                 legend.position="bottom",
                 legend.direction="horizontal",
                 legend.box="vertical",
                 legend.key=element_blank(),
                 legend.key.size = unit(2, 'lines'))
p3

## Plot together
ggarrange(p2,p3, ncol=1, nrow=2) 


#################################
# INFLUENCE OF MAMMAL BEHAVIOUR # 
#################################

## How does specific mammal behaviour and frequency influence read counts?
## Import data on behavioural observation of captive mammals
behaviour.metadata <- read.csv("../Data/Captive_mammal_behaviour_sample_related.csv", 
                               header=TRUE)
head(behaviour.metadata)
str(behaviour.metadata)
names(behaviour.metadata)
dim(behaviour.metadata)

## First, change POOL for lynx at Wildwood Trust to T01
behaviour.metadata$Sample <- gsub("POOL", "T01", behaviour.metadata$Sample)

## Duplicate the sample column
behaviour.metadata$Enclosure <- behaviour.metadata$Sample

## Remove abbreviated location, filtration day and round from sample ID
behaviour.metadata$Enclosure <- gsub("HWP|WWT", "", behaviour.metadata$Enclosure)
behaviour.metadata$Enclosure <- gsub("-D1|-D2|-D3|-D01|-D02|-D03", "", behaviour.metadata$Enclosure)
behaviour.metadata$Enclosure <- gsub("-Round01-|-Round02-|-Round03-|-Round04-", "", behaviour.metadata$Enclosure)

## Split Sample column into multiple columns
behaviour.metadata <- behaviour.metadata %>% separate(Enclosure, into=paste("id", 1:4, sep = ""))

## Remove redundant columns
behaviour.metadata <- behaviour.metadata[,-c(7:8)]

## Rename columns
names(behaviour.metadata) <- c("Location","ID","Behaviour","Frequency",
                               "Enclosure","Sample_no")

## Rename sample numbers to reflect directed or stratified sampling strategy
behaviour.metadata$Sample_no <- gsub("P01","STR01",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("P02","STR02",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("P03","STR03",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("P04","STR04",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("P05","STR05",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("P06","STR06",behaviour.metadata$Sample_no)

behaviour.metadata$Sample_no <- gsub("T01","DIR01",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("T02","DIR02",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("T03","DIR03",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("T04","DIR04",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("T05","DIR05",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("T06","DIR06",behaviour.metadata$Sample_no)

## Change names of water vole samples
behaviour.metadata$Sample_no <- gsub("1IND","QUAR1",behaviour.metadata$Sample_no)
behaviour.metadata$Sample_no <- gsub("4IND","QUAR2",behaviour.metadata$Sample_no)

## Swap beaver and zoo in enclosure and sample number columns
behaviour.metadata$Enclosure <- gsub("ZOO", "BEAV", behaviour.metadata$Enclosure)
behaviour.metadata$Sample_no <- gsub("BEAV", "ZOO", behaviour.metadata$Sample_no)

## Merge metadata and read count data
behaviour <- merge(cap.spp, behaviour.metadata, 
                   by=c("ID","Location","Enclosure","Sample_no"))

## Check structure of dataframe
str(behaviour)

## Create new dataframe containing only variables to be modeled
## Exclude observations that correspond to designated drinking sources, 
## stratified samples, or quarantine enclosures where behaviour was not 
## observed
m3.df <- behaviour[which(!grepl("Stratified|Other", behaviour$Sample_type)),]

## Reset row names
rownames(m3.df) <- NULL

## Drop unused factor levels
m3.df <- drop.levels(m3.df)

## Check structure of new dataframe
str(m3.df)

## Plot proportional read count data against new variables
plot(m3.df$Prop_reads ~ m3.df$Location)
plot(m3.df$Prop_reads ~ m3.df$Species)
plot(m3.df$Prop_reads ~ m3.df$Behaviour)
plot(m3.df$Prop_reads ~ m3.df$Frequency)

## We do not have enough data to run a model that nests frequency of each 
## behaviour within species, even though this would be the most appropriate
## model. This because not every species displayed each behaviour.
## Instead, we will first look at a global model for behaviour with no
## nested effect/interaction, then we will examine frequency of behaviours
## on a species-by-species basis.

## Run global model for behaviour
m3 <- glmmTMB(Prop_reads ~ (1|Location/Species) + Behaviour,
              data=m3.df,
              family=binomial(link="logit"))
summary(m3)
drop1(m3, test = "Chi")

## Calculate R-squared of model
r.squaredGLMM(m3)
## marginal R-squared = 7.23% (proportion of variance explained by fixed effects)
## conditional R-squared = 9.17% (proportion of variance explained by fixed + random effects)

## Examine difference between species and between sample types
summary(glht(m3, linfct=mcp(Behaviour="Tukey")))

## Test for overdispersion
overdisp(m3)

## Plot the fitted data against the observed data
plot(m3.df$Prop_reads ~ fitted(m3))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(m3.df$Prop_reads, fitted(m3))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(m3, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.004301

## Some deviation from normality as residuals not normally distributed
## therefore model may not be that reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ m3.df$Behaviour)

## Assumption 3: no collinearity
## Only one variable being modeled so no collinearity.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Plot boxplot against bubble plot of frequency of behaviours broken down 
## by species:
m3.df$fSpecies <- factor(m3.df$Species,
                         levels=c("Castor fiber","Lutra lutra","Meles meles",
                                  "Lynx lynx","Cervus elaphus","Martes martes"))

p4a <- ggplot(m3.df, aes(x=Behaviour, y=Prop_reads))
p4a <- p4a + geom_jitter(aes(fill=fSpecies, size=Frequency),
                         colour="black", pch=21, width=0.2, alpha=0.8)
p4a <- p4a + geom_boxplot(outlier.shape=NA, alpha=0.7)
p4a <- p4a + labs(title="(a)",
                  x="Behaviour\n", y="Proportional read counts")
p4a <- p4a + scale_y_continuous(limits=c(0,1),
                                breaks=c(0,0.25,0.50,0.75,1.0))
p4a <- p4a + scale_fill_manual(name="Species",
                               labels=c(expression(italic("Castor fiber"),
                                                   italic("Lutra lutra"),
                                                   italic("Meles meles"),
                                                   italic("Lynx lynx"),
                                                   italic("Cervus elaphus"),
                                                   italic("Martes martes"))),
                               values=c("dodgerblue","dodgerblue4","greenyellow",
                                        "limegreen","green4","goldenrod2"))
p4a <- p4a + guides(fill=guide_legend(override.aes=list(size=5),
                                      label.theme=element_text(size=16)))
p4a <- p4a + theme_bw()
p4a <- p4a + theme(panel.background = element_rect(fill = 'white'),
                   panel.grid.major = element_line(colour="white"),
                   panel.grid.minor = element_line(colour="white"),
                   plot.title = element_text(face = "bold"),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black", angle=45, hjust=1),
                   axis.text.y = element_text(colour="black"),
                   strip.text.x = element_text(face="italic"),
                   text = element_text(size=20),
                   legend.position="bottom",
                   legend.box="horizontal",
                   legend.key=element_blank(),
                   legend.key.size = unit(2, 'lines'),
                   legend.text.align = 0)
p4a <- p4a + facet_grid(. ~ fSpecies, scale="free_x", space="free_x")
p4a


## Now model frequency of behaviour for behaviours that have enough data. For
## example, drinking, immersed, sniffing and swimming.
## Create dataframes for different behaviours:
drinking <- subset(m3.df, Behaviour=="Drinking")
immersed <- subset(m3.df, Behaviour=="Immersed")
sniffing <- subset(m3.df, Behaviour=="Sniffing")
swimming <- subset(m3.df, Behaviour=="Swimming")

## Drop unused factor levels
drinking <- data.frame(lapply(drinking, function(x) if(is.factor(x)) factor(x) else x))
immersed <- data.frame(lapply(immersed, function(x) if(is.factor(x)) factor(x) else x))
sniffing <- data.frame(lapply(sniffing, function(x) if(is.factor(x)) factor(x) else x))
swimming <- data.frame(lapply(swimming, function(x) if(is.factor(x)) factor(x) else x))

## Run models for each behaviour
## Drinking:
m4a <- glmmTMB(Prop_reads ~ (1|Location) + Species/Frequency,
               data=drinking,
               family=binomial(link="logit"))
summary(m4a)
drop1(m4a, test = "Chi")
r.squaredGLMM(m4a)
overdisp(m4a)
plot(drinking$Prop_reads ~ fitted(m4a))
hoslem.test(drinking$Prop_reads, fitted(m4a))
sresid <- resid(m4a, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.3053
plot(sresid ~ drinking$Species)
acf(sresid, main = "Auto-correlation plot")

## Immersed:
m4b <- glmmTMB(Prop_reads ~ (1|Location) + Species/Frequency,
               data=immersed,
               family=binomial(link="logit"))
summary(m4b)
drop1(m4b, test = "Chi")
r.squaredGLMM(m4b)
overdisp(m4b)
plot(immersed$Prop_reads ~ fitted(m4b))
hoslem.test(immersed$Prop_reads, fitted(m4b))
sresid <- resid(m4b, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.4055
plot(sresid ~ immersed$Species)
acf(sresid, main = "Auto-correlation plot")

## Sniffing:
m4c <- glmmTMB(Prop_reads ~ (1|Location) + Species/Frequency,
               data=sniffing,
               family=binomial(link="logit"))
summary(m4c)
drop1(m4c, test = "Chi")
r.squaredGLMM(m4c)
overdisp(m4c)
plot(sniffing$Prop_reads ~ fitted(m4c))
hoslem.test(sniffing$Prop_reads, fitted(m4c))
sresid <- resid(m4c, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.511
plot(sresid ~ sniffing$Species)
acf(sresid, main = "Auto-correlation plot")

## Swimming:
m4d <- glmmTMB(Prop_reads ~ (1|Location) + Species/Frequency,
               data=swimming,
               family=binomial(link="logit"))
summary(m4d)
drop1(m4d, test = "Chi")
r.squaredGLMM(m4d)
overdisp(m4d)
plot(swimming$Prop_reads ~ fitted(m4d))
hoslem.test(swimming$Prop_reads, fitted(m4d))
sresid <- resid(m4d, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.0589
plot(sresid ~ swimming$Species)
acf(sresid, main = "Auto-correlation plot")

## No significant differences in proportional read counts for any behaviour 
## exhibited by any species. This may be due to insufficient behavioural data 
## for each species. We will group specific behaviours into broader groups 
## (animal contacted water versus no contact) to determine whether there
## is any variation in proportional read counts more generally.

## Duplicate dataframe used for previous model of specific behaviours
m5.df <- m3.df

## Create new column specifying whether a behaviour involved contact
## with water or not
m5.df$Behaviour_type <- as.factor(ifelse(grepl("Drinking|Swimming|Immersed|Sniffing|Urinating", m5.df$Behaviour), "Water contact",
                                         ifelse(grepl("None|Defecating|Feeding|Walking|Standing|Resting|Grooming", m5.df$Behaviour), "No water contact", "NA")))

## Plot the new factor variable against proportional read counts
plot(m5.df$Prop_reads ~ m5.df$Behaviour_type)

## Run previous behaviour model with new category 
m5 <- glmmTMB(Prop_reads ~ (1|Location/Species) + Behaviour_type,
              data=m5.df,
              family=binomial(link="logit"))
summary(m5)
drop1(m5, test = "Chi")

## Calculate R-squared of model
r.squaredGLMM(m5)
## marginal R-squared = <1% (proportion of variance explained by fixed effects)
## conditional R-squared = 8.50% (proportion of variance explained by fixed + random effects)

## Examine difference between species and between sample types
summary(glht(m5, linfct=mcp(Behaviour_type="Tukey")))

## Test for overdispersion
overdisp(m5)

## Plot the fitted data against the observed data
plot(m5.df$Prop_reads ~ fitted(m5))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(m5.df$Prop_reads, fitted(m5))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(m5, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.006646

## Some deviation from normality as residuals not normally distributed
## therefore model may not be that reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ m5.df$Behaviour_type)

## Assumption 3: no collinearity
## Only one variable being modeled so no collinearity.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Plot boxplot of behaviour type
p4b <- ggplot(m5.df, aes(x=Behaviour_type, y=Prop_reads))
p4b <- p4b + geom_jitter(aes(fill=fSpecies), colour="black",
                         pch=21, cex=5, width=0.2, alpha=0.8)
p4b <- p4b + geom_boxplot(outlier.shape=NA, alpha=0.7)
p4b <- p4b + labs(title="(b)",
                  x="Behaviour type", y="Proportional read counts")
p4b <- p4b + scale_y_continuous(limits=c(0,1),
                                breaks=c(0,0.25,0.50,0.75,1.0))
p4b <- p4b + scale_fill_manual(name="Species",
                               labels=c(expression(italic("Castor fiber"),
                                                   italic("Lutra lutra"),
                                                   italic("Meles meles"),
                                                   italic("Lynx lynx"),
                                                   italic("Cervus elaphus"),
                                                   italic("Martes martes"))),
                               values=c("dodgerblue","dodgerblue4","greenyellow",
                                        "limegreen","green4","goldenrod2"))
p4b <- p4b + guides(fill=guide_legend(override.aes=list(size=5),
                                      label.theme = element_text(face = "italic", size=16)))
p4b <- p4b + theme(panel.background = element_rect(fill = 'white'),
                   panel.grid.major = element_line(colour="white"),
                   panel.grid.minor = element_line(colour="white"),
                   plot.title = element_text(face = "bold"),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.margin = unit(c(0,0,0.5,0),"cm"),
                   text = element_text(size=20),
                   legend.position="none")
p4b

## Plot together
ggarrange(p4a,p4b, ncol=1, nrow=2, common.legend=TRUE, legend="bottom") 


################################
# EFFECTS OF SAMPLE PROCESSING #
################################

## Now model the effect of volume of water filters and number of filters used
## on read counts
m6 <- glmmTMB(Prop_reads ~ (1|Location/Species) + Volume + No_filters,
              data=cap.spp,
              family=binomial(link="logit"))
summary(m6)
drop1(m6, test = "Chi")

## Calculate R-squared of model
r.squaredGLMM(m6)
## marginal R-squared = 8.57% (proportion of variance explained by fixed effects)
## condition R-squared = 20.57% (proportion of variance explained by fixed + random effects)

## Test for overdispersion
overdisp(m6)  # not overdispersed

## Plot the fitted data against the observed data
plot(cap.spp$Prop_reads ~ fitted(m6))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(cap.spp$Prop_reads, fitted(m6))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(m6, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.1221

## Some deviation from normality as residuals not normally distributed
## therefore model may not be that reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ cap.spp$Volume)   
plot(sresid ~ cap.spp$No_filters)   

## Assumption 3: no collinearity
plot(cap.spp$Volume ~ cap.spp$No_filters)
cor(cap.spp[,8:9])
## Borderline collinearity

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag


## PLOT MODEL FIT
## Obtain predicted values for full data set: cap.spp
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
pdat <- expand.grid(Location = levels(cap.spp$Location),
                    Species = levels(cap.spp$Species),
                    Volume = seq(min(cap.spp$Volume), max(cap.spp$Volume), length.out = 500),
                    No_filters = seq(min(cap.spp$No_filters), max(cap.spp$No_filters), length.out=2))
head(pdat)

## Make a dataframe containing the predicted data
## re.form=NA means predictions are made for random effect as a whole
## se.fit=TRUE obtains the standard error for each of the predictions
predictions <- data.frame(predict(m6, newdata=pdat, se.fit=TRUE, 
                                  na.action=na.exclude, type="response"))

## Combine predictions with the prediction data into a dataframe for plotting
preds.for.plot <- data.frame(pdat, predictions)

## Calculate and add upper/lower 95% CIs from the SEs of the predictions
preds.for.plot$ciu <- (predictions$fit+1.96*predictions$se.fit) 
preds.for.plot$cil <- (predictions$fit-1.96*predictions$se.fit)

## Check dataframe for plotting
head(preds.for.plot) 

## Collapse predictions so there is only one prediction for each value of volume
coll.preds.for.plot1 <- ddply(preds.for.plot, .(Volume), summarise,
                              fit = mean(fit),
                              ciu = mean(ciu),
                              cil = mean(cil))

## Plot predicted data against observed data
p5a <- ggplot() + ggtitle("(a)")
p5a <- p5a + geom_jitter(aes(x=cap.spp$Volume, 
                             y=cap.spp$Prop_reads, 
                             colour=cap.spp$fSpecies), cex=4, width=0.2)
p5a <- p5a + geom_line(aes(x=coll.preds.for.plot1$Volume, 
                           y=coll.preds.for.plot1$fit))
p5a <- p5a + geom_ribbon(aes(x=coll.preds.for.plot1$Volume, 
                             ymin = coll.preds.for.plot1$cil, 
                             ymax = coll.preds.for.plot1$ciu), alpha = 0.25)
p5a <- p5a + labs(x="Volume (mL)", y="Reads")
p5a <- p5a + scale_y_continuous(limits=c(0,1), 
                                breaks=c(0,0.25,0.50,0.75,1.0))
p5a <- p5a + scale_colour_manual(guide=guide_legend(label.theme=element_text(size=16)),
                                 name="Species",
                                 breaks=c("Arvicola amphibius","Castor fiber",
                                          "Lutra lutra","Erinaceus europaeus",
                                          "Meles meles","Lynx lynx","Cervus elaphus",
                                          "Sciurus vulgaris","Martes martes"),
                                 labels=c(expression(italic("Arvicola amphibius"),
                                                     italic("Castor fiber"),
                                                     italic("Lutra lutra"),
                                                     italic("Erinaceus europaeus"),
                                                     italic("Meles meles"),
                                                     italic("Lynx lynx"),
                                                     italic("Cervus elaphus"),
                                                     italic("Sciurus vulgaris"),
                                                     italic("Martes martes"))),
                                 values=c("dodgerblue","dodgerblue3","dodgerblue4",
                                          "greenyellow","olivedrab3",
                                          "limegreen","green4",
                                          "goldenrod1","goldenrod4"))
p5a <- p5a + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0, colour="black"),
                   plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
                   text = element_text(size=20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   legend.key.size = unit(2, 'lines'),
                   legend.text.align = 0)
p5a

## Collapse predictions so there is only one prediction for each value of number
## of filters
coll.preds.for.plot2 <- ddply(preds.for.plot, .(No_filters), summarise,
                              fit = mean(fit),
                              ciu = mean(ciu),
                              cil = mean(cil))

## Plot predicted data against observed data
p5b <- ggplot() + ggtitle("(b)")
p5b <- p5b + geom_jitter(aes(x=cap.spp$No_filters, 
                             y=cap.spp$Prop_reads, 
                             colour=cap.spp$fSpecies), cex=4, width=0.2)
p5b <- p5b + geom_line(aes(x=coll.preds.for.plot2$No_filters, 
                           y=coll.preds.for.plot2$fit))
p5b <- p5b + geom_ribbon(aes(x=coll.preds.for.plot2$No_filters, 
                             ymin = coll.preds.for.plot2$cil, 
                             ymax = coll.preds.for.plot2$ciu), alpha = 0.25)
p5b <- p5b + labs(x="Number of filters used", y="")
p5b <- p5b + scale_x_continuous(breaks=c(1,2))
p5b <- p5b + scale_y_continuous(limits=c(0,1), 
                                breaks=c(0,0.25,0.50,0.75,1.0))
p5b <- p5b + scale_colour_manual(guide=guide_legend(label.theme=element_text(size=16)),
                                 name="Species",
                                 breaks=c("Arvicola amphibius","Castor fiber",
                                          "Lutra lutra","Erinaceus europaeus",
                                          "Meles meles","Lynx lynx","Cervus elaphus",
                                          "Sciurus vulgaris","Martes martes"),
                                 labels=c(expression(italic("Arvicola amphibius"),
                                                     italic("Castor fiber"),
                                                     italic("Lutra lutra"),
                                                     italic("Erinaceus europaeus"),
                                                     italic("Meles meles"),
                                                     italic("Lynx lynx"),
                                                     italic("Cervus elaphus"),
                                                     italic("Sciurus vulgaris"),
                                                     italic("Martes martes"))),
                                 values=c("dodgerblue","dodgerblue3","dodgerblue4",
                                          "greenyellow","olivedrab3",
                                          "limegreen","green4",
                                          "goldenrod1","goldenrod4"))
p5b <- p5b + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0, colour="black"),
                   plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
                   text = element_text(size=20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   legend.key.size = unit(2, 'lines'),
                   legend.text.align = 0)
p5b

## Show plots side by side:
mylegend <- g_legend(p5a)
p5 <- grid.arrange(arrangeGrob(p5a + theme(legend.position="none"),
                               p5b + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2, heights=c(5, 1))



#' ---
#' 
#' ## 5) Experiment 2: data analysis
#' 
#' The primary research questions for this experiment are:
#' 
#' - How does eDNA metabarcoding compare to camera trapping for species detection?
#' - How long does eDNA from mammals persist in a natural environment?
#' 

####################################################################
# PRESENCE-ABSENCE DETECTION BY DIFFERENT METHODS AT NATURAL SITES #
####################################################################

## Replicate data frame
wild.pa <- wild.df

## Remove empty taxonomic assignments (i.e. poss. assignments from wildlife parks)
wild.pa <- wild.pa[apply(wild.pa[,-1], 1, function(x) !all(x==0)),]

## Reset row names for further indexing
rownames(wild.pa) <- NULL

## Retain only the mammal assignments
wild.pa <- wild.pa[c(4:5,8:11,23:26,31,34,36:37,40),]

## Reset row names for further indexing
rownames(wild.pa) <- NULL

## Create column specifying which lifestyle of each mammal species
wild.pa$Lifestyle <- ifelse(wild.pa$Assignment %in% wild.pa[c(1,5,8,12),1], "Semi-aquatic",
                            ifelse(wild.pa$Assignment %in% wild.pa[c(11,13),1], "Arboreal", "Ground-dwelling"))

## Move to start of dataframe
wild.pa <- wild.pa [,c(1,142,2:141)]

## Manipulate dataframe
wild.pa <- melt(wild.pa, id=c("Assignment","Lifestyle"))
colnames(wild.pa) <- c("Assignment","Lifestyle","ID","Prop_reads")

## Remove empty assignments for each sample
wild.pa[wild.pa==0] <- NA
wild.pa <- wild.pa[complete.cases(wild.pa),]
rownames(wild.pa) <- NULL

## Create column specifying which site samples belong to
wild.pa$Site <- ifelse(grepl("BAM", wild.pa$ID), "Bamff Estate",
                       ifelse(grepl("TM", wild.pa$ID), "Thorne Moors", "Tophill Low Nature Reserve"))

## Create column specifying which pond samples belong to
wild.pa$Pond <- wild.pa$ID
wild.pa$Pond <- gsub("D01-|D02-|D03-|D04-|D05-", "", wild.pa$Pond)
wild.pa <- wild.pa %>% separate(Pond, into=paste("id", 1:2, sep = ""))
colnames(wild.pa)[6:7] <- c("Pond","Sample")

## Remove location name from Pond
## Replace 01 and 02 with Pond 1 and Pond 2
wild.pa$Pond <- gsub("BAM|TM|THL", "", wild.pa$Pond)
wild.pa$Pond <- gsub("01","Pond 1", wild.pa$Pond)
wild.pa$Pond <- gsub("02","Pond 2", wild.pa$Pond)

## Now average proportional reads for assignments from Tophill Low Nature Reserve
wild.pa <- aggregate(Prop_reads ~ Assignment + Lifestyle + Site + Pond + Sample, 
                     data=wild.pa, mean)

## Convert proportional read counts for each species to presence-absence 
## data for each pond
wild.pa <- aggregate(Prop_reads ~ Assignment + Lifestyle + Site + Pond, 
                     data=wild.pa, sum)
wild.pa$Prop_reads[wild.pa$Prop_reads > 0] <- 1

## Rename columns
colnames(wild.pa) <- c("Species","Lifestyle","Site","Pond","eDNA")

## Create a new dataframe for eDNA-based species detection at site-level
eDNA.dat <- aggregate(eDNA ~ Species + Lifestyle + Site, data=wild.pa, sum)
eDNA.dat$Pond <- "Site"
eDNA.dat <- eDNA.dat[,c(1:3,5,4)]
eDNA.dat$eDNA <- ifelse(eDNA.dat$eDNA > 0, 1, 0)

## Add to dataframe containing eDNA-based species detection at pond-level
wild.pa <- rbind(wild.pa, eDNA.dat)

## Reorder dataframe
wild.pa <- wild.pa[,c(3:4,1:2,5)]

## Import camera trap and field sign data
historic.dat <- read.csv("../Data/Historical_data.csv", header=TRUE)
camera.dat <- read.csv("../Data/Camera_trap_data.csv", header=TRUE)
field.dat <- read.csv("../Data/Field_signs_data.csv", header=TRUE)

## Replace periods with space in species name
colnames(historic.dat) <- gsub("[.]", " ", colnames(historic.dat))
colnames(camera.dat) <- gsub("[.]", " ", colnames(camera.dat))
colnames(field.dat) <- gsub("[.]", " ", colnames(field.dat))

## Melt dataframes
historic.dat <- melt(historic.dat, id=c("Site","Pond"))
camera.dat <- melt(camera.dat, id=c("Site","Pond"))
field.dat <- melt(field.dat, id=c("Site","Pond"))

## Rename columns
colnames(historic.dat)[3:4] <- c("Species","Cumulative survey data")
colnames(camera.dat)[3:4] <- c("Species","Camera traps")
colnames(field.dat)[3:4] <- c("Species","Field signs")

## Create new column specifying species lifestyle
historic.dat$Lifestyle <- ifelse(historic.dat$Species %in% c("Arvicola amphibius","Castor fiber","Lutra lutra","Neovison vison"), "Semi-aquatic",
                                 ifelse(historic.dat$Species %in% c("Sciurus carolinensis","Sciurus vulgaris"), "Arboreal", "Ground-dwelling"))
camera.dat$Lifestyle <- ifelse(camera.dat$Species == "Castor fiber", "Semi-aquatic","Ground-dwelling")
field.dat$Lifestyle <- ifelse(field.dat$Species == "Castor fiber", "Semi-aquatic","Ground-dwelling")

## Reorder dataframes
historic.dat <- historic.dat[,c(1:3,5,4)]
camera.dat <- camera.dat[,c(1:3,5,4)]
field.dat <- field.dat[,c(1:3,5,4)]

## Combine camera trap and field sign data with eDNA data to compare results
method.dat <- smartbind(historic.dat, field.dat, camera.dat, wild.pa)

## Replace NAs with 0
method.dat[is.na(method.dat)] <- 0

## Melt data frame again to have a column containing method and corresponding
## data in existing columns
method.dat <- melt(method.dat, id=c("Site","Pond","Species","Lifestyle"))

## Rename columns
colnames(method.dat)[5:6] <- c("Method", "Occupancy")

## Merge data by rows
method.dat <- aggregate(Occupancy ~ Site + Pond + Species + Lifestyle + Method, 
                        data=method.dat, sum)

## Ensure Occupancy column only contains exact 0s and 1s
method.dat$Occupancy <- ifelse(method.dat$Occupancy == 0, 0, 1)

## Make sure all variables in dataframe are factors
method.dat <- data.frame(lapply(method.dat, function(x) if(is.factor(x)) factor(x) else x))
method.dat$Occupancy <- as.factor(method.dat$Occupancy)

## Create factors specifying order of factor levels for plotting
method.dat$fMethod <- factor(method.dat$Method, 
                             levels=rev(unique(method.dat$Method)))

## Remove rows containing species at sites where they have no detections
## and are not expected
method.dat <- method.dat[!(method.dat$Site == "Bamff Estate" 
                           & method.dat$Species %in% c("Muntiacus reevesi",
                                                       "Mustela erminea",
                                                       "Mustela nivalis",
                                                       "Arvicola amphibius", 
                                                       "Neovison vison",
                                                       "Lepus europaeus",
                                                       "Oryctolagus cuniculus",
                                                       "Sciurus carolinensis",
                                                       "Bos taurus")),]

method.dat <- method.dat[!(method.dat$Site == "Thorne Moors" 
                           & method.dat$Species %in% c("Sciurus vulgaris",
                                                       "Castor fiber",
                                                       "Lutra lutra",
                                                       "Neovison vison",
                                                       "Sus scrofa",
                                                       "Oryctolagus cuniculus",
                                                       "Lepus europaeus",
                                                       "Sciurus carolinensis",
                                                       "Bos taurus")),]

method.dat <- method.dat[!(method.dat$Site == "Tophill Low Nature Reserve" 
                           & method.dat$Species %in% c("Sciurus vulgaris",
                                                       "Cervus elaphus",
                                                       "Muntiacus reevesi",
                                                       "Castor fiber",
                                                       "Sus scrofa")),]

## Remove water vole data for Thorne Moors Pond 1 as species was not 
## detected by any method here
method.dat <- method.dat[!(method.dat$Site == "Thorne Moors" 
                           & method.dat$Species == "Arvicola amphibius"
                           & method.dat$Pond == "Pond 1"),]

## Create factor to order species by lifestyle
method.dat$fSpecies <- factor(method.dat$Species, 
                              levels= c("Arvicola amphibius","Neomys fodiens",
                                        "Rattus norvegicus","Neovison vison",
                                        "Castor fiber","Lutra lutra",
                                        "Sorex araneus","Myodes glareolus",
                                        "Oryctolagus cuniculus","Lepus europaeus",
                                        "Mustela erminea",
                                        "Mustela nivalis","Meles meles",
                                        "Vulpes vulpes","Muntiacus reevesi",
                                        "Capreolus capreolus","Cervus elaphus",
                                        "Sus scrofa","Canis lupus familiaris",
                                        "Sus scrofa domesticus",
                                        "Ovis aries","Bos taurus",
                                        "Pipistrellus pipistrellus",
                                        "Sciurus carolinensis","Sciurus vulgaris"))

## Create factor to colour species by lifestyle
method.dat$fLifestyle <- factor(method.dat$Lifestyle,
                                levels=c("Semi-aquatic",
                                         "Ground-dwelling",
                                         "Arboreal"))

## Plot results
hm2 <- ggplot(method.dat, aes(x=fMethod, fct_rev(as_factor(fSpecies))))
hm2 <- hm2 + geom_tile(aes(fill=Occupancy, colour=fLifestyle), size=1.5)
hm2 <- hm2 + scale_colour_manual(name="Species lifestyle",
                                 values=c("dodgerblue","limegreen","goldenrod2"))
hm2 <- hm2 + scale_fill_manual(name="Record",
                               values=c("white","black"),
                               breaks=c(0,1),
                               labels=c("Absent","Present"))
hm2 <- hm2 + labs(x="Method", y="Species")
hm2 <- hm2 + theme_bw()
hm2 <- hm2 + theme(panel.grid.major = element_line(colour="white"),
                   panel.grid.minor = element_line(colour="white"), 
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black", angle=55, vjust=1, hjust=1),
                   axis.text.y = element_text(face="italic", colour="black"),
                   text = element_text(size=20),
                   legend.key = element_rect(size = 5),
                   legend.key.size = unit(1, 'lines'))
hm2 <- hm2 + guides(fill=guide_legend(keywidth=0.5,
                                      keyheight=0.5,
                                      default.unit="inch",
                                      override.aes=list(colour="black")),
                    colour=guide_legend(keywidth=0.5,
                                        keyheight=0.5,
                                        default.unit="inch"))
hm2 <- hm2 + facet_grid(Site ~ Pond, scales = "free", space = "free")
hm2


##########################################################
# READ PROPORTIONS FOR SPECIES DETECTED AT NATURAL SITES #
##########################################################

## Replicate data frame
wild.dat <- wild.df

## Remove empty taxonomic assignments (i.e. poss. assignments from wildlife parks)
wild.dat <- wild.dat[apply(wild.dat[,-1], 1, function(x) !all(x==0)),]

## Reset row names for further indexing
rownames(wild.dat) <- NULL

## Create column specifying which group species belong to
wild.dat$Group <- ifelse(wild.dat$Assignment %in% wild.dat[c(1,3,7,12:13,15:17,19,22,27,29:30,38:39,42),1], "Bird",
                         ifelse(wild.dat$Assignment %in% wild.dat[c(2,14,18,32,35),1], "Fish",
                                ifelse(wild.dat$Assignment %in% wild.dat[c(6,20:21,28,33,41),1], "Amphibian", "Mammal")))

## Move to start of dataframe
wild.dat <- wild.dat[,c(1,142,2:141)]

## Manipulate dataframe
wild.dat <- melt(wild.dat, id=c("Assignment","Group"))
colnames(wild.dat) <- c("Assignment","Group","ID","Prop_reads")

## Remove empty assignments for each sample
wild.dat[wild.dat==0] <- NA
wild.dat <- wild.dat[complete.cases(wild.dat),]
rownames(wild.dat) <- NULL

## Create column specifying which site samples belong to
wild.dat$Site <- ifelse(grepl("BAM", wild.dat$ID), "Bamff Estate",
                        ifelse(grepl("TM", wild.dat$ID), "Thorne Moors", "Tophill Low Nature Reserve"))

## Create column specifying which pond samples belong to
wild.dat$Pond <- wild.dat$ID
wild.dat$Pond <- gsub("D01-|D02-|D03-|D04-|D05-", "", wild.dat$Pond)
wild.dat <- wild.dat %>% separate(Pond, into=paste("id", 1:2, sep = ""))
colnames(wild.dat)[6:7] <- c("Pond","Sample")

## Remove location name from Pond
## Replace 01 and 02 with Pond 1 and Pond 2
wild.dat$Pond <- gsub("BAM|TM|THL", "", wild.dat$Pond)
wild.dat$Pond <- gsub("01","Pond 1", wild.dat$Pond)
wild.dat$Pond <- gsub("02","Pond 2", wild.dat$Pond)

## Now average proportional reads for assignments from Tophill Low Nature Reserve
wild.final <- aggregate(Prop_reads ~ Assignment + Group + Site + Pond + Sample, data=wild.dat, mean)

## Order assignments by vertebrate group
wild.final <- with(wild.final, wild.final[order(Group, Assignment),])
wild.final$Assignment <- as.factor(wild.final$Assignment)
wild.final$fAssignment <- factor(wild.final$Assignment,
                                 levels=c("Bufo bufo","Rana temporaria","Pelophylax ridibundus",
                                          "Lissotriton helveticus","Lissotriton vulgaris","Triturus cristatus",
                                          "Anguilla anguilla","Gasterosteus aculeatus",
                                          "Pungitius pungitius","Salmo trutta",
                                          "Ctenopharyngodon idella","Anas",
                                          "Gallinula chloropus","Larus",
                                          "Ardea cinerea","Parus major",
                                          "Emberiza","Turdus philomelos",
                                          "Sturnus vulgaris","Garrulus glandarius",
                                          "Pica pica","Strix aluco","Buteo buteo",
                                          "Columba livia","Columba oenas",
                                          "Phasianus colchicus","Meleagris gallopavo",
                                          "Arvicola amphibius","Neomys fodiens","Rattus norvegicus",
                                          "Castor fiber","Sorex araneus","Myodes glareolus",
                                          "Oryctolagus cuniculus","Capreolus capreolus",
                                          "Cervus elaphus","Canis lupus familiaris",
                                          "Sus scrofa domesticus","Ovis aries","Bos taurus",
                                          "Pipistrellus pipistrellus",
                                          "Sciurus carolinensis"))

## Plot proportional read counts
hm3 <- ggplot(wild.final, aes(x=Sample, y=fct_rev(as_factor(fAssignment))))
hm3 <- hm3 + geom_tile(aes(fill=Prop_reads, colour=Group), size=0.5)
hm3 <- hm3 + scale_colour_manual(values=c("limegreen","goldenrod1","dodgerblue","hotpink"))
hm3 <- hm3 + scale_fill_gradient(name="Proportional\nread counts", 
                                 limits=c(0,1),
                                 breaks=c(0,0.25,0.50,0.75,1.0),
                                 low="white", high="grey50",
                                 guide=guide_colourbar(frame.colour="black",
                                                       ticks.colour="black"))
hm3 <- hm3 + geom_text(aes(label=round(Prop_reads, 3)), cex=2.7)
hm3 <- hm3 + labs(x="Sample", y="Species")
hm3 <- hm3 + theme_bw()
hm3 <- hm3 + theme(panel.grid.major = element_line(colour="white"),
                   panel.grid.minor = element_line(colour="white"), 
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(face="italic", colour="black"),
                   text = element_text(size=20),
                   legend.key.size = unit(1, 'lines'))
hm3 <- hm3 + facet_grid(Site ~ Pond, scales = "free", space = "free")
hm3


###################################################
# SPATIOTEMPORAL VARIATION IN MAMMAL eDNA SIGNALS #
###################################################

## Create dataframe for Tophill Low experiment
THL.exp <- wild.df[,which(grepl("Assignment|THL", colnames(wild.df)))]

## Remove empty taxonomic assignments as these are problematic for vegan
THL.exp <- THL.exp[apply(THL.exp[,-1], 1, function(x) !all(x==0)),]

## Reset row names
rownames(THL.exp) <- NULL

## Create column specifying which group species belong to
THL.exp$Group <- ifelse(THL.exp$Assignment %in% THL.exp[c(1,3,8:9,11,13,17,19:20,26),1], "Bird",
                        ifelse(THL.exp$Assignment %in% THL.exp[c(2,10),1], "Fish",
                               ifelse(THL.exp$Assignment %in% THL.exp[c(5,12,18,25),1], "Amphibian", "Mammal")))

## More new column to start of dataframe
THL.exp <- THL.exp[,c(1,102,2:101)]

## Melt dataframe
THL.exp <- melt(THL.exp, id=c("Assignment","Group"))

## Rename new columns
colnames(THL.exp)[3:4] <- c("ID", "Prop_reads")

## Remove empty assignments for each sample
THL.exp[THL.exp==0] <- NA
THL.exp <- THL.exp[complete.cases(THL.exp),]
rownames(THL.exp) <- NULL

## Move ID to start of dataframe
THL.exp <- THL.exp[,c(3,1:2,4)]

## Split ID column in multiple columns containing Day, Pond and Sample
THL.exp$temp <- THL.exp$ID
THL.exp <- THL.exp %>% separate(temp, into=paste("id", 1:3, sep = ""))

## Move new columns to start of dataframe
THL.exp <- THL.exp[,c(1,5:7,2:4)]

## Rename new columns
colnames(THL.exp)[2:4] <- c("Day","Pond","Sample")

## Sort dataframe by Group
THL.exp <- with(THL.exp, THL.exp[order(Group, Assignment),])

## Create new column specifying order that factor levels are to be plotted in
THL.exp$Assignment <- as.factor(THL.exp$Assignment)
THL.exp$fAssignment <- factor(THL.exp$Assignment,
                              levels=c("Bufo bufo","Pelophylax ridibundus",
                                       "Lissotriton vulgaris","Triturus cristatus",
                                       "Anguilla anguilla","Gasterosteus aculeatus",
                                       "Anas","Gallinula chloropus",
                                       "Larus","Ardea cinerea",
                                       "Parus major","Emberiza",
                                       "Turdus philomelos","Pica pica",
                                       "Phasianus colchicus","Meleagris gallopavo",
                                       "Neomys fodiens","Rattus norvegicus",
                                       "Sorex araneus","Oryctolagus cuniculus",
                                       "Capreolus capreolus","Canis lupus familiaris",
                                       "Sus scrofa domesticus","Ovis aries",
                                       "Bos taurus","Sciurus carolinensis"))

## Plot spatiotemporal variation in mammal proportional read counts
hm4 <- ggplot(THL.exp, aes(x=Sample, 
                           y=fct_rev(as_factor(fAssignment)), 
                           fill=Prop_reads))
hm4 <- hm4 + geom_tile(aes(colour=Group), size=0.5)
hm4 <- hm4 + scale_colour_manual(values=c("limegreen","goldenrod1","dodgerblue","hotpink"))
hm4 <- hm4 + scale_fill_gradient(name="Proportional\nread counts",
                                 limits=c(0,1.0),
                                 breaks=c(0,0.25,0.50,0.75,1.0),
                                 low="white", high="grey50", 
                                 guide=guide_colourbar(frame.colour="black",
                                                       ticks.colour="black"))
hm4 <- hm4 + geom_text(aes(label=round(Prop_reads, 3)), cex=2.5)
hm4 <- hm4 + labs(x="Sample", y="Species")
hm4 <- hm4 + theme_bw()
hm4 <- hm4 + theme(panel.grid.major = element_line(colour="white"),
                   panel.grid.minor = element_line(colour="white"), 
                   plot.title = element_text(colour="black", face="bold", hjust=0),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(face="italic", colour="black"),
                   text = element_text(size=20),
                   legend.key = element_rect(fill = alpha("white", 0.0)),
                   legend.key.size = unit(1, 'lines'))
hm4 <- hm4 + facet_grid(Pond ~ Day,  scale="free_x", space="free_x")
hm4

