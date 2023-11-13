##########################################################################################
## Signor-Lipps corrected estimates of extinction of dwarf hippos & elephants in Cyprus ##
## Corey Bradshaw
## November 2023
##########################################################################################

library(dplyr) 
library(Rexinct)

# source functions
source("qualityRating.R")
source("EndRating.R")

## dwarf hippo Phanourios minor (from Zazzo et al. 2015 PLoS One doi:10.1371/journal.pone.0134429)
phanourios <- read.table("phanourios.txt", header = T, sep="\t")
phanourios[is.na(phanourios)] <- 'na' #replaces all the missing data with "na"to avoid TRUE/FALSE errors
phanourios$C14_CNRatioValue <-as.numeric(as.character(phanourios$C14_CNRatioValue)) #makes CNRatioValue numeric
phanourios$C14_CNRatioValue[is.na(phanourios$C14_CNRatioValue)] <- '0' #replaces missing data for CNRatioValue with 0
phanourios$C14_NPercentage <-as.numeric(as.character(phanourios$C14_NPercentage)) #makes NPercentage nurmeric
phanourios$C14_NPercentage[is.na(phanourios$C14_NPercentage)] <- '0' #replaces missing data for NPercentage with 0

# quality rating
phanourios.rated <- qualityRating(phanourios)
phanourios.endrated <- EndRating(phanourios.rated)

phanourios.merged <- merge(phanourios, phanourios.endrated, by="AgeID")

phanourios.out <- data.frame(phanourios.merged$AgeID, phanourios.merged$Age, phanourios.merged$Precision, phanourios.merged$preQuality)
colnames(phanourios.out) <- c("ID", "age", "err", "rating")
phanourios.out
phanourios.good <- subset(phanourios.out, rating!="C")
phanourios.good
phanourios.sort <- phanourios.good[order(phanourios.good[,2], decreasing=F),]
phanourios.sort

phanouriosTS <- data.frame(phanourios.sort$age, phanourios.sort$err)

# save time series to working directory
write.table(phanouriosTS, file = "phanouriosTS.txt", row.names = FALSE, col.names = FALSE)

# CRIWM
criwm(chrono_data = "phanouriosTS.txt", signor_lipps = "ext", biased=T, radiocarbon="all", cal_curve = "intcal20", cal_save=T, criwm_save = T)


# dwarf elephant Palaeoloxodon cypriotes (from Wigand and Simmons 1999. The dating of Akrotiri <em>Aetokremnos</em>, in
# <em>Faunal Extinction in an Island Society. Pygmy Hippopotamus Hunters of Cyprus</em>. A.H. Simmons (ed). 
# Kluwer Academic Publishers, New York. pp. 193-215)
palaeoloxodon <- read.table("palaeoloxodon.txt", header = T, sep="\t")
palaeoloxodon[is.na(palaeoloxodon)] <- 'na' #replaces all the missing data with "na"to avoid TRUE/FALSE errors
palaeoloxodon$C14_CNRatioValue <-as.numeric(as.character(palaeoloxodon$C14_CNRatioValue)) #makes CNRatioValue numeric
palaeoloxodon$C14_CNRatioValue[is.na(palaeoloxodon$C14_CNRatioValue)] <- '0' #replaces missing data for CNRatioValue with 0
palaeoloxodon$C14_NPercentage <-as.numeric(as.character(palaeoloxodon$C14_NPercentage)) #makes NPercentage nurmeric
palaeoloxodon$C14_NPercentage[is.na(palaeoloxodon$C14_NPercentage)] <- '0' #replaces missing data for NPercentage with 0

# quality rating
palaeoloxodon.rated <- qualityRating(palaeoloxodon)
palaeoloxodon.endrated <- EndRating(palaeoloxodon.rated)

palaeoloxodon.merged <- merge(palaeoloxodon, palaeoloxodon.endrated, by="AgeID")
palaeoloxodon.merged$qualRating <- ifelse(palaeoloxodon.merged$preRating == "C", "C", palaeoloxodon.merged$preQuality)

palaeoloxodon.out <- data.frame(palaeoloxodon.merged$AgeID, palaeoloxodon.merged$Age*1000, palaeoloxodon.merged$Precision*1000, palaeoloxodon.merged$qualRating)
colnames(palaeoloxodon.out) <- c("ID", "age", "err", "rating")
palaeoloxodon.out
palaeoloxodon.good <- subset(palaeoloxodon.out, rating=="B" | rating=="A")
palaeoloxodon.good
palaeoloxodon.sort <- palaeoloxodon.good[order(palaeoloxodon.good[,2], decreasing=F),]
palaeoloxodon.sort

palaeoloxodonTS <- data.frame(palaeoloxodon.sort$age, palaeoloxodon.sort$err)

# save time series to working directory
write.table(palaeoloxodonTS, file = "palaeoloxodonTS.txt", row.names = FALSE, col.names = FALSE)

# CRIWM
criwm(chrono_data = "palaeoloxodonTS.txt", signor_lipps = "ext", biased=T, radiocarbon=0, cal_save=T, criwm_save = T)
