library(adegenet)
library(tidyr)
library(dplyr)

setwd("C:/Users/mdzievit.IASTATE/Dropbox/Classes/EEOB_546X/PEB_Power_Final_Project/Structure_Analysis/")
data <- read.delim("genetic_map_data_EXPIM2012_wh.txt",header = TRUE)
dim(data)
colnames(data)[1:5]

genotype <- read.delim("Table_S1.txt",header = TRUE, skip = 1)
colnames(genotype)
genotype <- genotype %>% select(SolCAP.T.number,Market.Class)

##Identified incorrect classifications in the genotype file related to the original PCA figure. Need to correct these to make the right plot. Fresh market/Processing is one that was changed after. Vintage/processing were ones that were changed after analysis. I am changing these to their correct classification to remove the type. It causes problems down stream.
levels(genotype$Market.Class)
genotype[genotype == "Fresh market/Processing"] <- "Processing"
genotype[genotype == "Vintage/Processing"] <- "Processing"
genotype <- droplevels(genotype)
levels(genotype$Market.Class)

data_t <- t(data)

##Adding the columns and removing the 3 lines (markers, chr, cm)
colnames(data_t) <- data_t[1,]
data_t <- data_t[-(1:3),]

##Converting to dataframe
data_t <- as.data.frame(data_t)

##Adding the row names to the same column name as the genotype file
data_t$SolCAP.T.number <- rownames(data_t)
rownames(data_t) <- NULL

##Checking the levels of the market class
levels(genotype$Market.Class)

##Joining the market class info, filtering out hyrbids, dropping the market class, and transforming the table back to the original way
data_t_info <- genotype %>% 
  select(SolCAP.T.number, Market.Class) %>% 
  full_join(data_t, by = c("SolCAP.T.number")) %>%
  filter(Market.Class != "Hybrid") %>%
  select(everything(),-Market.Class)

data_t_info.proc <- genotype %>% 
  select(SolCAP.T.number, Market.Class) %>% 
  full_join(data_t, by = c("SolCAP.T.number")) %>%
  filter(Market.Class != "Hybrid" & 
           Market.Class == "Fresh market" | 
           Market.Class == "Processing" | 
           Market.Class == "Processing" |
           Market.Class == "Vintage") %>%
  select(everything(),-Market.Class)

##Converting the full dataset into a ped file here:
##Making the first line the header and removing that line
data_filt_info <- t(data_t_info)
data_filt_info <- cbind(rownames(data_filt_info),data_filt_info)
rownames(data_filt_info) <- NULL
data_filt_info <- as.data.frame(data_filt_info)
colnames(data_filt_info) <- unlist(data_filt_info[1,])
data_filt_info <- data_filt_info[-1,]

##Checking the col names, need to rename the marker column name
colnames(data_filt_info)[1:3]
colnames(data)[1:3]
colnames(data_filt_info)[1] <- "SNP.markers"

data.filt.info.pos <- data %>% 
  select(SNP.markers,Chromosome,Position..cM.) %>%
  right_join(data_filt_info,by = ("SNP.markers"))

data.filt.info.pos <- data.filt.info.pos %>% arrange(Chromosome,Position..cM.)
data.mat <- t(as.matrix(data.filt.info.pos))
colnames(data.mat) <- data.mat[1,]
data.mat <- data.mat[-(1:3),]
data.mat <- as.data.frame(data.mat)

data.mat.alleles <- df2genind(data.mat, ploidy = 2, sep = "")
data.mat.alleles <- genind2df(data.mat.alleles, oneColPerAll = TRUE)
data.mat.alleles[data.mat.alleles == "-"] <- 0
data.mat.alleles <- as.matrix(data.mat.alleles)
rownames(data.mat.alleles) <- NULL

names <- colnames(data.mat)
names <- cbind(1:length(names),names)
names <- rbind(names,names)
names <- as.data.frame(names)
names <- names[sort(names$V1),]
colnames(data.mat.alleles) <- names$names

data.sampleID <- rownames(data.mat)
data.famID <- c(1:length(data.sampleID))
data.famID <- paste("FAM",data.famID,sep = "")
data.patID <- rep(0,length(data.sampleID))
data.matID <- rep(0,length(data.sampleID))
data.sex <- rep(0,length(data.sampleID))
data.affection <- rep(0,length(data.sampleID))
data.snpmatrix <- cbind(data.famID,data.sampleID,data.patID,data.matID,data.sex,
                   data.affection,data.mat.alleles)

data.snps.orig <- data[1:3] %>% select(everything(),names = SNP.markers) %>% 
  mutate(BP = 0) %>%
  arrange(Chromosome,Position..cM.) %>%
  select(Chromosome,names,Position..cM., BP)

data.snps <- data.snps.orig %>%
  mutate(BP = 1:dim(data.snps.orig)[1], Position..cM. = 0,
         names = 1:dim(data.snps.orig)[1]) %>% 
  select(Chromosome,names,Position..cM., BP)


data.sampleID.info <- as.data.frame(data.sampleID)
colnames(data.sampleID.info) <- colnames(genotype)[1]
data.sampleID.info <- data.sampleID.info %>% 
  left_join(genotype,by = "SolCAP.T.number") %>%
  select(SolCAP.T.number,Market.Class)

write.table(data.snpmatrix,file = "tomato_SNPs.ped",row.names = FALSE, 
            col.names = FALSE,sep = " ", quote = FALSE)
write.table(data.snps, file = "tomato_SNPs.map", row.names = FALSE,
            col.names = FALSE,sep = " ", quote = FALSE)
write.table(data.sampleID,file = "tomato_sampleID.txt", row.names = FALSE, 
            col.names = FALSE,sep = " ", quote = FALSE)
write.table(data.sampleID.info,file = "tomato_sampID_info.txt",row.names = FALSE,
            col.names = TRUE, sep = "\t", quote = FALSE)

################################################################################
##Process the reduced file to create a PED file
data_filt_info <- t(data_t_info.proc)
data_filt_info <- cbind(rownames(data_filt_info),data_filt_info)
rownames(data_filt_info) <- NULL
data_filt_info <- as.data.frame(data_filt_info)
colnames(data_filt_info) <- unlist(data_filt_info[1,])
data_filt_info <- data_filt_info[-1,]

##Checking the col names, need to rename the marker column name
colnames(data_filt_info)[1:3]
colnames(data)[1:3]
colnames(data_filt_info)[1] <- "SNP.markers"

data.filt.info.pos <- data %>% 
  select(SNP.markers,Chromosome,Position..cM.) %>%
  right_join(data_filt_info,by = ("SNP.markers"))

data.filt.info.pos <- data.filt.info.pos %>% arrange(Chromosome,Position..cM.)
data.mat <- t(as.matrix(data.filt.info.pos))
colnames(data.mat) <- data.mat[1,]
data.mat <- data.mat[-(1:3),]
data.mat <- as.data.frame(data.mat)

data.mat.alleles <- df2genind(data.mat, ploidy = 2, sep = "")
data.mat.alleles <- genind2df(data.mat.alleles, oneColPerAll = TRUE)
data.mat.alleles[data.mat.alleles == "-"] <- 0
data.mat.alleles <- as.matrix(data.mat.alleles)
rownames(data.mat.alleles) <- NULL

names <- colnames(data.mat)
names <- cbind(1:length(names),names)
names <- rbind(names,names)
names <- as.data.frame(names)
names <- names[sort(names$V1),]
colnames(data.mat.alleles) <- names$names

data.sampleID <- rownames(data.mat)
data.famID <- c(1:length(data.sampleID))
data.famID <- paste("FAM",data.famID,sep = "")
data.patID <- rep(0,length(data.sampleID))
data.matID <- rep(0,length(data.sampleID))
data.sex <- rep(0,length(data.sampleID))
data.affection <- rep(0,length(data.sampleID))
data.snpmatrix <- cbind(data.famID,data.sampleID,data.patID,data.matID,data.sex,
                        data.affection,data.mat.alleles)

data.snps <- data.snps.orig %>%
  mutate(BP = 1:dim(data.snps.orig)[1], Position..cM. = 0,
         names = 1:dim(data.snps.orig)[1]) %>% 
  select(Chromosome,names,Position..cM., BP)

data.sampleID.info <- as.data.frame(data.sampleID)
colnames(data.sampleID.info) <- colnames(genotype)[1]
data.sampleID.info <- data.sampleID.info %>% 
  left_join(genotype,by = "SolCAP.T.number") %>%
  select(SolCAP.T.number,Market.Class)

write.table(data.snpmatrix,file = "tomato_SNPs_proc.ped",row.names = FALSE, 
            col.names = FALSE,sep = " ", quote = FALSE)
write.table(data.snps, file = "tomato_SNPs_proc.map", row.names = FALSE,
            col.names = FALSE,sep = " ", quote = FALSE)
write.table(data.sampleID,file = "tomato_sampleID_proc.txt", row.names = FALSE, 
            col.names = FALSE,sep = " ", quote = FALSE)
write.table(data.sampleID.info,file = "tomato_sampID_proc_info.txt",row.names = FALSE,
            col.names = TRUE, sep = "\t", quote = FALSE)

################################################################
##Input the data into PLINK to convert into .bed files
##Send files to linux server to run on admixture w/ 5 and 7
##subpopulations for both the full and reduced dataset
################################################################