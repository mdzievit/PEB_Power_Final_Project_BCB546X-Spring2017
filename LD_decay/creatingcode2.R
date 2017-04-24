setwd("C:/Users/mcnellie/Desktop/BCB_group") #setting working directory
#loading required packages
library(LDheatmap)
library(genetics)

start_time <- proc.time() #for use in calculating how long loop takes 

#reading in snp data
snp <- read.delim("C:/Users/mcnellie/Desktop/BCB_group/PEB_Power_Final_Project/Jinyu/genetic_map_data_EXPIM2012_wh.txt",header=T, stringsAsFactors = F) #physical map
snp[snp=="--"] <- NA #setting missing SNP data equal to NA
#reading in what I call the 'key' contains data needed to connect genotypes with market class (processing, fresh market, vintage, etc)
key <- read.delim("C:/Users/mcnellie/Desktop/BCB_group/PEB_Power_Final_Project/Supplemental_Data/Table_S1.txt", header=T, stringsAsFactors = F) #file that will be used to determine which genotypes belong in which catagoryu
colnames(key) <- key[1,]; key <- key[-1,] #formatting the 'key' df, setting column names
subpops <- c("Processing", "Fresh market", "Vintage") #listing the subpopulations that LD will be calculated for
final_out <- data.frame() #empty dataframe that will hold final data

#starting loop, will go through each market class named in the vector 'subpops'
for(i in 1:3){
  sub_genos <- key$`SolCAP T-number`[key$`Market Class` == subpops[i]] #getting the names of genotypes that are in the market class
  columnsToKeep <- which(colnames(snp) %in% sub_genos) #identifying the columns in the 'snp' dataframe that contain the genotypes of interest
  sub_snp <- snp[,c(1:3, columnsToKeep)] #subsetting 'snp' dataframe so it contains only the genotypes of interest
  sub_out <- data.frame() #empty data that will hold the output of each individual subpopulation
  
  #looping through each chromosome
  for(j in 1:12){
    sub_snp_chr <- sub_snp[sub_snp$Chromosome == j,] #subsetting for the chromosome of interest
    sub_snp_chr <- sub_snp_chr[order(sub_snp_chr$Position..cM.),] #ordering rows based upon marker position (cM)
    
    temp_new_snp <- data.frame() #temporary df that will hold the output of the following loop
    #loop will calculate the minor alelle frequency for each marker, need to drop markers with a MAF >= 10
    for(k in 1:nrow(sub_snp_chr)){
      snp_data <- sub_snp_chr[k,-c(1:3)] #pulling out only snp calls for each row
      snp_table <- table(unlist(snp_data)) #obtaining a count of each marker type
      snp_p_table <- prop.table(snp_table) #obtaining percentage of each marker type
      if(length(snp_p_table) == 2){ #if there are two unique alleles
        if(max(snp_p_table) < 0.90){
          temp_new_snp <- rbind.data.frame(temp_new_snp, sub_snp_chr[k,]) #if most frequent allele occurs less than 90% of time, save that row
        }
      }
      if(length(snp_p_table) == 3){ #if there are three unique alleles, meaning heterozygotes
        if(max(snp_p_table) + 0.5*(snp_p_table[2]) < 0.90){
          temp_new_snp <- rbind.data.frame(temp_new_snp, sub_snp_chr[k,]) #if most frequent allele occurs less than 90% of time, save that row
        }
      }
    }
    sub_snp_chr <- temp_new_snp #renaming the df created from the above loop
    genpos_cM <- sub_snp_chr$Position..cM. #saving the genetic position of each snp
    sub_snp_chr_org <- sub_snp_chr #saving orginal data for that subpopulation and chromosome, will be used later
    sub_snp_chr <- sub_snp_chr[,-c(1:3)] #removing columns that do not contain SNP information
    sub_snp_chr <- t(sub_snp_chr) #transposing, required for later steps
    sub_snp_chr <- paste(substring(sub_snp_chr,1,1), substring(sub_snp_chr,2,2), sep="/") #changing format of SNP calls from 'AA' to 'A/A'
    sub_snp_chr <- as.factor(sub_snp_chr) #going from charecter to factor
    
    sub_snp_chr <- as.data.frame(matrix(data=sub_snp_chr, ncol=nrow(sub_snp_chr_org))) #somewhere along the lines the SNP calls became a vector, I think during the change of format, this converts it back into a dataframe 
    sub_snp_chr[sub_snp_chr=="NA/NA"] <- NA #NAs became misnamed during the changing of format 
    colnames(sub_snp_chr) <- sub_snp_chr_org[,1] #setting column names equal to marker names. Since it has been transposed, the markers were orginally rows, but are now columns

    #the SNP calls are currently factors, need to be 'genotype', using command as.genotype from 'genetics' package
    num <- ncol(sub_snp_chr) #number of columns that need to be looped over
    for(l in 1:num){
      sub_snp_chr[,l] <- as.genotype(sub_snp_chr[,l]) #converting snp calls from factors to 'genotypes'
    }

    ld <- LDheatmap(sub_snp_chr, add.map=F, distances = "genetic", genetic.distances = genpos_cM, LDmeasure="r") #using package 'LDHeatmap' to calculate the r^2 value (LDmeasure='r'), using genetic distances ()
    temp_matrix <- as.matrix(ld$LDmatrix) #extracting the LD data
    rowCol <- expand.grid(rownames(temp_matrix), colnames(temp_matrix)) #the temp_matrix has the marker names for the row and column names, extracting these to be used to create a DF of pairwise marker LD
    labs <- rowCol[as.vector(upper.tri(temp_matrix,diag=F)),] #setting labels
    ld_df <- cbind(labs, temp_matrix[upper.tri(temp_matrix,diag=F)]) #creating dataframe with the pairwise markers and their LD
    dev.off()#the function 'LDheatmap' creates a graph, turning that graph off
    
    #need to loop through and calculate genetic distance between each set of pairwise markers
    for(m in 1:nrow(ld_df)){
      p1 <- sub_snp_chr_org$Position..cM.[sub_snp_chr_org$SNP.markers==ld_df[m,1]] #extracting genetic position of first marker
      p2 <- sub_snp_chr_org$Position..cM.[sub_snp_chr_org$SNP.markers==ld_df[m,2]] #extracting genetic position of second marker
      ld_df[m,4] <- round(abs(p1-p2),2) #getting the absolute difference between the two rounded to two decimal places
    }
    ld_df$market <- subpops[i] #adding a column listing what market class it is
    ld_df$chr <- j #adding a column listing what chromosome LD was just calculated for
    sub_out <- rbind.data.frame(sub_out, ld_df) #combing the LD results for individual chromosomes together
  }
  write.csv(sub_out, paste("subpop", subpops[i], "LD.csv", sep="_"), row.names=F) #saving results for the individual subpopulations
  final_out <- rbind.data.frame(final_out, sub_out) #combining the results of all subpopulations together
}
write.csv(final_out, "all_subpops_LD_combined.csv", row.names = F) #saving the LD of all subpopulations

end_time <- proc.time() - start_time; print(end_time) #calculates and prints (displays) how long the loop took to run
#2125 or 35min


