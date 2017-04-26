#This script is a long, inefficient way of creating the graphs. On git is another R script that is more compact, but contains an error in formatting. See the README file for more details
#setting working directory and loading required libraries
setwd("C:/Users/jmcne/Dropbox/Spring 2017/BCB 546X/project")
library(ggplot2)

#importing LD data created in previous script
LD <- read.csv("all_subpops_LD_combined.csv", stringsAsFactors = F)
#setting column names
colnames(LD) <- c("SNP1", "SNP2", "LD", "Distance", "Subpop", "Chr")
#changing name of Fresh market subpop, removes the " " between words, replaces with an "_"
LD$Subpop[LD$Subpop=="Fresh market"] <- "Fresh_market"

#Only interested in pairs of SNPs that are within cM (genetic distance) of each other
LD <- LD[LD$Distance <= 50,]
#Creating 3 dataframes, each containing only a single subpopulation/market type. 
proc <- LD[LD$Subpop=="Processing",]
fm <- LD[LD$Subpop=="Fresh_market", ]
vin <- LD[LD$Subpop=="Vintage",]

#As mentioned above, this code is inefficient. It has the same thing repeated 36 times, with different names.
#3 subpopulations * 12 chromosomes = 36 different graphs
#For the sake of brevity I will copy and paste, and annotate the code for the first graph, proc1 (processing subpop, chromosome 1).


proc1 <- proc[proc$Chr==1,] #Subsetting the dataframe proc, which only contains LD data for processing market types, for chromosome 1

#1 # '1' denotes the code below is for chromosome 1

# Below chunk of code is used to calculate the non-linear regression trendlines for the market type and chromosome being analyzed. 
# The output of the code, fpoints_proc1, is used in the ggplot code that immediately follows

distance <- proc1$Distance #extacting and saving genetic distance between marker pairs (for processing cultivars on chromosome 1)
LD_data <- proc1$LD #extracting and saving LD between marker pairs
n <- nrow(proc1) #the number of marker pairs
HW_st<-c(C=0.3) #A required value for the code below, do not know what it means or how it influences results

#Below is the code/formula to calculate the NLR

temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))

nlr_summary<-summary(temp_nlr) #saving summary data
new_rho<-nlr_summary$parameters[1] #saving estimation of population recombination parameter, for use in code below

#Generates new estimates (fpoints_proc1) of LD that correspond to genetic distance. fpoints_proc1 will be used to add a trendline in the ggplot

fpoints_proc1 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))


#generates and saves a ggplot item(?). THe output, proc1_gg, can be called again in the future to generate a graph as long as the 
#the following data is loaded into the R environment: proc1 and fpoints_proc1

proc1_gg<- ggplot(data=proc1, aes(x=proc1$Distance, y=proc1$LD)) + #telling ggplot which data points to look at
  geom_point(col="grey",size=.4) + #adds grey points of size .4 at the coordinates passed in previous line of code
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) + #adds a red LOESS curve, line width of 0.5, span (smoothing parameter) of 0.9
  geom_line(aes(x=proc1$Distance, y=fpoints_proc1), col="blue", size=.5) + #adds a blue NLR curve
  theme_bw() + #white background, grey background is the default
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), #what settings the plot title will have, size=font size, face='bold' means bold font. 
    #Starting line for adjusting the 'theme' of the plot
    legend.title=element_blank(), legend.position = "none", #no legend information
    axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), #setting y axis label parameters, rotated 90 degrees, horiztonally adjusted
    axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), #setting x axis label parameters
    axis.title.x = element_text(size=8), axis.title.y = element_text(size=8), #setting x and y axis title parameters
    axis.ticks.length=unit(1.5, "mm"), #adjusting tick length on the axis'
    plot.margin=unit(c(2,5,2,2),"mm"), #adjusting whole plot margins
    panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #setting the panel border  and grid parameters
    axis.line = element_line(colour = "black"), #Setting color of axis lines
    panel.background = element_rect(colour = "black", fill=NA, size=.5)) + #makes a black border go all around all sides of plot, default is only a line for x and y axis, nothing to the top and right
    #End theme parameters
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc1$Chr), sep=" "))) + #setting x & y, and title labels (actual text that will be produced)
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) + #changing the number of tick marks that appear along the y axis
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) + #changing the number of tick marks that appear along the x axis
  geom_hline(yintercept = 0.23, linetype="dashed") + #adding a dashed horizontal line at y=0.23
  geom_hline(yintercept = 0.2) #adding a solid horizontal line at y=0.2

###################################
#the code below appears again about 1000 lines later.
#It is used to generate plots that contain 12 images. 3 sub groups, 4 chromosomes. 
##################################
#loading packages required
#library(grid)
#library(gridExtra)
library("gridExtra", lib.loc="~/R/win-library/3.2")
library("grid", lib.loc="C:/Program Files/R/R-3.2.3/library")



#arrangeGrob generates an image that is not displayed in the Rstudio plot viewer
#the code below will paste together the images listed (proc1_gg then fm1_gg, etc), with 3 columns. So after vin1_gg a new row starts. For chromosomes 1 to 4, for the the three subpopulations of interesting (processing, fresh market and vintage)
g1 <- arrangeGrob(proc1_gg, fm1_gg, vin1_gg,proc2_gg, fm2_gg, vin2_gg, proc3_gg, fm3_gg, 
                  vin3_gg, proc4_gg, fm4_gg, vin4_gg, ncol=3)
ggsave(file="chr1_4.pdf", g1, width = 8, height = 12, units=c("in")) #saves the image 'created' by arrangeGrob as a pdf, 8x12 inches in size

#Creates an image for chromosomes 5 to 8
g2 <- arrangeGrob(proc5_gg, fm5_gg, vin5_gg,proc6_gg, fm6_gg, vin6_gg, proc7_gg, fm7_gg, vin7_gg, proc8_gg, fm8_gg, vin8_gg, ncol=3)
ggsave(file="chr5_8.pdf", g2, width = 8, height = 12, units=c("in")) #saving the image

#Creates an image for chromosomes 9 to 12
g3 <- arrangeGrob(proc9_gg, fm9_gg, vin9_gg,proc10_gg, fm10_gg, vin10_gg, proc11_gg, fm11_gg, vin11_gg, proc12_gg, fm12_gg, vin12_gg, ncol=3)
ggsave(file="chr9_12.pdf", g3, width = 8, height = 12, units=c("in")) #saving the image



####################processing####################
proc1 <- proc[proc$Chr==1,]
proc2 <- proc[proc$Chr==2,]
proc3 <- proc[proc$Chr==3,]
proc4 <- proc[proc$Chr==4,]
proc5 <- proc[proc$Chr==5,]
proc6 <- proc[proc$Chr==6,]
proc7 <- proc[proc$Chr==7,]
proc8 <- proc[proc$Chr==8,]
proc9 <- proc[proc$Chr==9,]
proc10 <- proc[proc$Chr==10,]
proc11 <- proc[proc$Chr==11,]
proc12 <- proc[proc$Chr==12,]


#1
distance <- proc1$Distance
LD_data <- proc1$LD
n <- nrow(proc1)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc1 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc1_gg<- ggplot(data=proc1, aes(x=proc1$Distance, y=proc1$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc1$Distance, y=fpoints_proc1), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc1$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#2

distance <- proc2$Distance
LD_data <- proc2$LD
n <- nrow(proc2)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc2 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc2_gg<- ggplot(data=proc2, aes(x=proc2$Distance, y=proc2$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc2$Distance, y=fpoints_proc2), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc2$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#3
distance <- proc3$Distance
LD_data <- proc3$LD
n <- nrow(proc3)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc3 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc3_gg<- ggplot(data=proc3, aes(x=proc3$Distance, y=proc3$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc3$Distance, y=fpoints_proc3), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc3$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#4
distance <- proc4$Distance
LD_data <- proc4$LD
n <- nrow(proc4)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc4 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc4_gg<- ggplot(data=proc4, aes(x=proc4$Distance, y=proc4$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc4$Distance, y=fpoints_proc4), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc4$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#5
distance <- proc5$Distance
LD_data <- proc5$LD
n <- nrow(proc5)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc5 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc5_gg<- ggplot(data=proc5, aes(x=proc5$Distance, y=proc5$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc5$Distance, y=fpoints_proc5), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc5$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#6
distance <- proc6$Distance
LD_data <- proc6$LD
n <- nrow(proc6)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc6 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc6_gg<- ggplot(data=proc6, aes(x=proc6$Distance, y=proc6$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc6$Distance, y=fpoints_proc6), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc6$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#7
distance <- proc7$Distance
LD_data <- proc7$LD
n <- nrow(proc7)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc7 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc7_gg<- ggplot(data=proc7, aes(x=proc7$Distance, y=proc7$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc7$Distance, y=fpoints_proc7), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc7$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#8


distance <- proc8$Distance
LD_data <- proc8$LD
n <- nrow(proc8)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc8 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc8_gg<- ggplot(data=proc8, aes(x=proc8$Distance, y=proc8$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc8$Distance, y=fpoints_proc8), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc8$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)



#9
distance <- proc9$Distance
LD_data <- proc9$LD
n <- nrow(proc9)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc9 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc9_gg<- ggplot(data=proc9, aes(x=proc9$Distance, y=proc9$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc9$Distance, y=fpoints_proc9), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc9$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#10

distance <- proc10$Distance
LD_data <- proc10$LD
n <- nrow(proc10)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc10 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc10_gg<- ggplot(data=proc10, aes(x=proc10$Distance, y=proc10$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc10$Distance, y=fpoints_proc10), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc10$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#11

distance <- proc11$Distance
LD_data <- proc11$LD
n <- nrow(proc11)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc11 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc11_gg<- ggplot(data=proc11, aes(x=proc11$Distance, y=proc11$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc11$Distance, y=fpoints_proc11), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc11$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)


distance <- proc12$Distance
LD_data <- proc12$LD
n <- nrow(proc12)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_proc12 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



proc12_gg<- ggplot(data=proc12, aes(x=proc12$Distance, y=proc12$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=proc12$Distance, y=fpoints_proc12), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(proc12$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.23, linetype="dashed") +
  geom_hline(yintercept = 0.2)




############################
#Fresh market

fm1 <- fm[fm$Chr==1,]
fm2 <- fm[fm$Chr==2,]
fm3 <- fm[fm$Chr==3,]
fm4 <- fm[fm$Chr==4,]
fm5 <- fm[fm$Chr==5,]
fm6 <- fm[fm$Chr==6,]
fm7 <- fm[fm$Chr==7,]
fm8 <- fm[fm$Chr==8,]
fm9 <- fm[fm$Chr==9,]
fm10 <- fm[fm$Chr==10,]
fm11 <- fm[fm$Chr==11,]
fm12 <- fm[fm$Chr==12,]

#1
distance <- fm1$Distance
LD_data <- fm1$LD
n <- nrow(fm1)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm1 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm1_gg<- ggplot(data=fm1, aes(x=fm1$Distance, y=fm1$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm1$Distance, y=fpoints_fm1), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm1$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#2

distance <- fm2$Distance
LD_data <- fm2$LD
n <- nrow(fm2)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm2 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm2_gg<- ggplot(data=fm2, aes(x=fm2$Distance, y=fm2$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm2$Distance, y=fpoints_fm2), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm2$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#3
distance <- fm3$Distance
LD_data <- fm3$LD
n <- nrow(fm3)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm3 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm3_gg<- ggplot(data=fm3, aes(x=fm3$Distance, y=fm3$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm3$Distance, y=fpoints_fm3), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm3$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#4
distance <- fm4$Distance
LD_data <- fm4$LD
n <- nrow(fm4)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm4 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm4_gg<- ggplot(data=fm4, aes(x=fm4$Distance, y=fm4$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm4$Distance, y=fpoints_fm4), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm4$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#5
distance <- fm5$Distance
LD_data <- fm5$LD
n <- nrow(fm5)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm5 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm5_gg<- ggplot(data=fm5, aes(x=fm5$Distance, y=fm5$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm5$Distance, y=fpoints_fm5), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm5$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#6
distance <- fm6$Distance
LD_data <- fm6$LD
n <- nrow(fm6)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm6 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm6_gg<- ggplot(data=fm6, aes(x=fm6$Distance, y=fm6$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm6$Distance, y=fpoints_fm6), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm6$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#7
distance <- fm7$Distance
LD_data <- fm7$LD
n <- nrow(fm7)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm7 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm7_gg<- ggplot(data=fm7, aes(x=fm7$Distance, y=fm7$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm7$Distance, y=fpoints_fm7), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm7$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#8


distance <- fm8$Distance
LD_data <- fm8$LD
n <- nrow(fm8)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm8 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm8_gg<- ggplot(data=fm8, aes(x=fm8$Distance, y=fm8$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm8$Distance, y=fpoints_fm8), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm8$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)



#9
distance <- fm9$Distance
LD_data <- fm9$LD
n <- nrow(fm9)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm9 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm9_gg<- ggplot(data=fm9, aes(x=fm9$Distance, y=fm9$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm9$Distance, y=fpoints_fm9), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm9$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#10

distance <- fm10$Distance
LD_data <- fm10$LD
n <- nrow(fm10)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm10 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm10_gg<- ggplot(data=fm10, aes(x=fm10$Distance, y=fm10$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm10$Distance, y=fpoints_fm10), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm10$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#11

distance <- fm11$Distance
LD_data <- fm11$LD
n <- nrow(fm11)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm11 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm11_gg<- ggplot(data=fm11, aes(x=fm11$Distance, y=fm11$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm11$Distance, y=fpoints_fm11), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm11$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)


distance <- fm12$Distance
LD_data <- fm12$LD
n <- nrow(fm12)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_fm12 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



fm12_gg<- ggplot(data=fm12, aes(x=fm12$Distance, y=fm12$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=fm12$Distance, y=fpoints_fm12), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(fm12$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.12, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#########################################
##vintage



vin1 <- vin[vin$Chr==1,]
vin2 <- vin[vin$Chr==2,]
vin3 <- vin[vin$Chr==3,]
vin4 <- vin[vin$Chr==4,]
vin5 <- vin[vin$Chr==5,]
vin6 <- vin[vin$Chr==6,]
vin7 <- vin[vin$Chr==7,]
vin8 <- vin[vin$Chr==8,]
vin9 <- vin[vin$Chr==9,]
vin10 <- vin[vin$Chr==10,]
vin11 <- vin[vin$Chr==11,]
vin12 <- vin[vin$Chr==12,]

#1
distance <- vin1$Distance
LD_data <- vin1$LD
n <- nrow(vin1)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin1 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin1_gg<- ggplot(data=vin1, aes(x=vin1$Distance, y=vin1$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin1$Distance, y=fpoints_vin1), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin1$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#2

distance <- vin2$Distance
LD_data <- vin2$LD
n <- nrow(vin2)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin2 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin2_gg<- ggplot(data=vin2, aes(x=vin2$Distance, y=vin2$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin2$Distance, y=fpoints_vin2), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin2$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#3
distance <- vin3$Distance
LD_data <- vin3$LD
n <- nrow(vin3)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin3 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin3_gg<- ggplot(data=vin3, aes(x=vin3$Distance, y=vin3$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin3$Distance, y=fpoints_vin3), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin3$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#4
distance <- vin4$Distance
LD_data <- vin4$LD
n <- nrow(vin4)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin4 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin4_gg<- ggplot(data=vin4, aes(x=vin4$Distance, y=vin4$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin4$Distance, y=fpoints_vin4), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin4$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#5
distance <- vin5$Distance
LD_data <- vin5$LD
n <- nrow(vin5)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin5 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin5_gg<- ggplot(data=vin5, aes(x=vin5$Distance, y=vin5$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin5$Distance, y=fpoints_vin5), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin5$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#6
distance <- vin6$Distance
LD_data <- vin6$LD
n <- nrow(vin6)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin6 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin6_gg<- ggplot(data=vin6, aes(x=vin6$Distance, y=vin6$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin6$Distance, y=fpoints_vin6), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin6$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#7
distance <- vin7$Distance
LD_data <- vin7$LD
n <- nrow(vin7)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin7 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin7_gg<- ggplot(data=vin7, aes(x=vin7$Distance, y=vin7$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin7$Distance, y=fpoints_vin7), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin7$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#8


distance <- vin8$Distance
LD_data <- vin8$LD
n <- nrow(vin8)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin8 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin8_gg<- ggplot(data=vin8, aes(x=vin8$Distance, y=vin8$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin8$Distance, y=fpoints_vin8), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin8$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)



#9
distance <- vin9$Distance
LD_data <- vin9$LD
n <- nrow(vin9)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin9 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin9_gg<- ggplot(data=vin9, aes(x=vin9$Distance, y=vin9$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin9$Distance, y=fpoints_vin9), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin9$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#10

distance <- vin10$Distance
LD_data <- vin10$LD
n <- nrow(vin10)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin10 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin10_gg<- ggplot(data=vin10, aes(x=vin10$Distance, y=vin10$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin10$Distance, y=fpoints_vin10), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin10$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)

#11

distance <- vin11$Distance
LD_data <- vin11$LD
n <- nrow(vin11)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin11 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin11_gg<- ggplot(data=vin11, aes(x=vin11$Distance, y=vin11$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin11$Distance, y=fpoints_vin11), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin11$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)


distance <- vin12$Distance
LD_data <- vin12$LD
n <- nrow(vin12)
HW_st<-c(C=0.3)
temp_nlr<-nls(LD_data ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))
              *(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
              start=HW_st, control=nls.control(maxiter=100))
nlr_summary<-summary(temp_nlr)
new_rho<-nlr_summary$parameters[1]
fpoints_vin12 <-((10+new_rho*distance)/((2+new_rho*distance)*(11+new_rho*distance)))*
  (1+((3+new_rho*distance)*(12+12*new_rho*distance+(new_rho*distance)^2))/(n*(2+new_rho*distance)*(11+new_rho*distance)))



vin12_gg<- ggplot(data=vin12, aes(x=vin12$Distance, y=vin12$LD)) + 
  geom_point(col="grey",size=.4) +
  geom_smooth(method = "loess", se=F, col="red", size=.5, span=.9) +
  geom_line(aes(x=vin12$Distance, y=fpoints_vin12), col="blue", size=.5) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.title=element_blank(), legend.position = "none", 
        axis.text.y = element_text(size=6, color = "black", angle = 90, margin=margin(1,5,1,1,"pt"), hjust=0.5), 
        axis.text.x = element_text(size=7, colour = "black", margin=margin(5,1,1,1, "pt")), 
        axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
        axis.ticks.length=unit(1.5, "mm"),
        plot.margin=unit(c(2,5,2,2),"mm"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_rect(colour = "black", fill=NA, size=.5))+
  labs(list(x="Genetic Distance",y="r-squared", title=paste("Chromosome", unique(vin12$Chr), sep=" "))) +
  scale_y_continuous(limits = c(0, 1.00), expand = c(0, 0), breaks = round(seq(0,1,0.1),1)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0), breaks=seq(0,50,5)) +
  geom_hline(yintercept = 0.11, linetype="dashed") +
  geom_hline(yintercept = 0.2)


#compiling images 36 images into 3 images of 12
library(grid)
library(gridExtra)
library("gridExtra", lib.loc="~/R/win-library/3.2")
library("grid", lib.loc="C:/Program Files/R/R-3.2.3/library")

g1 <- arrangeGrob(proc1_gg, fm1_gg, vin1_gg,proc2_gg, fm2_gg, vin2_gg, proc3_gg, fm3_gg, 
                  vin3_gg, proc4_gg, fm4_gg, vin4_gg, ncol=3)
ggsave(file="chr1_4.pdf", g1, width = 8, height = 12, units=c("in"))


g2 <- arrangeGrob(proc5_gg, fm5_gg, vin5_gg,proc6_gg, fm6_gg, vin6_gg, proc7_gg, fm7_gg, vin7_gg, proc8_gg, fm8_gg, vin8_gg, ncol=3)
ggsave(file="chr5_8.pdf", g2, width = 8, height = 12, units=c("in"))

g3 <- arrangeGrob(proc9_gg, fm9_gg, vin9_gg,proc10_gg, fm10_gg, vin10_gg, proc11_gg, fm11_gg, vin11_gg, proc12_gg, fm12_gg, vin12_gg, ncol=3)
ggsave(file="chr9_12.pdf", g3, width = 8, height = 12, units=c("in"))





