library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
library(lattice)

###Import the full dataset PC information and import the sample_ID file with market class
setwd("C:/Users/mdzievit.IASTATE/Dropbox/Classes/EEOB_546X/PEB_Power_Final_Project/Structure_Analysis/Admixture_Results/")
subpop2 <- read.table("tomato_proc.2.Q")
subpop3 <- read.table("tomato_proc.3.Q")
subpop4 <- read.table("tomato_proc.4.Q")
sampID_proc <- read.table("tomato_sampID_proc_info.txt", sep = "\t",header = TRUE)

##Figure out max ancestry before joining files
colnames(subpop2)
colnames(subpop2) <- c("Pop1","Pop2")
colnames(subpop2)

subpop2$max <- apply(subpop2,1,max)
subpop2$Pop <- colnames(subpop2)[apply(subpop2,1,which.max)]

colnames(subpop3)
colnames(subpop3) <- c("Pop1","Pop2","Pop3")
colnames(subpop3)

subpop3$max <- apply(subpop3,1,max)
subpop3$Pop <- colnames(subpop3)[apply(subpop3,1,which.max)]

colnames(subpop4)
colnames(subpop4) <- c("Pop1","Pop2","Pop3","Pop4")
colnames(subpop4)

subpop4$max <- apply(subpop4,1,max)
subpop4$Pop <- colnames(subpop4)[apply(subpop4,1,which.max)]

##Need to import these PCA information
proc_PC <- read.table("Proc_PC.txt", sep = "\t")

##Combine files
subpop2.info <- cbind(sampID_proc,subpop2)
subpop3.info <- cbind(sampID_proc,subpop3)
subpop4.info <- cbind(sampID_proc,subpop4)
rm(subpop2,subpop3,subpop4)

##Determining if there is admixture and classifying it accordingly
subpop2.info <- subpop2.info %>% mutate(NewPop = ifelse(max >= .51,Pop,"Admix"))
subpop3.info <- subpop3.info %>% mutate(NewPop = ifelse(max >= .51,Pop,"Admix"))
subpop4.info <- subpop4.info %>% mutate(NewPop = ifelse(max >= .51,Pop,"Admix"))

##Creating 5 and 7 subpop PC plot to compare with original
full_PC.subpop2 <- subpop2.info %>%
  select(SolCAP.T.number,NewPop) %>%
  left_join(proc_PC,by = "SolCAP.T.number")

colors4 <- c("Blue","Red","Green")
shapes4 <- c(1,2,3)

colors5 <- c("Gold","Blue","Red")
shapes5 <- c(0,1,2)

origPC_Plot <- ggplot(aes(x = PC1, y = PC2),data = full_PC.subpop2) + 
  geom_point(aes(colour = Market.Class, shape = Market.Class),size = 3) +
  scale_color_manual(values = colors4) +
  scale_shape_manual(values = shapes4) + 
  theme_bw() +
  theme(aspect.ratio = .75,legend.position = "bottom",plot.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(title = "Original processing figure, subpopulations assigned by breeders")

structure2_PC_Plot <- ggplot(aes(x = PC1, y = PC2),data = full_PC.subpop2) + 
  geom_point(aes(colour = NewPop, shape = NewPop),size = 3) +
  scale_color_manual(values = colors5) +
  scale_shape_manual(values = shapes5) + 
  theme_bw() +
  theme(aspect.ratio = .75,legend.position = "bottom",plot.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(title = "Subpopulation colored according to Admixture (K=2)")

pdf("Proc_PC_Comp-Subpop2.pdf",onefile = TRUE, paper = 'A4r', width = 11, height = 8) 
grid.arrange(origPC_Plot,structure2_PC_Plot,ncol = 2)
dev.off()

##############################
#subpop3 plot
full_PC.subpop3 <- subpop3.info %>%
  select(SolCAP.T.number,NewPop) %>%
  left_join(proc_PC,by = "SolCAP.T.number")

colors6 <- c("Gold","Blue","Green","Red")
shapes6 <- c(0,1,3,2)

structure3_PC_Plot <- ggplot(aes(x = PC1, y = PC2),data = full_PC.subpop3) + 
  geom_point(aes(colour = NewPop, shape = NewPop),size = 3) +
  scale_color_manual(values = colors6) +
  scale_shape_manual(values = shapes6) + 
  theme_bw() +
  theme(aspect.ratio = .75,legend.position = "bottom",plot.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(title = "Subpopulation colored according to Admixture (K=3)")

pdf("Proc_PC_Comp-Subpop3.pdf",onefile = TRUE, paper = 'A4r', width = 11, height = 8) 
grid.arrange(origPC_Plot,structure3_PC_Plot,ncol = 2)
dev.off()

##############################
#subpop4 plot
full_PC.subpop4 <- subpop4.info %>%
  select(SolCAP.T.number,NewPop) %>%
  left_join(proc_PC,by = "SolCAP.T.number")

colors7 <- c("Gold","Red","Green","Blue","Violet")
shapes7 <- c(0,2,3,1,4)

structure4_PC_Plot <- ggplot(aes(x = PC1, y = PC2),data = full_PC.subpop4) + 
  geom_point(aes(colour = NewPop, shape = NewPop),size = 3) +
  scale_color_manual(values = colors7) +
  scale_shape_manual(values = shapes7) + 
  theme_bw() +
  theme(aspect.ratio = .75,legend.position = "bottom",plot.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(title = "Subpopulation colored according to Admixture (K=4)")

pdf("Proc_PC_Comp-Subpop4.pdf",onefile = TRUE, paper = 'A4r', width = 11, height = 8) 
grid.arrange(origPC_Plot,structure4_PC_Plot,ncol = 2)
dev.off()


pdf("All4_Proc.pdf",onefile = TRUE, paper = 'A4r', width = 11, height = 8) 
grid.arrange(origPC_Plot,structure2_PC_Plot,structure3_PC_Plot,structure4_PC_Plot, ncol = 2)
dev.off()
