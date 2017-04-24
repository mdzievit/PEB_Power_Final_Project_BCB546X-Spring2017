library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
library(lattice)

###Import the full dataset PC information and import the sample_ID file with market class
setwd("C:/Users/mdzievit.IASTATE/Dropbox/Classes/EEOB_546X/PEB_Power_Final_Project/Structure_Analysis/Admixture_Results/")
subpop5 <- read.table("tomato.5.Q")
subpop7 <- read.table("tomato.7.Q")
sampID_full <- read.table("tomato_sampID_info.txt", sep = "\t",header = TRUE)

##Figure out max ancestry before joining files
colnames(subpop5)
colnames(subpop5) <- c("Pop1","Pop2","Pop3","Pop4","Pop5")
colnames(subpop5)

subpop5$max <- apply(subpop5,1,max)
subpop5$Pop <- colnames(subpop5)[apply(subpop5,1,which.max)]


colnames(subpop7)
colnames(subpop7) <- c("Pop1","Pop2","Pop3","Pop4","Pop5","Pop6","Pop7")
colnames(subpop7)

subpop7$max <- apply(subpop7,1,max)
subpop7$Pop <- colnames(subpop7)[apply(subpop7,1,which.max)]

##Need to import these PCA information
full_PC <- read.table("full_PC.txt", sep = "\t")

##Combine files
subpop5.info <- cbind(sampID_full,subpop5)
subpop7.info <- cbind(sampID_full,subpop7)
rm(subpop5,subpop7)

##Determining if there is admixture and classifying it accordingly
subpop5.info <- subpop5.info %>% mutate(NewPop = ifelse(max >= .51,Pop,"Admix"))
subpop7.info <- subpop7.info %>% mutate(NewPop = ifelse(max >= .51,Pop,"Admix"))

##Creating 5 and 7 subpop PC plot to compare with original
full_PC.subpop5 <- subpop5.info %>%
  select(SolCAP.T.number,NewPop) %>%
  left_join(full_PC,by = "SolCAP.T.number")

colors <- c("violet","Blue","Gold","Red","Green","Black","Gray")
shapes <- c(0,1,4,2,3,5,22)

colors2 <- c("Orange","Blue","Black","Green","Gray","Red")
shapes2 <- c(13,1,5,3,0,2)

origPC_Plot <- ggplot(aes(x = PC1, y = PC2),data = full_PC.subpop5) + 
  geom_point(aes(colour = Market.Class, shape = Market.Class),size = 3) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) + 
  theme_bw() +
  theme(aspect.ratio = 1,legend.position = "bottom",plot.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  labs(title = "Original figure, subpopulations assigned by breeders")

structure5_PC_Plot <- ggplot(aes(x = PC1, y = PC2),data = full_PC.subpop5) + 
  geom_point(aes(colour = NewPop, shape = NewPop),size = 3) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = shapes2) + 
  theme_bw() +
  theme(aspect.ratio = 1,legend.position = "bottom",plot.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  labs(title = "Subpopulation colored according to Admixture (K=5)")
grid.arrange(origPC_Plot,structure5_PC_Plot,ncol = 2)


pdf("PC_Comp-Subpop5.pdf",onefile = TRUE, paper = 'A4r', width = 11, height = 8) 
grid.arrange(origPC_Plot,structure5_PC_Plot,ncol = 2)
dev.off()

full_PC.subpop7 <- subpop7.info %>%
  select(SolCAP.T.number,NewPop) %>%
  left_join(full_PC,by = "SolCAP.T.number")

colors <- c("violet","Blue","Gold","Red","Green","Black","Gray")
shapes <- c(0,1,4,2,3,5,22)

colors3 <- c("Orange","Red","Gray","Black","Red","Blue","Green","Red")
shapes3 <- c(13,17,0,5,2,1,3,25)

structure7_PC_Plot <- ggplot(aes(x = PC1, y = PC2),data = full_PC.subpop7) + 
  geom_point(aes(colour = NewPop, shape = NewPop),size = 3) +
  scale_color_manual(values = colors3) +
  scale_shape_manual(values = shapes3) + 
  theme_bw() +
  theme(aspect.ratio = 1,legend.position = "bottom",plot.title = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  labs(title = "Subpopulation colored according to Admixture (K=7)")
grid.arrange(origPC_Plot,structure7_PC_Plot,ncol = 2)

pdf("PC_Comp-Subpop7.pdf",onefile = TRUE, paper = 'A4r', width = 11, height = 8) 
grid.arrange(origPC_Plot,structure7_PC_Plot,ncol = 2)
dev.off()

pdf("PC_Comp_both.pdf",onefile = TRUE, paper = 'A4r', width = 11, height = 8)
grid.arrange(structure5_PC_Plot,structure7_PC_Plot,ncol = 2)
dev.off()

#####Make the structure barplot figure for subpop5
subpop5.tidy <- subpop5.info %>%
  arrange(NewPop,desc(Pop1),desc(Pop2),desc(Pop3),desc(Pop4),desc(Pop5)) %>%
  select(-Market.Class,-max,-Pop) %>%
  mutate(ind = c(1:dim(subpop5.info)[1])) %>%
  gather(SubPop,Ancestry,-ind,-SolCAP.T.number,-NewPop)

bar_color <- c("Blue","Black","Green","Gray","Red")
five_subpop <- ggplot(aes(x = ind,y = Ancestry,fill = SubPop), data = subpop5.tidy) + 
  geom_bar(stat = "identity",width = 1) +
  scale_fill_manual(values = bar_color) +
  facet_wrap(~NewPop, strip.position = "left", scales = "free_y", ncol = 1) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") + 
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.text.y = element_text(angle = 180), panel.spacing = unit(0,"lines")) + 
  coord_flip() +
  labs(subtitle = "Admixture k = 5")

#####Make the structure barplot figure for subpop7
subpop7.tidy <- subpop7.info %>%
  arrange(NewPop,desc(Pop1),desc(Pop2),desc(Pop3),desc(Pop4),desc(Pop5),
          desc(Pop6),desc(Pop7)) %>%
  select(-Market.Class,-max,-Pop) %>%
  mutate(ind = c(1:dim(subpop5.info)[1])) %>%
  gather(SubPop,Ancestry,-ind,-NewPop,-SolCAP.T.number)


subpop7.tidy.noind <- subpop7.info %>%
  arrange(NewPop,desc(Pop1),desc(Pop2),desc(Pop3),desc(Pop4),desc(Pop5),
          desc(Pop6),desc(Pop7)) %>%
  select(-Market.Class,-max,-Pop) %>%
  gather(SubPop,Ancestry,-NewPop,-SolCAP.T.number)

subpop7.ordered5 <- subpop5.tidy %>%
  select(SolCAP.T.number,ind,NewPop) %>%
  left_join(subpop7.tidy.noind,by = c("SolCAP.T.number")) %>%
  unique()



bar_color2 <- c("Red","Gray","Black","RosyBrown","Blue","Green","FireBrick4")

seven_subpop <- ggplot(aes(x = ind,y = Ancestry,fill = SubPop), data = subpop7.tidy) + 
  geom_bar(stat = "identity",width = 1) +
  facet_wrap(~NewPop, strip.position = "left", scales = "free_y", ncol = 1) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") + 
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.text.y = element_text(angle = 180), panel.spacing = unit(0,"lines")) + 
  coord_flip()

seven_subpop2 <- ggplot(aes(x = ind,y = Ancestry,fill = SubPop), data = subpop7.ordered5) + 
  geom_bar(stat = "identity",width = 1) +
  facet_wrap(~NewPop.x,scales = "free_y", ncol = 1) +
  scale_fill_manual(values = bar_color2) +
  theme_classic() +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  coord_flip() +
  labs(subtitle = "Admixture k = 7")

pdf("Comp_5_7_admixture.pdf",onefile = TRUE, paper = 'A4r', width = 11, height = 8)
grid.arrange(five_subpop, seven_subpop2,ncol = 2,
             top = "Comparison of k = 5 to k = 7 classification")
dev.off()
