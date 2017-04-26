##README for creation of Figures 5 - 7

This file contains the code and results for re-creating images 5 to 7 in the Sim et al. (2012) paper. 

The code for calculating LD is in the .R file named 'creatingcode2.R' and uses the files Table_S1.txt and 'genetic map data EXPIM2012 wh.txt'. The code creates 4 .csv files as output. The LD calculations for three sub-populations of interest (processing, fresh market and vintage), and a .csv containing all sub-populations combined (named 'all subpops LD combined.csv')

The LD data saved in the .csv files is used to create the images. The code for image creation is in the .R file ''graphloop3.R'. The output from that code is three PDF files, each containing 12 LD decay images (4 images from each of the three sub-populations).