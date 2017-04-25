#Read me file for the Structure Analysis part of this project.
###This file will walk you through all files contained within this folder. Please see the readme for those folders.

This contains all the .ped format files that plink needs to to convert .bed files. There are 2 sets, the full data set and the processing snp dataset (has proc in the title)

All three files, sampleID, .map, and .ped need to be input into plink for conversion

Code for plink: 
`plink –noweb –file {file name no ext} –out {out file name no ext} –make-bed`