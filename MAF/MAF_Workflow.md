##Work on to generate the data for Fig.3A and Fig.4A. The data files and documentation for this part are in folder named MAF

#### Object for this step: To generate the figure of MAF patterns across each chromosome for the 4-sub populations. 
 
###1.Create files that will be used to calculate MAF for each of the sub-populations. 


Before move on to  do anything with the following code. First, I have checked that the inbred order in file Table\_S1.txt correspondind to the column order of inbreds in file physical\_map\_data\_wh.txt.


####(1) Creating two initial files. 

The following command lines are used to create file MAF\_f1.txt and MAF\_f2.txt and check their content as we expected. One file will be created from physical\_map\_data\_wh.txt and be named MAF\_f1.txt contains the physical SNP marker information for the 410 inbreds accessions. And the second file will be created from Table\_S1.txt named MAF\_f2.txt contains the information that connect each inbred to their corresponding subgroup

```
cut -f1-413 physical_map_data_wh.txt  > MAF_f1.txt
```

```
awk -F "\t" '{print NF; exit}' MAF_f1.txt
```

```
tail -n +2 Table_S1.txt | cut -f1,6 | head -n 411 > MAF_f2.txt
```

```
head MAF_f2.txt
```

```
tail MAF_f2.txt
```

####(2) Determine the minor allele for each of the SNP markers based on genotypic data of 410 inbred accessions with perl code. Then add this minor allele information to the extracted physical map data in file named MAF\_f1.txt.
**The perl code does the calculation is stored in a file named as Figure\_out\_minor\_allele.pl (there is detailed annotation about the perl code in the file). The following command line are used to run the perl code. And the output file from the perl code named as Minor\_allele.txt** 

```
perl Figure_out_minor_allele.pl
```
**Then the following command are used to combine the minor allele information from file Minor\_allele.txt to the file MAF\_f1.txt and generate a new file named MAF\_f1\_n.txt**

```
cut -f4-413 MAF_f1.txt > MAF_f1_1.txt
```

```
paste Minor_allele.txt MAF_f1_1.txt > MAF_f1_n.txt
```
####(3) Create the combined file for data extraction. We will create a file named MAF\_f2\_f1\_combine.txt which contain the group information for each inbred and it will be used to extract the SNP marker information for the next step.

**First, transpose the file MAF\_f1.txt with the **awk** command and the transpose function in file transpose.awk and named the transposed file as **MAF\_f1\_transposed.txt****. In this way, Inbred ID will be in one column and this data format can be used to merge with the group information extracted from MAF\_f2.txt. After that we inspect the transposed data.

```
awk -f transpose.awk MAF_f1_n.txt  > MAF_f1_transposed.txt
```

```
awk -F "\t" '{print NF; exit}' MAF_f1_transposed.txt
```

```
wc -l MAF_f1_transposed.txt
```

**Subset two files from the transposed data MAF\_f1\_transposed.txt.** One file contain the first four lines of the file and named MAF\_f1\_transposed\_m1.txt. And the second file named MAF\_f1\_transposed_m2.txt contains all the information without the first three lines**

```
head -n 4 MAF_f1_transposed.txt > MAF_f1_transposed_m1.txt
```

```
tail -n +5 MAF_f1_transposed.txt > MAF_f1_transposed_m2.txt
```

**Subset one file from the MAF\_f2.txt.** The following code first remove the header line of file MAF_f2.txt and then it pass the information to *cut** command and extract the second column and put it in a file named MAF\_f2\_m1.txt

```
tail -n +2 MAF_f2.txt | cut -f2 > MAF_f2_m1.txt
```

**Merge file MAF\_f2\_m1.txt and file MAF\_f1\_transposed\_m2.txt with paste command to get file MAF\_f1\_transposed\_m1.txt.** And this file will be used to extract the information for each subgroup

```
paste MAF_f2_m1.txt  MAF_f1_transposed_m2.txt > MAF_f2_f1_combine.txt
```

```
awk -F "\t" '{print NF; exit}'  MAF_f2_f1_combine.txt
```

####(4) Extract the genotype for information for each subgroup from file MAF\_f1\_transposed\_m1.txt. Here I only focus on to process 4 subgroups as the MAF of those four groups are plotted in the main paper.

**First, check how many subgroup are their using the folloing command line.** And this will help us for the subseting step.

```
sort MAF_f2_m1.txt | uniq| cat
```

**The above command lines shows there are 9 subgroups of inbreds in the data. And they are Processing, Cultivated Cherry,Fresh market,Fresh market/Processing, Latin American Cultivar, Vintage, Vintage/Processing, Wild, and Wild Cherry. Since the paper says there are in total 7 inbreds subgroups. And it did not provide further information how each inbreds are assigned. Here I think the Fresh market/Processing group belong to the Fresh market group   Vintage/Processing group belong to the Vintage group. Besides that the Wild group are called Pimp group from here**

First, extract information for the **Processing subgroup** with the **awk** command. After extractting the information pass the outpu to **cut** command to remove the first column as we will not need it.

```
 awk -F "\t" '$1 ~ /Processing/&&$1 !~ /Vintage/&& $1 !~/Fresh market/'   MAF_f2_f1_combine.txt | cut -f2-7312 > Processing.txt
```

With the same way, I extrac the **Fresh market, Vintage and the Pimp** subgroup information. 

```
 awk -F "\t" '$1 ~/Fresh market/'   MAF_f2_f1_combine.txt | cut -f2-7312 > Fresh_market.txt
```

```
 awk -F "\t" '$1 ~ /Vintage/'   MAF_f2_f1_combine.txt | cut -f2-7312 > Vintage.txt
```

```
awk -F "\t" '$1 ~ /Wild/&&$1 !~ /Cherry/'   MAF_f2_f1_combine.txt | cut -f2-7312 > Pimp.txt
```

**Inspect the extracted file with the wc -l command and the awk to check the row and columns for each of the file** 

```
wc -l Fresh_market.txt Processing.txt Pimp.txt Vintage.txt
```

Here is the output from the wc -l command

    110 Fresh_market.txt
    141 Processing.txt
     16 Pimp.txt
     61 Vintage.txt
    328 total

```
awk -F "\t" '{print NF; exit}' Fresh_market.txt Processing.txt Pimp.txt Vintage.txt
```

####(5) Merge the file MAF\_f1\_transposed_m1.txt with each of the subgroup information respectively. Then transpose the file back. Those files will be used to calculate the MAF information for each of the subgroup

**The folloing command first merge the file MAF\_f1\_transposed_m1.txt with Processing.txt and then pass the output to transpose to transpose the merged file. The new file will be named as Processing\_file.txt**

```
cat MAF_f1_transposed_m1.txt Processing.txt | awk -f transpose.awk > Processing_file.txt
```
**Inspect the file Processing\_file.txt and confirmed that there are 144 columns in the file which matches there are 141 inbreds + 3 columns of (marker, chromosome and position)**

```
awk -F "\t" '{print NF; exit}' Processing_file.txt
```

**Go through similar process for the other three groups and confirmed there are respectively 113, 64 and 19 columns for file Fresh\_market\_file.txt, Vintage\_file.txt and Pimp\_file.txt**

```
cat MAF_f1_transposed_m1.txt Fresh_market.txt | awk -f transpose.awk > Fresh_market_file.txt
```

```
 awk -F "\t" '{print NF; exit}' Fresh_market_file.txt
```

```
cat MAF_f1_transposed_m1.txt Vintage.txt | awk -f transpose.awk > Vintage_file.txt
```

```
awk -F "\t" '{print NF; exit}' Vintage_file.txt
```

```
cat MAF_f1_transposed_m1.txt Pimp.txt | awk -f transpose.awk > Pimp_file.txt
```

```
awk -F "\t" '{print NF; exit}' Pimp_file.txt
```

**The above file file Processing\_file.txt, Fresh\_market\_file.txt, Vintage\_file.txt and Pimp\_file.txt will be processed in perl later on to calculate the MAF**


###2. Run the perl code to calculate MAF for each of the four sub-populations. 

####(1) The perl code is stored in file named as Calculate\_MAF\_from\_SNP_matrix.pl

**The perl code is run by the following command line**

```
perl Calculate_MAF_from_SNP_matrix.pl 0 3 &
```

####(2) The output file is named as SNP_MAF.txt and it will be used to plot out the figure3A and 4A in R program. 



###3. Plot the MAF for each of the 4 subpopulations of 12 chromosomes using ggplot2.  The R code is documented in file named MAF\_analysis.RMD within folder named MAF



###4. Check the correlation between my calculation of MAF for those 4 subgroups and those from the paper that is provided in Table\_S4. The R code is documented in file named MAF\_analysis.RMD within folder named MAF. Before load the data into R need to format the Table\_S4 in unix.

####(1) The following code is used to remove the first two rows of Table\_S4

```
tail -n +3 Table_S4.txt > Table_S4_m1.txt
```