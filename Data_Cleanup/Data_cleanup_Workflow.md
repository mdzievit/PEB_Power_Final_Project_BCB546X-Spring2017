##Clean the data and create two files for the later on analysis from file Table_S2.txt

###Object for this step: Keep only those SNPs that have missing rate smaller than 10% and polymorphic and remove those SNPs that has missing rate larger or equal to 10% and those that noted as non_polymorphic SNPs

###1. Insepect the data. Notice that the file Table_S2.txt was converted from original excel table which can be downloaded from the website.

**- Remove the first three header lines of the file Table_S2.txt and create Tm1.txt file to make it easier to work with, and check the number of SNPs within the file**

```
tail -n +4 Table_S2.txt > Tm1.txt
```

```
wc -l Tm1.txt
```

The above command shows in total there are **7720 SNPs** in the file. Which agrees with the number of SNPs presented in the paper.

**- check to see how many categories are there for 2nd and 3rd columns that we are interested to work with**

```
cut -f2 Tm1.txt | sort | uniq | cat
```

```
cut -f3 Tm1.txt | sort | uniq | cat
```

From the above code we could see each of those two columns have four columns. 1.The 2nd column contains information for missing rate, four categories are respectively are: **< 10%, > 20%,  100%, 10-20%**. 2.The 3rd columns contains information indicating whether the SNPs are polymorphic or not, four categories are respectively : **MONO, No Call, POLY, Undetermined**

###2. Clean the data based on last step's inspection

**- Using the grep command to filter the data and keep those SNPs that have missing rate smaller than 10% and polymorphic. And then check the number of SNPs left after filtering and make sure there are only 1 category left for 2nd and 3rd column**

```
grep -e '< 10%' Tm1.txt | grep -e 'POLY' > Tm2.txt
```

```
wc -l Tm2.txt
```

```
cut -f2 Tm2.txt | sort | uniq | cat
```

```
cut -f3 Tm2.txt | sort | uniq | cat
```

The above command shows after filtering there are **7375 SNPs** left in the **Tm2.txt** file. And cnfirmed there are only one category left for each of those two columns.

###3. Create those two files. One file with the physical map information and genotype information, the  other file contains genetic map information and genotype information.

**-Check how many columns of file Tm2.txt which will be useful for later on data extraction**

```
awk -F "\t" '{print NF; exit}' Tm2.txt
```

The above line of code shows there are in total **433 columns** of file Tm2.txt, which also agrees with 7 columns + 426 genotypes presented in the file


**-Get the physical map and genotype information and put it to file named physical_map_data.txt. Then check the number of SNPs in the file**

```
awk -F "\t" '$4 !~ /Unknown/ && $5 !~ /Unknown/' Tm2.txt| wc -l
```

From this above line of code we could see that there are **7323 SNPs** that have known physical position 

```
awk -F "\t" '$4 !~ /Unknown/ && $5 !~ /Unknown/' Tm2.txt | awk -F "\t" '$4 == $6 || $6 ~ /Unknown/' | cut -f1,4-5,8-433 > physical_map_data.txt 
```

The above line of code first use the **awk** command to remove those SNPs that have unknown physical map position and then pass it to next **awk** command which remove SNPs that have inconsistent chromosome assignment between  physical map (column 4) and EXPIM 2012 (column 6). After that the output from second **awk** command was passed to **cut** command to extract only the marker, physical map and the genotype information.

```
wc -l physical_map_data.txt
```

The **wc** command shows that there are 7310 SNPs left in the physical_map_data.txt file

**- First, filter those SNPs that have unknown genetic map position and remove those 13 SNPs that have inconsistent chromosome assignment between column 4 and column 6. After that extract columns that corresponding to genetic map and genotype information of EXPIM2012. Then put it to file named genetic_map_data_EXPIM2012.txt. The last step check the number of SNPs left after filtering**

```
awk -F "\t" '$6 !~ /Unknown/ && $7 !~ /Unknown/' Tm2.txt | awk -F "\t" '$4 == $6 || $4 ~ /Unknown/' | cut -f1,6-433 > genetic_map_data_EXPIM2012.txt
```

The above line of code first use the **awk** command to remove those SNPs that have unknown EXPIM2012 genetic map position and then pass it to next **awk** command which remove those SNPs that have inconsistent chromosome assignment between  physical map (column 4) and EXPIM 2012 (column 6). After that the output from second **awk** command was passed to **cut** command to extract only the marker, genetic map and the genotype information. 

```
wc -l genetic_map_data_EXPIM2012.txt
```


The above code shows that there are 4393 polymorphism SNPs present in the file genetic\_map\_data\_EXPIM2012.txt . And this agrees with what it presented in the paper. 

###4. create header lines and add it to the physical\_map\_data.txt and genetic\_map\_data\_EXPIM2012.txt to create another two files physical\_map\_data\_wh.txt, genetic\_map\_data\_EXPIM2012\_wh.txt.

**-The following command lines are used to create the header files for physical\_file\_header.txt and the genetic\_file\_header.txt**

```
head -n 2 Table_S2.txt | tail -n 1 | cut -f1 > header1.txt
```

```
head -n 3 Table_S2.txt | tail -n 1 | cut -f4-5 > header_physical.txt
```

```
head -n 3 Table_S2.txt | tail -n 1 | cut -f6-7 > header_genetic.txt
```

```
head -n 3 Table_S2.txt | tail -n 1 | cut -f8-433 > header_genotype.txt
```

The above code create the component of the header file

```
paste header1.txt header_physical.txt header_genotype.txt > physical_file_header.txt
```

```
paste header1.txt header_genetic.txt header_genotype.txt > genetic_file_header.txt
```

Using **paste** command to create the header file physical\_file\_header.txt and genetic\_file\_header.txt

```
cat physical_file_header.txt physical_map_data.txt > physical_map_data_wh.txt
```

```
cat genetic_file_header.txt genetic_map_data_EXPIM2012.txt > genetic_map_data_EXPIM2012_wh.txt
```

Using  **cat** command to combine the header file with the previously created the physical\_map\_data.txt and the 
genetic\_map\_data\_EXPIM2012.tx to create the physical\_map\_data\_wh.txt and the genetic\_map\_data\_EXPIM2012\_wh.txt


**The file physical\_map\_data\_wh.txt will be used for the minor allele frequency and refraction analysis. And the file named genetic\_map\_data\_EXPIM2012\_wh.txt will be used for (PCA), pairwise Fst, descriptive statistics**






