###Workflow for Generating SNP Data for  minor allele frequency (MAF), PCA and rarefaction analyses


 
1. Start on the 4th row of the TableS2 data to get rid of the column titles and table header
> tail -n +4 Table_S2.txt > T1.txt

2. Filter data based on the data being "POLY" AND "10" FYI I changed "< 10%" to "10" in the raw data file (Table_S2.txt) prior to performing Step 1 above.
> $ awk '$3 ~ /POLY/&& $2 ~ /10/' T1.txt > T2.txt

3. Calling a "wc" command on the file T2.txt shows that there are 7,375 rows in T2.txt file
> wc T2. txt

4. Remove all rows in Column 5 (Position) of T2.txt file whose value is "Unknown" since the paper specifically mentions that SNP's with "physical positions" were filtered for the aforementioned ("10" and "POLY") characteristics.
> $ awk -F'\t' '$5!="Unknown"' T2.txt > T3.txt

5. Calling a "wc" command on the file T3.txt shows that there are 7,323 rows in the T3.txt file which is consistent with the number of rows in the paper. 
6. Removed the 13 markers with inconsistent chromosome assignments between "Physical" (column 4 of T3.txt) and "EXPIM 2012" (column 6 of T3.txt) maps.
> $ awk '{ if($4 = = $6 || $6 = = "Unknown") {print}  }'  T3.txt > T4.txt

7. Calling a "wc" command on the file T4.txt revealed that the 13 markers were removed and the no. of rows was now 7,310 as found in the research article.
> wc T4.txt

** T4.txt SHOULD BE USED FOR MINOR ALLELE FREQUENCY (MAF) AND RAREFACTION ANALYSES

**NEXT: FILTER DATA IN T4.TXT BASED ON "< 10%" DATA THAT ARE GENETICALLY MAPPED IN THE Moneymaker x LA0121 F2 population of 184 plants**




    
  









