###Workflow for Generating SNP Data for  minor allele frequency (MAF) and rarefaction analyses


 
1. Start on the 4th row of the TableS2 data
> touch tables2_notitles.txt
> tail -n +4 Table_S2.txt > tables2_notitles.txt

2. Get rid of 'Chromosome' and 'Position' columns for EXPIM 2012 Map (columns 6 & 7) 
> touch tables3_notitles.txt
> cut -d, -f-5,8- tables2_notitles.txt > tables3_notitles.txt



3. Filter based on polymorphisms with < 10% missing data
> $ grep '< 10%\|POLY' tables2_notitles.txt > percent_poly.txt

4. Filter the data in "percent_poly.txt" so that it only contains SNPs with physical positions.

	> $ awk '$5 != "Unknown"' percent_poly.txt > poly_snps.txt





