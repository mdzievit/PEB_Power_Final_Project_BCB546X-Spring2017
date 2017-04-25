#!/usr/local/bin/perl -w

##This code calculate the minor allele frequency for 7310 SNPs for each of the four subgroup. 
##And here need to be clear for those SNPs denote as -- or heterozygous we do not consider it when calculate MAF.

use strict;

my $file_dir = '/media/lixr/_media_disk_1_/Jinyu/BCB546Project/';         ##Specify the directory for the input file and output file                      
my $output_file = $file_dir.'SNP_MAF.txt';                                ##Point out the output file

my @subgroup = ("Processing", "Fresh_market", "Vintage", "Pimp");         ##Decalare an array named subgroup and within the array there are four variable corresponing to the four subgroups name.
my ($sub_s,$sub_e) = ($ARGV[0], $ARGV[1]);                                ##Decalre two input varaible and this will allow us to control which group we want to calculate from the keyboard

#####################################################################
open (OUT, '>'.$output_file) || die;	                                    ##open the output file handle
print OUT "Group\tSNP\tChromosome\tPosition\tMAF\n";                      ##print out header line to the output file

for (my $sub = $sub_s; $sub <= $sub_e; $sub++){                           ##Using a forloop to loop through each of the four subgroups
	                                        
	my $subg = $subgroup[$sub];
	my $input_file = $file_dir.$subg.'_file.txt';                           ##open output file
	next unless -e $input_file;
	open (F, $input_file) || die;
			

		while (<F>) {                                                         ##Using a while loop to loop through each line of the file
			chomp;	
			my @t = split /\t/;	                                                ##split each line of the file with tab
			next if $t[1] =~ /Chromosome/;                                      ##Skip the header line of the file
			my @allele_info = @t[0..2];                                         ##Put the SNP information which stored in the first three column in ecah line to an array named allele_info
			my @strain_allele;                                                  ##Decalare a new array variable and it will be used to store the SNP genotype information for the SNP across each of the strain
			for ( my $index =3;  $index <=$#t; $index=$index+1){                ##loopo through each of the genotype of that SNP
			my $genotype = $t[$index];
			my $gg;
			if ($genotype =~ /[ACGTN]/){                                        ##Using if else statement to check whether the genotype is missing --, and if it is missing replace the -- with NN. 
			$gg = $genotype;					
				}
			else {$gg = 'NN';}
			push @strain_allele, $gg;                                           ##store the genotype information fo a SNP to the strain_allele arrary
			}
			my $total=0;                                                        ##Declare an scalar variable named total with initial value 0, and it will be used to count the number of homozygous genotype
			my %hash;                                                           ##Declare an hash variable which will be used to store the homozygous genotype (key) and their count (value)
			foreach my $a (@strain_allele) {                                    ##loop through each of the genotype stored in strain_allele array
				next if $a =~ /N/;                                                ##skip the genotype if it is missing and matches NN
				my @snp_info = split('', $a);                                     ##split each genotype to and check whether it is an heterozygous 
				next if $snp_info[0] !~ $snp_info[1];                             ##If the genotype is heterozygous, skip it and do not count it to the total 
				$total ++;
				$hash{$a} ++;
			}
			my @sorted_allele	= sort { $hash{$a} <=> $hash{$b} } keys %hash;    ##sort the hash variable according to the key value, order the key based on their value from small to large. And store the sorted key in and array named sorted_allele
			my $maf;                                                            ##Decalre an scalar varaible named maf which will be used to store minor allele frequency for each of the snp.
			if($#sorted_allele == 1){                                           ##Using a if else statement to control the category of the homozygous genotype, and only consider those with two homozygous genotypes,, if this condition not met, the maf will be assgined as 0.
				$maf = sprintf "%.2f", $hash{$sorted_allele[0]} / $total;         ##The minor allele frequency is calculated as the number of minor homozygous allele divide by total number of homozygous allele
				}
			else {$maf = sprintf "%.2f",0 ;}
			print OUT $subg."\t";                                               ##Print out the information to the outout file.
			print OUT $_."\t" foreach (@allele_info);
			print OUT $maf."\n";			
		}
}


