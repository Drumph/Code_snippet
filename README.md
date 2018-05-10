# Description

## Script_for_HardyWeinberg_frame.R
* Plink "hardy" option was run all TGP population (i.e AFR, AMR, EAS, EUR, SAS) vcf files to output files with Hardy Weinberg (HWE) values for each snp. These files post plink command were used as input for the script. (NOTE: only founders were kept when running the plink command) <br />
* This script is used to create a dataframe containing all snp and it's corresponding HWE values in each chromosome and a column describing whether the snp has passed or failed HWE, indicated by F(fail) or T(pass), if the SNP failed in any of the population it was termed to be failed.
* The purpose of this exercise was to keep tabs how many SNP's were removed and the reason for it's removal. A similar script was written for other SNP QC filtering such as missingness, mendel, frequency etc.

## Script_genesis_null_model.R
* This script was used to create a null model for cases/control phenotype in Barbados Asthma studies
* The input to this script will be plink binary data with SNP's removed that failed HWE, missingness and frequency filters and then LD pruned.
* The null model was used for association analysis of all SNP's.
