# Description
## Script_for_HardyWeinberg_frame.R
* Plink "hardy" option was run all TGP population (i.e AFR,AMR,EAS,EUR,SAS) vcf files to output files with Hardy Weinberg (HWE) values for each snp. These files were used as input for script.(NOTE: only founders were kept when run the plink command) <br />
* This script is used to create a dataframe containing for every snp and it's corresponding HWE values in each chromosome and a columnn describing whether the snp has passed or failed HWE , indicated by F(fail) or T(pass), if the SNP failed in any of the population it was termed to failed SNP.

## Script_genesis_null_model.R
This script was used to create a null model for cases/control phenotype in Barbados Asthma studies
