#! Rscript --vanilla --default-packages=utils
args <- commandArgs()
bim_file=args[7]       #### Plink binary bim file
fam_file=args[9]       #### Plink binary family file
bed_file=args[11]      #### Plink binary bed file
gds_file=args[13]       ### Name for gds file being created
kin0_file=args[15]      #### kinship matrix created by KING software
kin_file=args[17]       #### kinship matrix created by KING software
library(GWASTools, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("gdsfmt", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("SNPRelate", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("SeqVarTools", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("GENESIS", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
snpgdsBED2GDS(bed.fn = bed_file, bim.fn = bim_file, fam.fn =bed_file, out.gdsfn = gds_file) #### Convert the plink binary file into GDS format used by GENESIS
file.kin0=kin0_file
file.kin=kin_file
geno <- GdsGenotypeReader(filename = gds_file)
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)
###Run PC analysis
Kingmat <- king2mat(file.kin0=file.kin0,file.kin=file.kin,type="kinship",iids = iids)
mypcair <- pcair(genoData = genoData,kinMat = Kingmat,divMat = Kingmat)
mypcrel <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:2],training.set = mypcair$unrels)
### make sure the phenoytpe data and family data are in same order
fam=read.delim(fam_file,header=F)
pheno=read.delim("analysis_phenotypes.txt",header=T)
k1=merge(fam,pheno,by.x=c("V2"),by.y=c("SAMPLE"))
pheno.tmp = k1[match(fam$V2,k1$V2),]
#### Getting the phenotype and converting the negative values into NA's
pheno=as.vector(pheno.tmp$AFFECTION_STATUS)
pheno=pheno-1
pheno[pheno==-1]=NA
scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = mypcrel$sample.id,pc1 = mypcair$vectors[,1],pc2 = mypcair$vectors[,2],pc3 = mypcair$vectors[,3],pc4=mypcair$vectors[,4],pc5 = mypcair$vectors[,5],pc6 = mypcair$vectors[,6],pc7 = mypcair$vectors[,7],pc8 = mypcair$vectors[,8],pc9 = mypcair$vectors[,9],pc10 = mypcair$vectors[,10],pc11 = mypcair$vectors[,11],pc12 = mypcair$vectors[,12],pc13 = mypcair$vectors[,13],pc14=mypcair$vectors[,14],pc15 = mypcair$vectors[,15],pc16 = mypcair$vectors[,16],pc17 = mypcair$vectors[,17],pc18 = mypcair$vectors[,18],pc19 = mypcair$vectors[,19],pc20 = mypcair$vectors[,20],pheno = pheno))
covMatList <- list("Kin" = pcrelateMakeGRM(mypcrel))
##### Determining which PC's are signifcant apart from the first 3 PC's for detecting variation in the dataset
pc.list <- c("pc1", "pc2", "pc3")
for (i in 4:20) {
  pc.str <- paste0("pc", i)
  nullmod <- fitNullMM(scanData = scanAnnot,outcome = "pheno", covars = pc.str,covMatList = covMatList,family=binomial(link = "logit"))
  p.val <- nullmod$fixef[2, "pval"]
  if (p.val < 0.05) {
    pc.list <- c(pc.list, pc.str)
  }
}
#Run the final model
nullmod <- fitNullMM(scanData = scanAnnot,outcome = "pheno", covars = pc.list,covMatList = covMatList,family=binomial(link = "logit"))
save(nullmod, file = "TOPMED_GENESIS_ASTHMA_SETUP.Rdata")
