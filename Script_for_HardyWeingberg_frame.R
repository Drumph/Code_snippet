Rscript --vanilla --default-packages=utils
args <- commandArgs()
path_file=args[7] ######## Path where the files are stored
file_name=args[9] ######## This file contains the chromosome number of specific sub population
file_list=list.files(path=path_file,pattern = "\\.hwe")
filenames1=file1[grep(pattern=file_name,file_list)]
split_data_hwe <- function(x){
        l2=read.table(x,header=T,stringsAsFactors=F)
        s<- strsplit(as.character(l2$SNP), ';')
        tt=data.frame(CHR = rep(l2$CHR, sapply(s, length)), SNP=unlist(s),P= rep(l2$P, sapply(s, length)))
        return(tt)
}
datalist1 = lapply(filenames1, function(x){split_data_hwe(x)})
r1=Reduce(function(x,y) {merge(x,y,by=c("CHR","SNP"))}, datalist1)
colnames(r1)[3:28]= Pval_cols
r1[29:54]=r1[3:28]>0.000001
r1[29:54][r1[3:28]>0.000001]=0
r1[29:54][r1[3:28]<=0.000001]=1
colnames(r1)[29:54]= Pval_cols
r1$TRUTH="F"
r1$TRUTH[rowSums(r1[29:54])==0]="T"
write.table(r1,file=paste(file_name, "_matrix_HWE.txt", sep=""), sep="\t", row.names=F, quote=F)
