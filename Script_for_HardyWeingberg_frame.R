Rscript --vanilla --default-packages=utils
args <- commandArgs()
path_file=args[7] ######## Path where the files are stored
file_name=args[9] ######## This file contains the chromosome number of specific sub population
file_list=list.files(path=path_file,pattern = "\\.hwe") ##### list of all files ending with ".hwe" extension
filenames1=file1[grep(pattern=file_name,file_list)]     ###### grepping only specific chromosome
col_names=do.call(rbind, strsplit(filenames1, '\\_'))[,7] ###### getting the specific population name to used as colname later
Pval_cols=paste0("P_",gsub(".","_",col_names,fixed=T))

split_data_hwe <- function(x){
        l2=read.table(x,header=T,stringsAsFactors=F)
        s<- strsplit(as.character(l2$SNP), ';')
        tt=data.frame(CHR = rep(l2$CHR, sapply(s, length)), SNP=unlist(s),P= rep(l2$P, sapply(s, length)))
        return(tt)
}
datalist1 = lapply(filenames1, function(x){split_data_hwe(x)}) #### Creating a list of dataframes
r1=Reduce(function(x,y) {merge(x,y,by=c("CHR","SNP"))}, datalist1) ##### Merging the list of dataframe 
colnames(r1)[3:(length(filenames1)+2)]= Pval_cols #### pasting the  specific colnames for each population
r1[((length(filenames1)+2)+1):((length(filenames1)*2)+2)][r1[3:(length(filenames1)+2)]>0.000001]=0  #### Extending the dataframe with zero's for hwe values less than "0.000001" 
r1[((length(filenames1)+2)+1):((length(filenames1)*2)+2)][r1[3:(length(filenames1)+2)]<=0.000001]=1 ###Extending the dataframe with ones' for hwe values more than "0.000001" 
colnames(r1)[((length(filenames1)+2)+1):((length(filenames1)*2)+2)]= Pval_cols
r1$TRUTH="F"
r1$TRUTH[rowSums(r1[((length(filenames1)+2)+1):((length(filenames1)*2)+2)])==0]="T" #### if the rowsums aren't summing to zero that means the snp is failed in specific population
write.table(r1,file=paste(file_name, "_matrix_HWE.txt", sep=""), sep="\t", row.names=F, quote=F)
