#/share/apps/R/R-2.12.2_gnu412_x86_64_vanilla/bin/Rscript barplot_auto_QF_multi_result_gathering.R file.names="list of files to be merged" file.matrix= file.out="out name index"

#check arguments
for (e in commandArgs()) {
        ta = strsplit(e,"=",fixed=TRUE)
        if(! is.na(ta[[1]][2])) {
                temp = ta[[1]][2]
                if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
                temp = as.integer(temp)
                }
        if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
                temp = as.numeric(temp)
                }
        assign(ta[[1]][1],temp)
        } else {
        assign(ta[[1]][1],TRUE)
        }
}

#check whether files exist
if (!exists("file.names")) {
	stop("\n\nUsage: file.names is not exist, please check the path. \n\n")
}

file_list=as.matrix(read.delim(file.names,header=FALSE))
matrix_merge <- as.matrix(read.delim(file_list[1],header=TRUE,stringsAsFactors=FALSE,blank.lines.skip=FALSE,fill=TRUE))
file_num=dim(file_list)[1]

PF_col=which(colnames(matrix_merge)=="NUM.of.perfectly.found")
F10_col=which(colnames(matrix_merge)=="NUM.of.fusion.found.within.10bp.of.the.simulated.breakpoint")
FP_col=which(colnames(matrix_merge)=="NUM.of.False.Positive.found")

matrix_table<-read.delim2(file=`file.matrix`,header=FALSE,stringsAsFactors=FALSE,blank.lines.skip=FALSE,fill=TRUE)
matrix_size=dim(matrix_table)[1]
rowname_matrix=paste("_",matrix_table[,1],"_",matrix_table[,2],sep="")

precision_fusion_table=matrix(0,ncol=3,nrow=matrix_size)
colnames(precision_fusion_table)=c("precision_fusion","Standard Deviation","95%CI")
rownames(precision_fusion_table)=rowname_matrix
recall_fusion_table=matrix(0,ncol=3,nrow=matrix_size)
colnames(recall_fusion_table)=c("recall_fusion","Standard Deviation","95%CI")
rownames(recall_fusion_table)=rowname_matrix
precision_perfect_table=matrix(0,ncol=3,nrow=matrix_size)
colnames(precision_perfect_table)=c("precision_perfect","Standard Deviation","95%CI")
rownames(precision_perfect_table)=rowname_matrix
recall_perfect_table=matrix(0,ncol=3,nrow=matrix_size)
colnames(recall_perfect_table)=c("recall_perfect","Standard Deviation","95%CI")
rownames(recall_perfect_table)=rowname_matrix

#build the 3D array to hold all the data
matrix_merge_array=array(0,dim=c(file_num,dim(matrix_merge)))
for (i in 1:file_num){
    matrix_merge_array[i,,]<-as.matrix(read.delim(file_list[i],header=TRUE,stringsAsFactors=FALSE,blank.lines.skip=FALSE,fill=TRUE))
}

for (row_num in 1:matrix_size){
    precision_fusion_table[row_num,1]=mean(1/(1+as.numeric(matrix_merge_array[,row_num,FP_col])/(as.numeric(matrix_merge_array[,row_num,PF_col])+as.numeric(matrix_merge_array[,row_num,F10_col]))),na.rm = TRUE)
    recall_fusion_table[row_num,1]=mean((as.numeric(matrix_merge_array[,row_num,PF_col])+as.numeric(matrix_merge_array[,row_num,F10_col]))/as.numeric(num_fusion.in),na.rm = TRUE)
    precision_perfect_table[row_num,1]=mean(as.numeric(matrix_merge_array[,row_num,PF_col])/(as.numeric(matrix_merge_array[,row_num,PF_col])+as.numeric(matrix_merge_array[,row_num,F10_col])+as.numeric(matrix_merge_array[,row_num,FP_col])),na.rm = TRUE)
    recall_perfect_table[row_num,1]=mean(as.numeric(matrix_merge_array[,row_num,PF_col])/as.numeric(num_fusion.in),na.rm = TRUE)
    
    precision_fusion_table[row_num,2]=sd(1/(1+as.numeric(matrix_merge_array[,row_num,FP_col])/(as.numeric(matrix_merge_array[,row_num,PF_col])+as.numeric(matrix_merge_array[,row_num,F10_col]))),na.rm = TRUE)
    recall_fusion_table[row_num,2]=sd((as.numeric(matrix_merge_array[,row_num,PF_col])+as.numeric(matrix_merge_array[,row_num,F10_col]))/as.numeric(num_fusion.in),na.rm = TRUE)
    precision_perfect_table[row_num,2]=sd(as.numeric(matrix_merge_array[,row_num,PF_col])/(as.numeric(matrix_merge_array[,row_num,PF_col])+as.numeric(matrix_merge_array[,row_num,F10_col])+as.numeric(matrix_merge_array[,row_num,FP_col])),na.rm = TRUE)
    recall_perfect_table[row_num,2]=sd(as.numeric(matrix_merge_array[,row_num,PF_col])/as.numeric(num_fusion.in),na.rm = TRUE)
    
    precision_num=sum(!is.na(1/(1+as.numeric(matrix_merge_array[,row_num,FP_col])/(as.numeric(matrix_merge_array[,row_num,PF_col])+as.numeric(matrix_merge_array[,row_num,F10_col])))))
    precision_perf_num=sum(!is.na(as.numeric(matrix_merge_array[,row_num,PF_col])/(as.numeric(matrix_merge_array[,row_num,PF_col])+as.numeric(matrix_merge_array[,row_num,F10_col])+as.numeric(matrix_merge_array[,row_num,FP_col]))))
    precision_fusion_table[row_num,3]=qnorm(0.975)*precision_fusion_table[row_num,2]/sqrt(precision_num)
    recall_fusion_table[row_num,3]=qnorm(0.975)*recall_fusion_table[row_num,2]/sqrt(file_num)
    precision_perfect_table[row_num,3]=qnorm(0.975)*precision_perfect_table[row_num,2]/sqrt(precision_perf_num)
    recall_perfect_table[row_num,3]=qnorm(0.975)*recall_perfect_table[row_num,2]/sqrt(file_num)
}

png(paste(file.out,"_fusion_precision.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
mp <- barplot(precision_fusion_table[,1], ylim=c(0, 1), main="precision of simulated fusion called by deFuse", xlab="Number of Simulated Split and Span Reads",beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
box()
segments(mp, precision_fusion_table[,1] - precision_fusion_table[,3], mp, precision_fusion_table[,1] + precision_fusion_table[,3], lwd=2)
segments(mp - 0.1, precision_fusion_table[,1] - precision_fusion_table[,3], mp + 0.1, precision_fusion_table[,1] - precision_fusion_table[,3], lwd=2)
segments(mp - 0.1, precision_fusion_table[,1] + precision_fusion_table[,3], mp + 0.1, precision_fusion_table[,1] + precision_fusion_table[,3], lwd=2)
dev.off()

png(paste(file.out,"_fusion_recall.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
mp <- barplot(recall_fusion_table[,1], ylim=c(0, 1), main="recall of simulated fusion called by deFuse", xlab="Number of Simulated Split and Span Reads",beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
box()
segments(mp, recall_fusion_table[,1] - recall_fusion_table[,3], mp, recall_fusion_table[,1] + recall_fusion_table[,3], lwd=2)
segments(mp - 0.1, recall_fusion_table[,1] - recall_fusion_table[,3], mp + 0.1, recall_fusion_table[,1] - recall_fusion_table[,3], lwd=2)
segments(mp - 0.1, recall_fusion_table[,1] + recall_fusion_table[,3], mp + 0.1, recall_fusion_table[,1] + recall_fusion_table[,3], lwd=2)
dev.off()

png(paste(file.out,"_perfect_fusion_precision.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
mp <- barplot(precision_perfect_table[,1], ylim=c(0, 1), main="precision of perfect match fusion called by deFuse", xlab="Number of Simulated Split and Span Reads",beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
box()
segments(mp, precision_perfect_table[,1] - precision_perfect_table[,3], mp, precision_perfect_table[,1] + precision_perfect_table[,3], lwd=2)
segments(mp - 0.1, precision_perfect_table[,1] - precision_perfect_table[,3], mp + 0.1, precision_perfect_table[,1] - precision_perfect_table[,3], lwd=2)
segments(mp - 0.1, precision_perfect_table[,1] + precision_perfect_table[,3], mp + 0.1, precision_perfect_table[,1] + precision_perfect_table[,3], lwd=2)
dev.off()

png(paste(file.out,"_perfect_fusion_recall.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
mp <- barplot(recall_perfect_table[,1], ylim=c(0, 1), main="recall of perfect match fusion called by deFuse", xlab="Number of Simulated Split and Span Reads",beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
box()
segments(mp, recall_perfect_table[,1] - recall_perfect_table[,3], mp, recall_perfect_table[,1] + recall_perfect_table[,3], lwd=2)
segments(mp - 0.1, recall_perfect_table[,1] - recall_perfect_table[,3], mp + 0.1, recall_perfect_table[,1] - recall_perfect_table[,3], lwd=2)
segments(mp - 0.1, recall_perfect_table[,1] + recall_perfect_table[,3], mp + 0.1, recall_perfect_table[,1] + recall_perfect_table[,3], lwd=2)
dev.off()


#output the file
write.table(precision_fusion_table, file=paste(file.out,"_detected_fusion_precision.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(recall_fusion_table, file=paste(file.out,"_detected_fusion_recall.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(precision_perfect_table, file=paste(file.out,"_perfect_fusion_precision.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(recall_perfect_table, file=paste(file.out,"_perfect_fusion_recall.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
