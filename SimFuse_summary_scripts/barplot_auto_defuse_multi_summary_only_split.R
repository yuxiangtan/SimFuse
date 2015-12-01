#/share/apps/R/R-2.12.2_gnu412_x86_64_vanilla/bin/Rscript barplot_auto_QF_multi_result_gathering.R file.in= file.matrix= file.out= num_fusion.in=
#The input should be the count matrix
#count the precision (TP/all_found) and recall (TP/simulated) for fusions [this one the key useful one], splits, span and query
#output barplot with colors.

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
if (!exists("file.in")) {
	stop("\n\nUsage: file.in is not exist, please check the path. \n\n")
}

if (!exists("file.matrix")) {
	stop("\n\nUsage: file.matrix is not exist, please check the path. \n\n")
}

#read in the file
sum_table<-as.matrix(read.delim2(file=`file.in`,header=TRUE,stringsAsFactors=FALSE,blank.lines.skip=FALSE,fill=TRUE))
matrix_table<-read.delim2(file=`file.matrix`,header=FALSE,stringsAsFactors=FALSE,blank.lines.skip=FALSE,fill=TRUE)

PF_col=which(colnames(sum_table)=="NUM.of.perfectly.found")
F10_col=which(colnames(sum_table)=="NUM.of.fusion.found.within.10bp.of.the.simulated.breakpoint")
FP_col=which(colnames(sum_table)=="NUM.of.False.Positive.found")
matrix_size=dim(matrix_table)[1]
rowname_matrix=paste("_",matrix_table[,1],"_",matrix_table[,2],sep="")

precision_fusion_table=matrix(0,ncol=1,nrow=matrix_size)
colnames(precision_fusion_table)="precision_fusion"
rownames(precision_fusion_table)=rowname_matrix
recall_fusion_table=matrix(0,ncol=1,nrow=matrix_size)
colnames(recall_fusion_table)="recall_fusion"
rownames(recall_fusion_table)=rowname_matrix
precision_perfect_table=matrix(0,ncol=1,nrow=matrix_size)
colnames(precision_perfect_table)="precision_perfect"
rownames(precision_perfect_table)=rowname_matrix
recall_perfect_table=matrix(0,ncol=1,nrow=matrix_size)
colnames(recall_perfect_table)="recall_perfect"
rownames(recall_perfect_table)=rowname_matrix


for (row_num in 1:matrix_size){
    precision_fusion_table[row_num]=1/(1+as.numeric(sum_table[row_num,FP_col])/(as.numeric(sum_table[row_num,PF_col])+as.numeric(sum_table[row_num,F10_col])))
    recall_fusion_table[row_num]=(as.numeric(sum_table[row_num,PF_col])+as.numeric(sum_table[row_num,F10_col]))/as.numeric(num_fusion.in)
    precision_perfect_table[row_num]=as.numeric(sum_table[row_num,PF_col])/(as.numeric(sum_table[row_num,PF_col])+as.numeric(sum_table[row_num,F10_col])+as.numeric(sum_table[row_num,FP_col]))
    recall_perfect_table[row_num]=as.numeric(sum_table[row_num,PF_col])/as.numeric(num_fusion.in)
}

png(paste(file.out,"_detected_fusion_precision_split.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(precision_fusion_table,main="precision of simulated fusion called",xlab="Number of Simulated Split and Span Reads", ylim=c(0,1), beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
dev.off()
png(paste(file.out,"_detected_fusion_recall_split.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(recall_fusion_table,main="recall of simulated fusion called",xlab="Number of Simulated Split and Span Reads", ylim=c(0,1),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
dev.off()
png(paste(file.out,"_perfect_fusion_precision_split.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(precision_perfect_table,main="precision of simulated fusions perfectly called",xlab="Number of Simulated Split and Span Reads", ylim=c(0,1),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
dev.off()
png(paste(file.out,"_perfect_fusion_recall_split.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(recall_perfect_table,main="recall of simulated fusions perfectly called",xlab="Number of Simulated Split and Span Reads", ylim=c(0,1),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
dev.off()



#output the file
write.table(precision_fusion_table, file=paste(file.out,"_detected_fusion_precision.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(recall_fusion_table, file=paste(file.out,"_detected_fusion_recall.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(precision_perfect_table, file=paste(file.out,"_perfect_fusion_precision.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(recall_perfect_table, file=paste(file.out,"_perfect_fusion_recall.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
