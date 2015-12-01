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

precision_fusion_table=matrix(0,ncol=matrix_size,nrow=matrix_size)
colnames(precision_fusion_table)=matrix_table[,1]
rownames(precision_fusion_table)=matrix_table[,1]
recall_fusion_table=matrix(0,ncol=matrix_size,nrow=matrix_size)
colnames(recall_fusion_table)=matrix_table[,1]
rownames(recall_fusion_table)=matrix_table[,1]
precision_perfect_table=matrix(0,ncol=matrix_size,nrow=matrix_size)
colnames(precision_perfect_table)=matrix_table[,1]
rownames(precision_perfect_table)=matrix_table[,1]
recall_perfect_table=matrix(0,ncol=matrix_size,nrow=matrix_size)
colnames(recall_perfect_table)=matrix_table[,1]
rownames(recall_perfect_table)=matrix_table[,1]


for (split_num in 1:matrix_size){
    #print(split_num)
    for (span_num in 1:matrix_size){
        #print(span_num)
        precision_fusion_table[split_num,span_num]=1/(1+as.numeric(sum_table[matrix_size*(split_num-1)+span_num,FP_col])/(as.numeric(sum_table[matrix_size*(split_num-1)+span_num,PF_col])+as.numeric(sum_table[matrix_size*(split_num-1)+span_num,F10_col])))
        recall_fusion_table[split_num,span_num]=(as.numeric(sum_table[matrix_size*(split_num-1)+span_num,PF_col])+as.numeric(sum_table[matrix_size*(split_num-1)+span_num,F10_col]))/as.numeric(num_fusion.in)
        precision_perfect_table[split_num,span_num]=as.numeric(sum_table[matrix_size*(split_num-1)+span_num,PF_col])/(as.numeric(sum_table[matrix_size*(split_num-1)+span_num,PF_col])+as.numeric(sum_table[matrix_size*(split_num-1)+span_num,F10_col])+as.numeric(sum_table[matrix_size*(split_num-1)+span_num,FP_col]))
        recall_perfect_table[split_num,span_num]=as.numeric(sum_table[matrix_size*(split_num-1)+span_num,PF_col])/as.numeric(num_fusion.in)
    }
}

png(paste(file.out,"_detected_fusion_precision_split.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(precision_fusion_table,main="precision of simulated fusion called by QueryFuse",xlab="Number of Simulated Splitting Reads", ylim=c(0,1), col=seq(2,657,length.out=matrix_size),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
legend("topright", inset=c(-0.2,0), legend=rownames(precision_fusion_table), pch=15, col=seq(2,657,length.out=matrix_size), title="Span Reads")
dev.off()
png(paste(file.out,"_detected_fusion_recall_split.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(recall_fusion_table,main="recall of simulated fusion called by QueryFuse",xlab="Number of Simulated Splitting Reads", ylim=c(0,1),col=seq(2,657,length.out=matrix_size),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
legend("topright", inset=c(-0.2,0), legend=rownames(precision_fusion_table), pch=15, col=seq(2,657,length.out=matrix_size), title="Span Reads")
dev.off()
png(paste(file.out,"_perfect_fusion_precision_split.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(precision_perfect_table,main="precision of simulated fusions perfectly called by QueryFuse",xlab="Number of Simulated Splitting Reads", ylim=c(0,1),col=seq(2,657,length.out=matrix_size),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
legend("topright", inset=c(-0.2,0), legend=rownames(precision_fusion_table), pch=15, col=seq(2,657,length.out=matrix_size), title="Span Reads")
dev.off()
png(paste(file.out,"_perfect_fusion_recall_split.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(recall_perfect_table,main="recall of simulated fusions perfectly called by QueryFuse",xlab="Number of Simulated Splitting Reads", ylim=c(0,1),col=seq(2,657,length.out=matrix_size),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
legend("topright", inset=c(-0.2,0), legend=rownames(precision_fusion_table), pch=15, col=seq(2,657,length.out=matrix_size), title="Span Reads")
dev.off()


png(paste(file.out,"_detected_fusion_precision_span.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(precision_fusion_table),main="precision of simulated fusion called by QueryFuse",xlab="Number of Simulated Spaning Reads", ylim=c(0,1), col=seq(2,657,length.out=matrix_size),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
legend("topright", inset=c(-0.2,0), legend=rownames(precision_fusion_table), pch=15, col=seq(2,657,length.out=matrix_size), title="Split Reads")
dev.off()
png(paste(file.out,"_detected_fusion_recall_span.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(recall_fusion_table),main="recall of simulated fusion called by QueryFuse",xlab="Number of Simulated Spaning Reads", ylim=c(0,1),col=seq(2,657,length.out=matrix_size),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
legend("topright", inset=c(-0.2,0), legend=rownames(precision_fusion_table), pch=15, col=seq(2,657,length.out=matrix_size), title="Split Reads")
dev.off()
png(paste(file.out,"_perfect_fusion_precision_span.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(precision_perfect_table),main="precision of simulated fusions perfectly called by QueryFuse",xlab="Number of Simulated Spaning Reads", ylim=c(0,1),col=seq(2,657,length.out=matrix_size),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
legend("topright", inset=c(-0.2,0), legend=rownames(precision_fusion_table), pch=15, col=seq(2,657,length.out=matrix_size), title="Split Reads")
dev.off()
png(paste(file.out,"_perfect_fusion_recall_span.png",sep=""),width=750,height=750)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(recall_perfect_table),main="recall of simulated fusions perfectly called by QueryFuse",xlab="Number of Simulated Spaning Reads", ylim=c(0,1),col=seq(2,657,length.out=matrix_size),beside=TRUE, cex.axis=0.01)
axis(2, at = seq(0, 1, by = 0.05),tck=1)
legend("topright", inset=c(-0.2,0), legend=rownames(precision_fusion_table), pch=15, col=seq(2,657,length.out=matrix_size), title="Split Reads")
dev.off()

#output the file
write.table(precision_fusion_table, file=paste(file.out,"_detected_fusion_precision.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(recall_fusion_table, file=paste(file.out,"_detected_fusion_recall.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(precision_perfect_table, file=paste(file.out,"_perfect_fusion_precision.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
write.table(recall_perfect_table, file=paste(file.out,"_perfect_fusion_recall.txt",sep=""), sep = "\t", row.names = TRUE , col.names = TRUE, quote=FALSE)
