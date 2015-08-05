#   Copyright {2015} Yuxiang Tan
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

#input file do not need to be sorted
#seperate the file into different matrix (record the number of matrix)
##when seperate, if a group has fewer than read_bound reads, merge it with next group.
#for each matrix:
#Generate N random number within range as on the left, and N random number in each other matrix as on the right.
#generate random location of each pair rows
#check the chr, distance(>1000000) for each pair to avoid overlap on same gene. then get the top (num_pair_def) pairs from the rest. (num_pair_def<=100)
#raw_count_based_random_background_exon_list_generator.R file.in="TL_03_coverage_on_HG.GRCh37.64.exons.txt" num_pair=N read_bound=M(M>100*N) file.out=(output name index, including the bed file that fastafromBed can use) frag_len=(length of fragment on each end)

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

#check whether file in is exist
if (!exists("file.in")) {
	stop("\n\nUsage: raw count matrix is not exist, please check the path. \n\n")
}

#to avoid exponential notation
options("scipen"=1000000)

#read in file
file_bed_sort<-read.delim2(file=paste(file.in,"_sort",sep=""),sep='\t',header=FALSE,stringsAsFactors=FALSE)

row_num<-dim(file_bed_sort)[1]
num_matrix<-nchar(file_bed_sort[row_num,6])

up_bound<-1
count<-1
low_row<-1
up_row<-1
#group the sorted file into subgroups
group_matrix<-matrix(ncol=3,nrow=num_matrix+1)

while(up_row<row_num){
	up_row<-(which(file_bed_sort[,6]<up_bound))[length(which(file_bed_sort[,6]<up_bound))]
	group_matrix[count,1]<-up_bound
	group_matrix[count,2]<-low_row
	group_matrix[count,3]<-up_row
	low_row<-up_row+1
	up_bound<-up_bound*10
	count=count+1
}
frag_len=as.numeric(frag_len)
num_pair=as.numeric(num_pair)
#read_bound means at least this number of genes in a expression bin, if not, merge them.
read_bound=as.numeric(read_bound)
if (read_bound<num_pair*100){
	read_bound=num_pair*100
}

#in case the read_bound it bigger than row_num            
if (read_bound>row_num){
    read_bound=row_num-1
}
num_pair_def=num_pair*2

#new added:
#reassign the subgroups:
#check row number of each group, if less than M,then merge with next group. 
#count the number of groups in all at the end and generate new matrix. Use the new matrix to generate pairs
###be careful about the para name changes of the code!!!!!

check_list=rep(0,(num_matrix+1))
fr=1

for (i in 1:(num_matrix+1)){
    ii=1
    
    if (check_list[i]==1){
        next
    }
    while ((group_matrix[i,3]-group_matrix[i,2])<read_bound){
        if (i+ii<=num_matrix+1) {
            group_matrix[i,3]=group_matrix[i+ii,3]
            group_matrix[i,1]=group_matrix[i+ii,1]
            check_list[i+ii]=1
            ii=ii+1
            if (i+ii<=num_matrix+1) {
                group_matrix[i+ii,2]=group_matrix[i,3]+1
            }
        } else {
            im=1
            while (check_list[i-im] == 1){
                im=im+1
            }
            group_matrix[i-im,3]=group_matrix[i,3]
            group_matrix[i-im,1]=group_matrix[i,1]
            if ((group_matrix[i,3]-group_matrix[i,2])<read_bound){
                check_list[i]=1
            }
                
            break
        }
    }
}
    
group_matrix=group_matrix[which(check_list==0),]
#count-1 can be replace by just count? Fix it later when to clean up the code.
count=length(which(check_list==0))+1

if (count<=2) {
    group_matrix=t(group_matrix)
}

#add the checkpoint to make sure the matrix is following that have enough gene ID for each group and each gene ID showed up only once.
checkpoint=1
checkpoint_count=0
while (checkpoint==1) {
    checkpoint=0
    if (checkpoint_count>50) {
        stop("\n\nThis data set resampled 50 times already but was not able to get a good sampling.\n\n")
    }
    #in order to protect from picking the same number for more than a time. For each group, I will just sample one time.
    row_matrix_whole<-matrix(ncol=num_pair_def*count,nrow=count-1)
    for (j in 1:(count-1)){
            row_matrix_whole[j,]<-sample(group_matrix[j,2]:group_matrix[j,3],size=num_pair_def*count,replace=FALSE)
    }
    
    gene_ID_list=list()
    #build the gene_list of all the genesID sampled. However, because there might be not enough genes to samples, I will only keep the query genes to be unique
    for (k in 1:(num_pair_def)){
        for (j in 1:(count-1)){
            anno<-strsplit(as.character(file_bed_sort[row_matrix_whole[j,k],5]),";")
            gene_id<-strsplit(anno[[1]][1]," ")[[1]][3]
            if (length(gene_ID_list[[gene_id]])==0){
                gene_ID_list[[gene_id]]=1
            } else {
                row_matrix_whole[j,k]=0
            }
        }
    }
    
    for (k in (num_pair_def+1):(num_pair_def*count)){
        for (j in 1:(count-1)){
            anno<-strsplit(as.character(file_bed_sort[row_matrix_whole[j,k],5]),";")
            gene_id<-strsplit(anno[[1]][1]," ")[[1]][3]
            if (length(gene_ID_list[[gene_id]])>0){
                row_matrix_whole[j,k]=0
            }
        }
    }
    
    row_matrix_left_temp<-matrix(ncol=num_pair_def,nrow=count-1)
    for (j in 1:(count-1)){
        row_matrix_left_temp[j,]<-row_matrix_whole[j,1:num_pair_def]
        if (table(row_matrix_left_temp[j,])[1]>num_pair){
            print("query genes matrix")
            print(table(row_matrix_left_temp[j,])[1])
            checkpoint=1
        }
    }
    
    row_matrix_right_temp<-matrix(ncol=num_pair_def,nrow=(count-1)*(count-1))
    for (i in 1:(count-1)){
        for (k in 1:(count-1)){
            row_matrix_right_temp[(count-1)*(i-1)+k,]<-row_matrix_whole[k,(num_pair_def*i+1):(num_pair_def*(i+1))]
            if (table(row_matrix_right_temp[(count-1)*(i-1)+k,])[1]>round(num_pair/2)){
                print("partner genes matrix")
                print(table(row_matrix_right_temp[(count-1)*(i-1)+k,])[1])
                checkpoint=1
            }
        }
    }
    print(checkpoint_count)
    checkpoint_count=checkpoint_count+1
}

row_matrix_left<-matrix(ncol=num_pair,nrow=count-1)
for (j in 1:(count-1)){
    shift_base=0
    for (k in 1:num_pair){
        while (row_matrix_left_temp[j,k+shift_base]==0){
            shift_base=shift_base+1
        }
        row_matrix_left[j,k]<-row_matrix_left_temp[j,k+shift_base]
    }
}

row_matrix_right<-matrix(ncol=round(num_pair*3/2),nrow=(count-1)*(count-1))
for (i in 1:(count-1)){
    for (k in 1:(count-1)){
        shift_base=0
        for (l in 1:round(num_pair*3/2)){
            while (row_matrix_right_temp[(count-1)*(i-1)+k,l+shift_base]==0){
                shift_base=shift_base+1
            }
            row_matrix_right[(count-1)*(i-1)+k,l]<-row_matrix_right_temp[(count-1)*(i-1)+k,l+shift_base]
        }
    }
}



#build bed file to extract sequence
bed_matrix<-matrix(ncol=6,nrow=num_pair*2*(count-1)*(count-1))

#matrix for supporting summary tables
pair_full_matrix<-matrix(ncol=12, nrow=num_pair*(count-1)*(count-1))
pair_exon_matrix<-matrix(ncol=2, nrow=num_pair*(count-1)*(count-1))
pair_bp_matrix<-matrix(ncol=6, nrow=num_pair*(count-1)*(count-1))
count_row=0

for (ii in 1:(count-1)){
	#for each left matrix
	for (ij in 1:(count-1)){
		#for each right matrix
		for (ik in 1:num_pair){
			#for each pair
			#check whether they are too close to each other in the while loop. If it is too close, use the next one instead.
			#loca<-as.character(file_bed_sort[row_matrix_left[ii,ik],3])
			chr_left<-as.character(file_bed_sort[row_matrix_left[ii,ik],1])
			range_left<-as.numeric(file_bed_sort[row_matrix_left[ii,ik],3])
			chr_right<-as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik],1])
                        #this is for right end, should start from the smaller number.Also, because of the fastqfrombed problem, if I want to have a correct start, I need to minus 1.
			range_right<-as.numeric(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik],2])-1
			il=0
			while (chr_left==chr_right) {
				if (abs(range_left-range_right)<100000){
					il=il+1
				} else {break}
				if (ik+il<=round(num_pair*3/2)){
					chr_right<-as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il],1])
					range_right<-as.numeric(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il],2])-1
				} else {
					chr_right<-as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il-round(num_pair*3/2)],1])
					range_right<-as.numeric(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il-round(num_pair*3/2)],2])-1
				}
			}
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,1]<-as.character(chr_left)
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,2]<-as.character(range_left-frag_len)
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,3]<-as.character(range_left)
			anno<-strsplit(as.character(file_bed_sort[row_matrix_left[ii,ik],5]),";")
			gene_id<-strsplit(anno[[1]][1]," ")[[1]][3]
			tanscript_id<-strsplit(anno[[1]][2]," ")[[1]][3]
			exon_id<-strsplit(anno[[1]][3]," ")[[1]][3]
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,4]<-as.character(paste(gene_id,tanscript_id,exon_id,ii,sep="_"))
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,5]<-as.character(file_bed_sort[row_matrix_left[ii,ik],4])
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,6]<-as.numeric(file_bed_sort[row_matrix_left[ii,ik],6])
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,1]<-as.character(chr_right)
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,2]<-as.character(range_right)
			bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,3]<-as.character(range_right+frag_len)
			if (ik+il<=num_pair){
				anno<-strsplit(as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il],5]),";")
				gene_id<-strsplit(anno[[1]][1]," ")[[1]][3]
				tanscript_id<-strsplit(anno[[1]][2]," ")[[1]][3]
				exon_id<-strsplit(anno[[1]][3]," ")[[1]][3]
				bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,4]<-as.character(paste(gene_id,tanscript_id,exon_id,ij,sep="_"))
				bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,5]<-as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il],4])
				bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,6]<-as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il],6])
			} else {
				anno<-strsplit(as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il-num_pair],5]),";")
				gene_id<-strsplit(anno[[1]][1]," ")[[1]][3]
				tanscript_id<-strsplit(anno[[1]][1]," ")[[1]][3]
				exon_id<-strsplit(anno[[1]][1]," ")[[1]][3]
				bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,4]<-as.character(paste(gene_id,tanscript_id,exon_id,ij,sep="_"))
				bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,5]<-as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il-num_pair],4])
				bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,6]<-as.character(file_bed_sort[row_matrix_right[(ii-1)*(count-1)+ij,ik+il-num_pair],6])
			}
                        #since they are annotate here already, should try to find a way to count and summarize the simulated fusions.
                        count_row=count_row+1
                        
                        pair_full_matrix[count_row,1:6]<-bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,1:6]
                        pair_full_matrix[count_row,7:12]<-bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,1:6]
                        
                        pair_exon_matrix[count_row,1]<-bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,4]
                        pair_exon_matrix[count_row,2]<-bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,4]
                        
                        pair_bp_matrix[count_row,1]<-bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,1]
                        pair_bp_matrix[count_row,2]<-paste("bp",bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,3],sep="")
                        pair_bp_matrix[count_row,3]<-strsplit(bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2-1,4],"_")[[1]][1]
                        pair_bp_matrix[count_row,4]<-bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,1]
                        pair_bp_matrix[count_row,5]<-paste("bp",bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,3],sep="")
                        pair_bp_matrix[count_row,6]<-strsplit(bed_matrix[((ii-1)*(count-1)*num_pair+(ij-1)*num_pair+ik)*2,4],"_")[[1]][1]
		}
	}
}

#since they are annotate here already, should try to find a way to count and summarize the simulated fusions.


#need to output all matrix for trace back.
write.table(group_matrix, file=paste(file.out,".expression_groups",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
write.table(row_matrix_left, file=paste(file.out,".row_matrix_left",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
write.table(row_matrix_right, file=paste(file.out,".row_matrix_right",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
write.table(bed_matrix, file=paste(file.out,".bed",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
write.table(pair_bp_matrix, file=paste(file.out,"_ref_bp.txt",sep=""),quote=F,sep="_",row.names=F,col.names=F)
write.table(pair_full_matrix, file=paste(file.out,"_full_pairs.bed",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
write.table(pair_exon_matrix, file=paste(file.out,"_exon_pairs.bed",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
