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
file_bed<-read.delim2(file=file.in,sep='\t',header=FALSE,stringsAsFactors=FALSE)
file_bed[,6]<-as.numeric(file_bed[,6])
#sort the file
file_bed_sort<-file_bed[order(file_bed[,6]),]
write.table(file_bed_sort,file=paste(file.in,"_sort",sep=""),quote=F,sep="\t",row.names=F,col.names=F)


