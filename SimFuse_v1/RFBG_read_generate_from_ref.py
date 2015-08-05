#dictionary & file I/O are in lab08

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


import sys
import string
import time
import random

#Another way to do the simulation by using ref sequence and add mutation by mutation rate estimate.
#by this method, a fix number of splitting read and spanning read will be generated.
##need the mutation rate file (liye's method)
##for each read generated, randomly pick N locations and then random pick mutated or not by rate. If mutated, then mark 1, not then 0.
##in each support read ID, put the mutation/mismatch number at the end of the read ID.
##when generate reads, need to find a way to extract seq efficiently each time, and then do str substitution.

#note, fraglen (means the input range of reads) must >=read_len*3+insertsize.

#import the reverse complementary function.	
def rc(dna):
    import string
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq


if __name__ == "__main__":
	#check num of parameter
	if len(sys.argv) < 11:
		#print "Usage: python RFBG_read_generate_from_ref.py random_pair_matrix seq.bed read_len split_num span_num mutation_rate_file mismatch_num=N pair_size pair_std outputfile_name_index"
		print "Usage: python RFBG_read_generate_from_ref.py random_pair_matrix seq.bed read_len split_num span_num mismatch_num=N pair_size pair_std outputfile_name_index exclud_len"
		sys.exit(1)
		
	
	#this version is with no mutation_rate_file (not used anywhere)
	data_pair = sys.argv[1]
	data_seq = sys.argv[2]
	read_len = int(sys.argv[3])
	split_num = int(sys.argv[4])
	span_num = int(sys.argv[5])
	#data_mutate_rate = sys.argv[6] this is not used at all.
	mismatch_num = int(sys.argv[6])
	pair_size = int(sys.argv[7])
	pair_std = int(sys.argv[8])
	outfile = sys.argv[9]
	exclud_len = int(sys.argv[10])
	
	if pair_size-read_len*2<5:
		print "The pair_size is too small to generate spanning reads efficiently."
		sys.exit(1)
	
	if read_len-exclud_len*2<10:
		print "The exclud_len is too large comparing to the read length. (This means read_len-exclud_len*2<10)"
		sys.exit(1)
		
	#open the files
	h_pair=open(data_pair,"r")
	h_seq=open(data_seq,"r")
	#h_mutate=open(data_mutate_rate,"r")
	h_out1=open(outfile+".fa1","w")
	h_out2=open(outfile+".fa2","w")
	
	pair_matrix=[line.strip().split("\t") for line in h_pair]
	SEQ_matrix=[line.strip().split("\t") for line in h_seq]
	
	if len(pair_matrix) != len(SEQ_matrix) :
		print "The random_pair_matrix may be a wrong pair with seq.bed."
		sys.exit(1)
	
	#the key is how to mimic the normal distribution or insert size. Get the breakpoint and use normal to get the size, then uniformly pick the breaklocation on the insert region.
	
	#number of pairs to generate
	PAIRS=len(pair_matrix)/2
	print(PAIRS)
	print(split_num)
	
	for row_num in range(PAIRS):
		#print(row_num)
		END1=row_num*2
		END2=row_num*2+1
		#check whether the row and ID match between two files.
		if pair_matrix[END1][3]!= SEQ_matrix[END1][0] or pair_matrix[END2][3]!= SEQ_matrix[END2][0]:
			print "The random_pair_matrix may be a wrong pair with seq.bed."
			print ("row_num is:"+str(row_num))
			sys.exit(1)
		#get info from matrix:
		chr_end1= pair_matrix[END1][0]
		#breakpoint on end1 is bigger side
		bp_end1= pair_matrix[END1][2]
		str_end1=pair_matrix[END1][4]
		chr_end2= pair_matrix[END2][0]
		str_end2=pair_matrix[END2][4]
		#breakpoint on end2 is smaller side
		bp_end2= pair_matrix[END2][1]
		frag_len1=len(SEQ_matrix[END1][1])
		frag_len2=len(SEQ_matrix[END2][1])
		#generate splitting read first
		#split reads must be at least exclud_len across breakpoint
		for n_pair in range(split_num):		
			#generate the breakpoint loc on a read from end1
			#print(n_pair)
			BP=random.randint(exclud_len,read_len-exclud_len)
			BP_END1=frag_len1-BP
			BP_END2=read_len-BP
			#because of the fastqfrombed problem,the correct breakpoint mark on end2, must +1 to show the correct number.
			READ_ID = ">HWI_"+pair_matrix[END1][0]+"_"+pair_matrix[END1][1]+"_"+"bp"+pair_matrix[END1][2]+"_"+pair_matrix[END1][3]+"_"+pair_matrix[END1][4]+"_"+pair_matrix[END1][5]+"_"+pair_matrix[END2][0]+"_"+"bp"+str(int(pair_matrix[END2][1])+1)+"_"+str(int(pair_matrix[END2][2])+1)+"_"+pair_matrix[END2][3]+"_"+pair_matrix[END2][4]+"_"+pair_matrix[END2][5]+"_"+"breaklen"+str(BP)+"_"+str(n_pair)
			SPLIT_END1=SEQ_matrix[END1][1][BP_END1:frag_len1]
			SPLIT_END2=SEQ_matrix[END2][1][0:BP_END2]
			SPLIT_READ=SPLIT_END1+SPLIT_END2
			#randomly pick which end the span come from
			END_SPAN_PICK=random.randint(1,2)
			
			if END_SPAN_PICK==1:
				random_insert_size=int(round(random.gauss((pair_size-read_len*2),pair_std)))
				START_loc=BP_END1-random_insert_size-read_len
				if START_loc<0:
					random_insert_size=BP_END1-read_len
					START_loc=0
				if START_loc>(frag_len1-read_len):
					START_loc=frag_len1-read_len
				SPAN_READ=SEQ_matrix[END1][1][START_loc:(read_len+START_loc)]
				#because all ref from + strand, if end1 gene should from + stand, then no rc (reverse complement), but if not, then need to rc, end2 should always on the opposite strand of end1. because of the sequencing potocal
				#if str(str_end1)=="+":
				h_out1.write(READ_ID+"_SPLITEND2_insertsize_"+str(random_insert_size)+"\n"+SPAN_READ+"\n")
				h_out2.write(READ_ID+"_SPLITEND2_insertsize_"+str(random_insert_size)+"\n"+rc(SPLIT_READ)+"\n")
				
			else:
				random_insert_size=int(round(random.gauss((pair_size-read_len*2),pair_std)))
				START_loc=BP_END2+random_insert_size
				if START_loc<0:
					START_loc=0
				if START_loc>(frag_len2-read_len):
					START_loc=(frag_len2-read_len)
				SPAN_READ=SEQ_matrix[END2][1][START_loc:read_len+START_loc]
				#if str(str_end1)=="+":
				h_out1.write(READ_ID+"_SPLITEND1_insertsize_"+str(random_insert_size)+"\n"+SPLIT_READ+"\n")
				h_out2.write(READ_ID+"_SPLITEND1_insertsize_"+str(random_insert_size)+"\n"+rc(SPAN_READ)+"\n")
				
		#generate spanning read 
		for n_span in range(span_num):
			#print(n_span)
			#new way to do this, get the insertsize first
			#To keep the insert size >0 because this is a spanning pair. However, in this case, it is not full follow normal.
			random_insert_size=-1
			#use this to bound the distribution and try to control the insertsize distribute as expected
			while random_insert_size<0 or random_insert_size>(pair_size-read_len*2)*2:
				random_insert_size=int(round(random.gauss((pair_size-read_len*2),pair_std)))
			#pick uniformly a loc on this insert_size string as the breakpoint
			BP=random.randint(0,random_insert_size)
			#the SPANFRAG indicate the insert_size on end1=BP
			READ_SPAN_ID = ">HWI_"+pair_matrix[END1][0]+"_"+pair_matrix[END1][1]+"_"+"bp"+pair_matrix[END1][2]+"_"+pair_matrix[END1][3]+"_"+pair_matrix[END1][4]+"_"+pair_matrix[END1][5]+"_"+pair_matrix[END2][0]+"_"+"bp"+str(int(pair_matrix[END2][1])+1)+"_"+pair_matrix[END2][2]+"_"+pair_matrix[END2][3]+"_"+pair_matrix[END2][4]+"_"+pair_matrix[END2][5]+"_"+"SPANFRAG_"+str(BP)+"_insertsize_"+str(random_insert_size)+"_"+str(n_span)
			#if frag_len1 is too short, SPAN_END_LOC1 will be smaller than 0, then it will get nothing, so need to avoid it.
			SPAN_END1_LOC=frag_len1-BP-read_len
			if SPAN_END1_LOC<0:
				print("warning: input fragments are too short on end1, insert_size adjusted for "+READ_SPAN_ID)
				SPAN_END1_LOC=0
				BP=frag_len1-read_len
			#if frag_len2 is too short, SPAN_END_LOC2 will be bigger than frag_len2-read_len, then it will get truncated reads, so need to avoid it.
			SPAN_END2_LOC=random_insert_size-BP
			if SPAN_END2_LOC>frag_len2-read_len:
				print("warning: input fragments are too short on end2, insert_size adjusted for "+READ_SPAN_ID)
				SPAN_END2_LOC=frag_len2-read_len
			SPAN_END1_READ=SEQ_matrix[END1][1][SPAN_END1_LOC:(read_len+SPAN_END1_LOC)]
			SPAN_END2_READ=SEQ_matrix[END2][1][SPAN_END2_LOC:(read_len+SPAN_END2_LOC)]
			#if str(str_end1)=="+":
			h_out1.write(READ_SPAN_ID+"\n"+SPAN_END1_READ+"\n")
			h_out2.write(READ_SPAN_ID+"\n"+rc(SPAN_END2_READ)+"\n")
			
				
				
	h_pair.close()
	h_seq.close()
	h_out1.close()
	h_out2.close()