#dictionary & file I/O are in lab08

import sys

#compare tophatfusion result to simulated fusion bp to check TP and FP.
#From tophatfusion_result_summary_on_simulation, a few columns needed: chr(3,6) , bp(4,7), ID(2,5), split_num (8), span_num (9).

#Things to consider:
#1\same simulation can have multiple output in tophatfusion. (sort 2ID, left is the small one and right is the bigger one, then use a list to eliminate duplicates)
#2\ differences between duplicate supporting reads. (just take the first one, not doing deep analysis.)

#limiations
#now, I only consider location, but not consider direction at all.
#Only the first event in a duplicate series will be consider.
#because I am using ID as the key, as a result, genes with same location but different strand will be mark as FP: As a result, why do I care about ID?!
#example: simulated: chr10_bp92638891_ENSG00000148688_chr9_bp100112691_ENSG00000197816; tophatfusion: chr10_92638891_ENSG00000148688_chr9_100112691_ENSG00000255036
#So in this version, no ID is used as key.


if __name__ == "__main__":
	#check num of parameter
	if len(sys.argv) < 6:
		#print "Usage: python RFBG_read_generate_from_ref.py random_pair_matrix seq.bed read_len split_num span_num mutation_rate_file mismatch_num=N pair_size pair_std outputfile_name_index"
		print "Usage: python tophatfusion_result_summary_on_simulation.py tophatfusion_result_file fusion_bp_file simulated_split_num simulated_span_num summary_tophatfusion_out"
		print "Example:"
		print "tophatfusion_result_file=result.txt"
		print "data_sim=coverage_on_exons.txt_50_50_ref_bp.txt"
		print "split_num=50"
		print "span_num=50"
		print "sum_out=tophatfusion_on_50_50_ref_bp_summary"
		sys.exit(1)
		
	
	#this version is with no mutation_rate_file (not used anywhere)
	data_TF = sys.argv[1]
	data_sim = sys.argv[2]
	split_num_in = int(sys.argv[3])
	span_num_in = int(sys.argv[4])
	sum_out = sys.argv[5]
	
	
	#open the files
	h_def=open(data_TF,"r")
	h_sim=open(data_sim,"r")
	h_def_dedup_out=open(data_TF+"_dedup.txt","w")
	h_out=open(sum_out,"w")
	h_out_TP=open(sum_out+"_TP.txt","w")
	h_out_FP=open(sum_out+"_FP.txt","w")
	
	#there is no header
	#def_hearder=h_def.readline()
	
	def_matrix=[line.strip().split("\t") for line in h_def]
	sim_matrix=[line.strip().split("_") for line in h_sim]
	h_def.close()
	h_sim.close()
	
	fusion_list=[]
	
	print "Generating deduplicated tophatfusion result"
	#deduplicate for tophatfusion result and extract informations needed:
	for row_num in range(len(def_matrix)):
		#chr(3,6) , bp(4,7), ID(2,5), split_num (8), span_num (9).
		chr_end1=def_matrix[row_num][2]
		chr_end2=def_matrix[row_num][5]
		bp_end1=def_matrix[row_num][3]
		bp_end2=def_matrix[row_num][6]
		ID_end1=def_matrix[row_num][1]
		ID_end2=def_matrix[row_num][4]
		split_num=def_matrix[row_num][7]
		span_num=def_matrix[row_num][8]
		
		#reorder then fusion pairs to make ID1 smaller than ID2
		if chr_end1<chr_end2:
			fusion_pair=chr_end1+"_"+bp_end1+"_"+ID_end1+"_"+chr_end2+"_"+bp_end2+"_"+ID_end2
		else: 	
			if chr_end1>chr_end2:
				fusion_pair=chr_end2+"_"+bp_end2+"_"+ID_end2+"_"+chr_end1+"_"+bp_end1+"_"+ID_end1
			else:
				if bp_end1<bp_end2:
					fusion_pair=chr_end1+"_"+bp_end1+"_"+ID_end1+"_"+chr_end2+"_"+bp_end2+"_"+ID_end2
				else:
					if bp_end1>bp_end2:
						fusion_pair=chr_end2+"_"+bp_end2+"_"+ID_end2+"_"+chr_end1+"_"+bp_end1+"_"+ID_end1
					else:
						print "Warning: there is a fusion pair from tophatfusion has both breakpoint the same"
						continue
		
		#check whether this fusion_pair already existed before
		if not fusion_pair in fusion_list:
			fusion_list.append(fusion_pair)
			h_def_dedup_out.write(fusion_pair+'\t'+split_num+"\t"+span_num+"\n")
	
	h_def_dedup_out.close()
	
	#read back in the dedup file
	h_def_dedup=open(data_TF+"_dedup.txt","r")
	def_dedup_matrix=[line.strip().split("\t") for line in h_def_dedup]
	h_def_dedup.close()
	
	perfect_match=0
	perfect_match_split=0
	perfect_match_span=0
	mismatch_but_within_range_10=0
	mismatch_but_within_range_10_split=0
	mismatch_but_within_range_10_span=0
	FP_temp_list=[]
	#for each row in def_dedup,
	for row_def in range(len(def_dedup_matrix)):
		split_num=int(def_dedup_matrix[row_def][1])
		span_num=int(def_dedup_matrix[row_def][2])
		fusion_pair=def_dedup_matrix[row_def][0]
		chr_end1=fusion_pair.split("_")[0]
		chr_end2=fusion_pair.split("_")[3]
		bp_end1=fusion_pair.split("_")[1]
		bp_end2=fusion_pair.split("_")[4]
		ID_end1=fusion_pair.split("_")[2]
		ID_end2=fusion_pair.split("_")[5]
		#use chr and bp to check whether any match, if no match (count_match==0),then it is a false positive
		count_match=0
		for row_sim in range(len(sim_matrix)):
			if chr_end1 in sim_matrix[row_sim]:
				#check which end it has the same chr combination
				if chr_end1==sim_matrix[row_sim][0] and chr_end2==sim_matrix[row_sim][3]:
					if bp_end1==sim_matrix[row_sim][1].split("bp")[1]:
						#check whether the 2nd end is perfect, withinrange or wrong position.
						if bp_end2==sim_matrix[row_sim][4].split("bp")[1]:
							perfect_match+=1
							perfect_match_split+=split_num
							perfect_match_span+=span_num
							h_out_TP.write("\t".join([fusion_pair,str(split_num),str(span_num),"perfect"])+"\n")
							#check whether it is already in the FP_temp_list
							count_match+=1
							#once a match is found, then go on to next tophatfusion record.
							break
						else:
							if int(bp_end2)>=(int(sim_matrix[row_sim][4].split("bp")[1])-10) and int(bp_end2)<=(int(sim_matrix[row_sim][4].split("bp")[1])+10):
								mismatch_but_within_range_10+=1
								mismatch_but_within_range_10_split+=split_num
								mismatch_but_within_range_10_span+=span_num
								list1=[fusion_pair,str(split_num),str(span_num),"within_range_10"]
								#to merge two list together
								list1+=sim_matrix[row_sim]
								h_out_TP.write("\t".join(list1)+"\n")
								count_match+=1
								break
					else:
						#check two situation: within range or totally another location
						if int(bp_end1)>=(int(sim_matrix[row_sim][1].split("bp")[1])-10) and int(bp_end1)<=(int(sim_matrix[row_sim][1].split("bp")[1])+10):
							if int(bp_end2)>=(int(sim_matrix[row_sim][4].split("bp")[1])-10) and int(bp_end2)<=(int(sim_matrix[row_sim][4].split("bp")[1])+10):
								mismatch_but_within_range_10+=1
								mismatch_but_within_range_10_split+=split_num
								mismatch_but_within_range_10_span+=span_num
								list1=[fusion_pair,str(split_num),str(span_num),"within_range_10"]
								#to merge two list together
								list1+=sim_matrix[row_sim]
								h_out_TP.write("\t".join(list1)+"\n")
								count_match+=1
								break
						else:
							if bp_end1==sim_matrix[row_sim][4].split("bp")[1]:
								if bp_end2==sim_matrix[row_sim][1].split("bp")[1]:
									perfect_match+=1
									perfect_match_split+=split_num
									perfect_match_span+=span_num
									h_out_TP.write("\t".join([fusion_pair,str(split_num),str(span_num),"perfect"])+"\n")
									#check whether it is already in the FP_temp_list
									count_match+=1
									#once a match is found, then go on to next tophatfusion record.
									break
								else:
									if int(bp_end2)>=(int(sim_matrix[row_sim][1].split("bp")[1])-10) and int(bp_end2)<=(int(sim_matrix[row_sim][1].split("bp")[1])+10):
										mismatch_but_within_range_10+=1
										mismatch_but_within_range_10_split+=split_num
										mismatch_but_within_range_10_span+=span_num
										list1=[fusion_pair,str(split_num),str(span_num),"within_range_10"]
										#to merge two list together
										list1+=sim_matrix[row_sim]
										h_out_TP.write("\t".join(list1)+"\n")
										count_match+=1
										break
							else:
								#check two situation: within range or totally another location
								if int(bp_end1)>=(int(sim_matrix[row_sim][4].split("bp")[1])-10) and int(bp_end1)<=(int(sim_matrix[row_sim][4].split("bp")[1])+10):
									if int(bp_end2)>=(int(sim_matrix[row_sim][1].split("bp")[1])-10) and int(bp_end2)<=(int(sim_matrix[row_sim][1].split("bp")[1])+10):
										mismatch_but_within_range_10+=1
										mismatch_but_within_range_10_split+=split_num
										mismatch_but_within_range_10_span+=span_num
										list1=[fusion_pair,str(split_num),str(span_num),"within_range_10"]
										#to merge two list together
										list1+=sim_matrix[row_sim]
										h_out_TP.write("\t".join(list1)+"\n")
										count_match+=1
										break
								
				else:
					if chr_end1==sim_matrix[row_sim][3] and chr_end2==sim_matrix[row_sim][0]:
								if bp_end1==sim_matrix[row_sim][4].split("bp")[1]:
									#check whether the 2nd end is perfect, withinrange or wrong position.
									if bp_end2==sim_matrix[row_sim][1].split("bp")[1]:
										perfect_match+=1
										perfect_match_split+=split_num
										perfect_match_span+=span_num
										h_out_TP.write("\t".join([fusion_pair,str(split_num),str(span_num),"perfect"])+"\n")
										#check whether it is already in the FP_temp_list
										count_match+=1
										#once a match is found, then go on to next tophatfusion record.
										break
									else:
										if int(bp_end2)>=(int(sim_matrix[row_sim][1].split("bp")[1])-10) and int(bp_end2)<=(int(sim_matrix[row_sim][1].split("bp")[1])+10):
											mismatch_but_within_range_10+=1
											mismatch_but_within_range_10_split+=split_num
											mismatch_but_within_range_10_span+=span_num
											list1=[fusion_pair,str(split_num),str(span_num),"within_range_10"]
											#to merge two list together
											list1+=sim_matrix[row_sim]
											h_out_TP.write("\t".join(list1)+"\n")
											count_match+=1
											break
										
								else:
									#check two situation: within range or totally another location
									if int(bp_end1)>=(int(sim_matrix[row_sim][4].split("bp")[1])-10) and int(bp_end1)<=(int(sim_matrix[row_sim][4].split("bp")[1])+10):
											if int(bp_end2)>=(int(sim_matrix[row_sim][1].split("bp")[1])-10) and int(bp_end2)<=(int(sim_matrix[row_sim][1].split("bp")[1])+10):
												mismatch_but_within_range_10+=1
												mismatch_but_within_range_10_split+=split_num
												mismatch_but_within_range_10_span+=span_num
												list1=[fusion_pair,str(split_num),str(span_num),"within_range_10"]
												#to merge two list together
												list1+=sim_matrix[row_sim]
												h_out_TP.write("\t".join(list1)+"\n")
												count_match+=1
												break
		
		if count_match==0:
			FP_str="\t".join([fusion_pair,str(split_num),str(span_num)])
			#add it in only when this str in not already in the list.
			if not FP_str in FP_temp_list:
				FP_temp_list.append(FP_str)


	FP_split=0
	FP_span=0
	#output the FP from the FP_temp_list
	FP=len(FP_temp_list)
	for row_num in range(FP):
		span_num=int(FP_temp_list[row_num].split("\t")[2])
		split_num=int(FP_temp_list[row_num].split("\t")[1])
		FP_split+=split_num
		FP_span+=span_num
		h_out_FP.write(FP_temp_list[row_num]+"\n")
	

	h_out.write("NUM of perfectly found: "+str(perfect_match)+"\n" )
	h_out.write("NUM of fusion found within 10bp of the simulated breakpoint: "+str(mismatch_but_within_range_10)+"\n"  )
	h_out.write("NUM of False Positive found: "+str(FP)+"\n"  )
	if perfect_match!=0:
		if split_num_in!=0:
			h_out.write("Average split reads for perfectly found fusions: "+str(perfect_match_split/perfect_match)+"; detected/simulated%: "+str(perfect_match_split/perfect_match/float(split_num_in))+"\n"  )
		else:
			h_out.write("Average split reads for perfectly found fusions: "+str(perfect_match_split/perfect_match)+"; detected/simulated%: "+"0\n"  )
		if span_num_in!=0:
			h_out.write("Average span reads for perfectly found fusions: "+str(perfect_match_span/perfect_match)+"; detected/simulated%: "+str(perfect_match_span/perfect_match/float(span_num_in))+"\n"  )
		else:
			h_out.write("Average span reads for perfectly found fusions: "+str(perfect_match_span/perfect_match)+"; detected/simulated%: "+"0\n"  )
	else:
		h_out.write("Average split reads for perfectly found fusions: 0; detected/simulated%: NA\n"  )
		h_out.write("Average span reads for perfectly found fusions: 0; detected/simulated%: NA\n"  )
	
	if mismatch_but_within_range_10!=0:
		if split_num_in!=0:
			h_out.write("Average split reads for within 10bp found fusions: "+str(mismatch_but_within_range_10_split/mismatch_but_within_range_10)+"; detected/simulated%: "+str(mismatch_but_within_range_10_split/mismatch_but_within_range_10/float(split_num_in))+"\n"  )
		else:
			h_out.write("Average split reads for within 10bp found fusions: "+str(mismatch_but_within_range_10_split/mismatch_but_within_range_10)+"; detected/simulated%: "+"0\n"  )
		if span_num_in!=0:
			h_out.write("Average span reads for within 10bp found fusions: "+str(mismatch_but_within_range_10_span/mismatch_but_within_range_10)+"; detected/simulated%: "+str(mismatch_but_within_range_10_span/mismatch_but_within_range_10/float(span_num_in))+"\n"  )
		else:
			h_out.write("Average span reads for within 10bp found fusions: "+str(mismatch_but_within_range_10_span/mismatch_but_within_range_10)+"; detected/simulated%: "+"0\n"  )
	else:
		h_out.write("Average split reads for within 10bp found fusions: 0; detected/simulated%: NA\n"  )
		h_out.write("Average span reads for within 10bp found fusions: 0; detected/simulated%: NA\n"  )
	if FP!=0:
		if split_num_in!=0:
			h_out.write("Average split reads for False Positive found fusions: "+str(FP_split/FP)+"; detected/simulated%: "+str(FP_split/FP/float(split_num_in))+"\n"  )
		else:
			h_out.write("Average split reads for False Positive found fusions: "+str(FP_split/FP)+"; detected/simulated%: "+"0\n"  )
		if span_num_in!=0:
			h_out.write("Average span reads for False Positive found fusions: "+str(FP_span/FP)+"; detected/simulated%: "+str(FP_span/FP/float(span_num_in))+"\n"  )
		else:
			h_out.write("Average span reads for False Positive found fusions: "+str(FP_span/FP)+"; detected/simulated%: "+"0\n"  )
	else:
		h_out.write("Average split reads for False Positive found fusions: 0; detected/simulated%: NA\n"  )
		h_out.write("Average span reads for False Positive found fusions: 0; detected/simulated%: NA\n"  )

	
	h_out.close()
	h_out_TP.close()
	h_out_FP.close()