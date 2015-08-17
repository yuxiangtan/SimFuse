# -*- coding: cp936 -*-

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

"""
This script is a fusion read simulator for multi random repeat
However, because the merged file is large, this simulation will have space issue in that case. As a result, in this multirun script, no mearge file will be generated.
Merge file needed for deFuse, tophat, etc need to be generated at the run and then removed.

This script will finally generate an branch of simulated reads with/out merging with proper aligned background in fa and fq formats.
=============================
Usage: python SimFuse.py
-h help

-i file contains a number each row, which will be used to generate split/span read follow read length distribution  *[No default value]

-b Targeted bam file location									                    *[No default value]

-o working/output directory [all folders should have / at the end]                                                  *[No default value]

-c bam filter script to get background reads with no potential fusion reads (only properly paired reads)            [default value is clear_bg_filter in the SF folder]

-e The txt file with all exon annotations for a genome, can be generated from biomart		                    *[No default value]

-G tophat_genome_reference_fa - the path of the genome fa file (such as hg19.fa)                                    *[No default value]

-g LOG_folder                                                                                                       *[No default value]

-F simulator path                                                                                                   *[No default value]

-l Read_length - length of reads		                                                                    [default value 99]

-r resume_status: check whether user want to skip finished step or start over                                       [default value 0, not resume]

-p number of simulation gene pairs in each expression group.                                                        [default value is 100]

-d minimum number of exons for each expression group to sample from 		                                    [default value is 100*100(from -p)]

-m number of mutation allowed in each read (this is a in progress function                                          [default value is 0]

-M fragment size from library preparation                                                                           [default value is 0]

-s standard deviation of fragment size                                                                              [default value is 0]

-E the length on each end of read you want to exclude when doing the simulation					    [default value is 30 because blat will no align read <=30]

-t number of simulation rounds 											    [default value is 100]
============================

Python & Module requirement:
Versions: 2.7 or above
Module: No additional Python Module is required.

============================
Library file requirement:
Not Standalone version, few library file is required.
bowtie2 is need to reinstall and preloaded.
============================

"""

##By Yuxiang Tan
##Contact: yuxiang.tan@gmail.com
##Compatible Python Version:2.7 or above

###Code Framework

#resume the run function
def resume_func(cmd, resume_stat, step_name, next_step_name, LOG_OUT):
    if resume_stat == 1:
        h_LOG_OUT=open(LOG_OUT,"r")
        count_res = 0
        for line in h_LOG_OUT:
	    if next_step_name in line:
                count_res +=1
	h_LOG_OUT.close()            
        if count_res==0:
            subprocess.call(cmd,shell=True)
	    resume_stat=0
        else:
            print step_name+" pass"
    else:
        log_outf=open(LOG_OUT,"a")
        log_outf.write(step_name+"\n")
        subprocess.call(cmd,shell=True)
        log_outf.close()
    return resume_stat

if __name__ == "__main__":
	###Python General Module Import	
	import sys, csv, getopt, re, subprocess, time 
	import os
	import math
	from itertools import ifilter
	
	OUTPUT_SEP_CHAR='\t'
	
	
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 9:
		print __doc__
		sys.exit(0)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	work_folder=None
	bam_name=None
        LOG_F=None
        script_path=None
        exon_txt=None
        genome_fa=None
        read_matrix_list=None
       
        read_len=99
	resume_stat=0
        num_pair_sim="100"
	boundary_read="100000"
	Align_percent=98
        mutate_num="0"
	frag_mean="150"
	frag_std="50"
	exclud_len="30"
	times=100
	clear_bg_filter=None
	
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:b:o:c:g:e:G:F:l:r:p:d:s:M:m:E:t:z')
        for opt in optlist:
            if opt[0] == '-h':
                print __doc__; sys.exit(0)
            elif opt[0] == '-i': read_matrix_list =opt[1]
	    elif opt[0] == '-b': bam_name = opt[1]
            elif opt[0] == '-o': work_folder = opt[1]
            elif opt[0] == '-c': clear_bg_filter = opt[1]
            elif opt[0] == '-g': LOG_F = opt[1]
            elif opt[0] == '-e': exon_txt =opt[1]
            elif opt[0] == '-G': genome_fa = opt[1]
            elif opt[0] == '-F': script_path = opt[1]
            elif opt[0] == '-l': read_len = int(opt[1])
            elif opt[0] == '-r': resume_stat = int(opt[1])
            elif opt[0] == '-p': num_pair_sim =opt[1]
            elif opt[0] == '-d': boundary_read =opt[1]
            elif opt[0] == '-m': mutate_num =opt[1]
	    elif opt[0] == '-M': frag_mean =opt[1]
	    elif opt[0] == '-s': frag_std =opt[1]
	    elif opt[0] == '-E': exclud_len = opt[1]
	    elif opt[0] == '-t': times = int(opt[1])
	    
	if not os.path.exists(work_folder):
	    os.mkdir(work_folder)
	
        if LOG_F==None:
            print "LOG_F is not provided"; sys.exit(0)
        
        if not os.path.exists(LOG_F):
	    os.mkdir(LOG_F)
                
        LOG_ERR=LOG_F+"/error.log"
        LOG_OUT=LOG_F+"/out.log"
	LOG_BAM_FQ=LOG_F+"/bamtofastq.log"
        
        log_error=open(LOG_ERR,"a")
	log_out=open(LOG_OUT,"a")
	log_bam_fq=open(LOG_BAM_FQ,"a")
       
        #for parameter input needed.
        if bam_folder==None:
            print "bam_folder is not provided"
            log_error.write("bam_folder is not provided\n"); sys.exit(1)
        
        if work_folder==None:
            print "work_folder is not provided"
            log_error.write("work_folder is not provided\n"); sys.exit(1)
        
        if bam_name==None:
            print "bam_name is not provided"
            log_error.write("bam_name is not provided\n"); sys.exit(1)
                
        if clear_bg_filter==None:
            clear_bg_filter=script_path+"/clear_bg_filter"
	    
        if exon_txt==None:
            print "exon_txt is not provided"
            log_error.write("exon_txt is not provided\n"); sys.exit(1)
                
        if genome_fa==None:
            print "genome_fa path is not provided"
            log_error.write("genome_fa path is not provided\n"); sys.exit(1)
                
        if script_path==None:
            print "simulator script_path info is not provided"
            log_error.write("simulator script_path is not provided\n"); sys.exit(1)   
        
        if read_matrix_list==None:
            print "read_matrix_list is not provided"
            log_error.write("read_matrix_list is not provided\n"); sys.exit(1)   
        
        
	#check whether the files provide is there.
	
        if not os.path.exists(genome_fa):
	    print "genome_fa not found"
            log_error.write("genome_fa not found\n"); sys.exit(1)
        
        if not os.path.exists(bam_name):
	    print "bam not found"
            log_error.write("bam not found\n"); sys.exit(1)
        
	if not os.path.exists(clear_bg_filter):
	    print "clear_bg_filter script not found"
            log_error.write("clear_bg_filter script not found\n"); sys.exit(1)
        
        if not os.path.exists(exon_txt):
	    print "exon_txt of the target genome not found"
            log_error.write("exon_txt of the target genome not found\n"); sys.exit(1)
        
        if not os.path.exists(read_matrix_list):
	    print "read_matrix_list not found"
            log_error.write("read_matrix_list not found\n"); sys.exit(1)
                 
                
        if not os.path.exists(script_path):
	    print "simulator script_path not found"
            log_error.write("simulator script_path not found\n"); sys.exit(1)

	resume_stat_loc=resume_stat
	coverage_on_exons=work_folder+"coverage_on_exons.txt"
	search_region_len=int(read_len)+int(frag_mean)           
		
	
	
	
	#preprocess files
        step_name="split bam into 3 files, in fusion_simulator.py"
        print "===================="
        print step_name
        next_step_name="Filter bam to get clear bg, in fusion_simulator.py"
        splitbam_cmd="python "+script_path+"/split_bam_no_align.py "+bam_name+" -o "+work_folder
	resume_stat_loc=resume_func(splitbam_cmd, 1, step_name, next_step_name, LOG_OUT)
        print "finished splitting bam into 3 files"
        print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
        print "===================="

        step_name="Filter bam to get clear bg, in fusion_simulator.py"
        print "===================="
        print step_name
        next_step_name="Get raw count expr, in fusion_simulator.py"
        filter_bam_cmd="bamtools filter -in "+work_folder+"paired.bam -out "+work_folder+"proper_pair_no_skip.bam -script "+clear_bg_filter+" ;"
        h_LOG_OUT=open(LOG_OUT,"r")
        count_res = 0
        for line in h_LOG_OUT:
	    if next_step_name in line:
                count_res +=1
	h_LOG_OUT.close()            
        if count_res==0:
            subprocess.call(filter_bam_cmd,shell=True,stdout=log_bam_fq, stderr=log_bam_fq)
        else:
            print step_name+" pass"
        print "finished filterring bam to get clear bg"
        print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
        print "===================="

        step_name="Get raw count expr, in fusion_simulator.py"
        print "===================="
        print step_name
        next_step_name="Get fastq from bam, in fusion_simulator.py"
        get_expr_cmd="coverageBed -abam "+work_folder+"proper_pair_no_skip.bam -b "+exon_txt+" > "+coverage_on_exons+"; "+"Rscript "+script_path+"sort_raw_count_based_random_background.R file.in="+coverage_on_exons
        resume_stat_loc=resume_func(get_expr_cmd, 1, step_name, next_step_name, LOG_OUT)
        print "finished getting raw count expr"
        print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
        print "===================="
	
	step_name="Get fastq from bam, in fusion_simulator.py"
        print "===================="
        print step_name
        next_step_name="Generate different simulation supporting reads (3,5,10)^2 and merge with bg, in fusion_simulator.py"
	get_fastq_cmd="bam2fastq "+work_folder+"proper_pair_no_skip.bam --force -o "+work_folder+"proper_pair_no_skip#.fq"
	resume_stat_loc=resume_func(get_fastq_cmd, 1, step_name, next_step_name, LOG_OUT)
        print "finished getting fastq from bam"
        print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
        print "===================="
	
	#for loop to do multiple random simulations from the same data
	for round_time in range(1,times+1):
	    ##need to reorder folder name and logs.
	    folder_name=work_folder+"/"+str(round_time)+"/"
	    LOG_F_round=LOG_F+"/"+str(round_time)+"/"
	    if not os.path.exists(folder_name):
		os.mkdir(folder_name)
	    if not os.path.exists(LOG_F_round):
		os.mkdir(LOG_F_round)
	    LOG_OUT_round=LOG_F_round+"/out.log"
	    LOG_ERR_round=LOG_F_round+"/error.log"
	    
	    log_error_round=open(LOG_ERR_round,"a")
	    log_out_round=open(LOG_OUT_round,"a")
	    
	    #check existance of final output, if exist, remove it then can generate a new one
	    out_all_merge_fq1=folder_name+"coverage_on_exons.txt"+"all_merge.1.fastq"
	    out_all_merge_fq2=folder_name+"coverage_on_exons.txt"+"all_merge.2.fastq"
	    if os.path.exists(out_all_merge_fq1):
		print out_all_merge_fq1+" already exist, removing"
		os.remove(out_all_merge_fq1)
		
	    if os.path.exists(out_all_merge_fq2):
		print out_all_merge_fq2+" already exist, removing"
		os.remove(out_all_merge_fq2)
		
	    read_rate=round(int(frag_mean))/round(read_len)/2-1
	    if read_rate<=0:
		read_rate=round(int(frag_mean)+int(frag_std))/round(read_len)/2-1
	    
	    if read_rate<=0:
		read_rate=0.05
	    
	    #for loop to submit multiple groups.
	    step_name="Generate different simulation supporting reads (3,5,10)^2 and merge with bg, in fusion_simulator.py, round "+str(round_time)
	    print "===================="
	    print step_name
	    next_step_name="This is the last step in this script"
	    h_read_matrix_list=open(read_matrix_list,"r")
	    h_read_matrix_out=open(read_matrix_list[:-3]+"sim_out.txt","w")
	    read_matrix=[line.strip().split("\t") for line in h_read_matrix_list]  
	    for i_split in read_matrix:
		i_span=int(round(round(int(i_split[0]))*read_rate))
		
		#file name formats:
		out_name=folder_name+"coverage_on_exons.txt"+"_"+str(i_split[0])+"_"+str(i_span)
		out_bed=folder_name+"coverage_on_exons.txt"+"_"+str(i_split[0])+"_"+str(i_span)+".bed"
		ref_bed=folder_name+"coverage_on_exons.txt"+"_"+str(i_split[0])+"_"+str(i_span)+"_ref.bed"
		out_index=folder_name+"coverage_on_exons.txt"+"_"+str(i_split[0])+"_"+str(i_span)+"_ref"
		out_merge_fq1=folder_name+"coverage_on_exons.txt"+"_"+str(i_split[0])+"_"+str(i_span)+"_ref_merged.1.fastq"
		out_merge_fq2=folder_name+"coverage_on_exons.txt"+"_"+str(i_split[0])+"_"+str(i_span)+"_ref_merged.2.fastq"
		
		print "for "+str(i_split[0])+" splits and "+str(i_span)+" spans in round "+str(round_time)
		#first round to generate all simulated reads into fa format
		list_gen_cmd="Rscript "+script_path+"raw_count_based_random_background_exon_list_generator.R file.in="+coverage_on_exons+" num_pair="+str(num_pair_sim)+" read_bound="+str(boundary_read)+" file.out="+out_name+" frag_len="+str(search_region_len)
		fasta_from_bed_cmd="fastaFromBed -name -tab -fi "+genome_fa+" -bed "+out_bed+" -fo "+ref_bed
		read_simulation_cmd="python "+script_path+"RFBG_read_generate_from_ref.py "+out_bed+" "+ref_bed+" "+str(read_len)+" "+str(i_split[0])+" "+str(i_span)+" "+str(mutate_num)+" "+str(frag_mean)+" "+str(frag_std)+" "+out_index+" "+exclud_len
		first_round_cmd=list_gen_cmd+"; "+fasta_from_bed_cmd+"; "+read_simulation_cmd
		resume_stat_loc=resume_func(first_round_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT_round)
		
		#second round, fa to fq and merge.
		fa_to_fq1_cmd="perl "+script_path+"fasta_to_fastq.pl "+out_index+".fa1 > "+out_index+".fq1"
		fa_to_fq2_cmd="perl "+script_path+"fasta_to_fastq.pl "+out_index+".fa2 > "+out_index+".fq2"
		#merge each with proper reads
		merge1_cmd="cat "+out_index+".fq1 "+work_folder+"proper_pair_no_skip_1.fq > "+out_merge_fq1
		merge2_cmd="cat "+out_index+".fq2 "+work_folder+"proper_pair_no_skip_2.fq > "+out_merge_fq2
		#merge all
		merge_1_all_cmd="cat "+out_index+".fq1 >> "+out_all_merge_fq1
		merge_2_all_cmd="cat "+out_index+".fq2 >> "+out_all_merge_fq2
		
		#get some summary file for each profile
		summary_cmd="sh "+script_path+"fusion_bp_extract.sh "+out_index+".fa1"+" "+out_name
		#in order to save space skip merge
		second_round_cmd=fa_to_fq1_cmd+"; "+fa_to_fq2_cmd+"; "+summary_cmd+";"
		resume_stat_loc=resume_func(second_round_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT_round)
		    
		h_read_matrix_out.write(str(i_split[0])+"\t"+str(i_span)+"\n")
		    
		    
	    h_read_matrix_list.close()
	    h_read_matrix_out.close()
	    #merge all with proper reads
	    merge_1_to_all_cmd="cat "+work_folder+"proper_pair_no_skip_1.fq >> "+out_all_merge_fq1
	    merge_2_to_all_cmd="cat "+work_folder+"proper_pair_no_skip_2.fq >> "+out_all_merge_fq2
	    final_cmd=merge_1_to_all_cmd+"; "+merge_2_to_all_cmd
	    #in order to save space skip this merge step
	    #resume_stat_loc=resume_func(final_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT_round)
	    print "finished generating different simulation supporting reads (3,5,10)^2 and merging with bg in round "+str(round_time) 
	    print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
	    print "===================="
	
	