#!/bin/bash

#in this sh, the job are:
#1\a while loop to create path to each sub simulation
#2\in each sub simulation, cp the auto_QF_multi_result_gathering.sh there and run
#3\output the result in to the whole summary for all the sub simulations.
#4\use an R to do calculation on the whole summary and output the result into a table and a few related figures(bar plot) in PDF.
#5\input the subgroup data into corresponding tables.

if [ $# -ne 3 ]
then
  echo ""
    echo "Usage: auto_QF_multi_result_gathering_wrapper.sh simulation_matrix_num_path fusion_simulator_path folder_name_extention"
    echo ""
    echo "simulation_matrix_num_path - matrix of simulation groups."
    echo "fusion_simulator_path - path to fusion simulator"
    echo "folder_name_extention - _QF or _QF_5k_expanded in most cases"    
    echo "example: auto_QF_multi_result_gathering_wrapper.sh /restricted/projectnb/montilab-p/LinGA_unprotected/ytan/ENCODE_simulation_data/simulation_only_split_matrix.sim_out.txt /usr3/graduate/ytan7/CBMrepository/utilities/tags/Fusion_Simulator/ _QF"
    exit 1
fi

if [ ! -f $1 ] 
then
  echo ""
  echo "The file $1 does not exist"
  echo ""
  exit 1
fi

simulation_matrix_num_path=$1
fusion_simulator_path=$2"/"
folder_name_extention=$3
#if [ -s QF_all_sims_summary_table.txt ] 
#then
#	rm -rf QF_all_sims_summary_table.txt
#fi
#echo header 
#echo -e "split_span\tquery_num\tquery_found\tfusion_num\tfusion_all_found\tfusion_TP\tsplit_all_found\tsplit_TP\tspan_all_found\tspan_TP" > QF_all_sims_summary_table.txt
#output header
echo -e "NUM of perfectly found\tNUM of fusion found within 10bp of the simulated breakpoint\tNUM of False Positive found" > QF_simulation_all_summary.txt


ID1="$(head -1 $simulation_matrix_num_path | cut -f1)"; 
ID2="$(head -1 $simulation_matrix_num_path | cut -f2)"; 
first_ID="_"$[ID1]"_"$[ID2]; 
#group_num=`wc -l "coverage_on_exons.txt"$first_ID".expression_groups" | cut -f1 -d" "`
#build the headers first
#for i in `seq 1 $group_num`;
#do
#    for j in `seq 1 $group_num`; 
#    do
#        echo -e "GroupID\tTruePos\tTotalFusion" > "QF_all_sims_summary_"$i"_"$j"_table.txt"
#    done
#done

#for job 1,2,3
while read myfile1; 
    do 
    ID1="$(echo $myfile1 | cut -f1 -d" ")"; 
    ID2="$(echo $myfile1 | cut -f2 -d" ")"; 
    ID="_"$[ID1]"_"$[ID2]; 
    folder_path="coverage_on_exons.txt"$ID$folder_name_extention
    #job1
    cd $folder_path"/results"
    pwd
    #job2
    $fusion_simulator_path"auto_QF_multi_result_gathering.sh" $ID1 $ID2 $fusion_simulator_path
    num_query="$(wc -l temp_folder_list |cut -f1 -d" ")"
    
    #job5
    #input data to subgroup table
    #for i in `seq 1 $group_num`;
    #do
    #    for j in `seq 1 $group_num`; 
    #    do
    #        row_num=`expr $i \* $group_num + $j - $group_num `
    #        TP=`sed -n "${row_num}p" QF_subgroup_summary_table.txt | cut -f2`
    #        TF=`sed -n "${row_num}p" QF_subgroup_summary_table.txt | cut -f3`
    #        echo -e "$ID\t$TP\t$TF" >> "../../QF_all_sims_summary_"$i"_"$j"_table.txt"
    #    done
    #done
    
    #job3
    cd ../../
    
    #grep the needed information
    PF=`sed -n '1p' "QF_on"$ID"_ref_bp_summary" |cut -f2 -d":" | cut -f2 -d" "`
    F10=`sed -n '2p' "QF_on"$ID"_ref_bp_summary" |cut -f2 -d":" | cut -f2 -d" "`
    FP=`sed -n '3p' "QF_on"$ID"_ref_bp_summary" |cut -f2 -d":" | cut -f2 -d" "`
    echo -e $ID"\t"$PF"\t"$F10"\t"$FP >> QF_simulation_all_summary.txt
    echo $i"/QF_simulation_all_summary.txt">>../file_defuse_summary_list.txt
    
    
    #line1=`sed -ne 's/\t/\\\t/g; 2p' $folder_path"/results/QF_sim_summary_table.txt"`
    #echo -e $ID"\t"$line1 >> QF_all_sims_summary_table.txt

done < $simulation_matrix_num_path

#job4
Rscript $fusion_simulator_path"barplot_auto_defuse_multi_summary_only_split.R" file.in=QF_simulation_all_summary.txt file.matrix=$simulation_matrix_num_path file.out=QF_simulation_all_summary num_fusion.in=$num_query

#Rscript $fusion_simulator_path"barplot_auto_QF_multi_result_gathering_only_split.R" file.matrix=$simulation_matrix_num_path file.in=QF_all_sims_summary_table.txt file.out=QF_all_sims_summary

#for i in `seq 1 $group_num`;
#do
#    for j in `seq 1 $group_num`; 
#    do
#        Rscript $fusion_simulator_path"barplot_auto_QF_subgroup_wrapup_only_split.R" file.matrix=$simulation_matrix_num_path file.in="QF_all_sims_summary_"$i"_"$j"_table.txt" file.out="QF_all_sims_summary_"$i"_"$j
#    done
#done
