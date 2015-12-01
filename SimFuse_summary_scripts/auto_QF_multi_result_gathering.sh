#!/bin/bash

#in this sh, the job are:
#1\get the folder list and use it into while loop to work for each folder
#2\check whehter whole_fusion_sum_filtered.txt is exist and not empty (this will give information about whether the query is checked)
#3\output all the whole_fusion_sum_filtered.txt into one file then we can use the breakpoint information to check how many of them are TP or FP
#4\split fusion pair list into sub expression groups
#5\get the number of true positive for each subgroups by chr and bp-3digits.

if [ $# -ne 3 ]
then
  echo ""
    echo "Usage: auto_QF_multi_result_gathering.sh split_num span_num"
    echo ""
    echo "split_num"
    echo "span_num"
    echo "fusion_simulator_path"
    echo "example: uto_QF_multi_result_gathering.sh 50 50 /usr3/graduate/ytan7/CBMrepository/utilities/tags/Fusion_Simulator"
    exit 1
fi

split_num=$1
span_num=$2
fusion_simulator_path=$3

if [ ! -d $fusion_simulator_path ] 
then
	echo "The folder $fusion_simulator_path does not exist"
        exit 1
fi

#for job1:
find . -name "ENS*" > temp_folder_list

#rm the old/existed output file
if [ -s query_found.txt ] 
then
	rm -rf query_found.txt
fi

if [ -s all_sum_file.txt ] 
then
	rm -rf all_sum_file.txt
fi

if [ -s TP_fusion_list.txt ] 
then
	rm -rf TP_fusion_list.txt
fi

if [ -s QF_sim_summary_table.txt ] 
then
	rm -rf QF_sim_summary_table.txt
fi

#echo header of QF_sim_summary_table.txt
echo -e "query_num\tquery_found\tfusion_num\tfusion_all_found\tfusion_TP\tsplit_all_found\tsplit_TP\tspan_all_found\tspan_TP" > QF_sim_summary_table.txt
echo -e "GroupID\tTruePos\tTotalFusion" > QF_subgroup_summary_table.txt

#for job2 and 3:
while read myfile;
do
    #now, whole_fusion_sum_filtered.txt is the final output, if it is change, this should be changed.
    whole_sum_file=$myfile"/whole_fusion_sum_filtered.txt"

    #check whether the file exist and not empty
    if [ -s $whole_sum_file ] 
    then
            echo $myfile >> query_found.txt
            #cat $whole_sum_file >> all_sum_file.txt
            #because I have header in the output now
            sed '1d' $whole_sum_file >> all_sum_file.txt
    fi
done < temp_folder_list
#can use the following script to find what queries were missed.
#cut -f2 -d"/" coverage_on_exons.txt_50_50_QF/results/query_found.txt > coverage_on_exons.txt_50_50_QF/results/query_found.txt_cp
#Rscript /protected/individuals/ytan7/CBMrepository/utilities/tags/Fusion_detect_tools_scripts/setdiff_ID.R file.listA="coverage_on_exons.txt_50_50_queryID_list.txt" file.listB="coverage_on_exons.txt_50_50_QF/results/query_found.txt_cp" file.out="coverage_on_exons.txt_50_50_QF/results/query_missed.txt"

#for job4
python $fusion_simulator_path"full_pair_bed_split_groups.py" "../../coverage_on_exons.txt_"$split_num"_"$span_num"_full_pairs.bed" "../../coverage_on_exons.txt_"$split_num"_"$span_num"_ref_bp"

echo `date`
if [ -s all_sum_file.txt ] 
then
    #while read myfile1; 
    #do
    #chr1="$(echo $myfile1 | cut -f1 -d"_")"; 
    #chr2="$(echo $myfile1 | cut -f4 -d"_")"; 
    #BP1="$(echo $myfile1 | cut -f2 -d"_" | cut -f2 -d"p")";
    ##get the length of the string
    #BP1_LEN=${#BP1}
    #BP2="$(echo $myfile1 | cut -f5 -d"_" | cut -f2 -d"p")"; 
    ##get the length of the string
    #BP2_LEN=${#BP2}
    #BP1_use="$(echo $BP1 | cut -b -$[$BP1_LEN-3])";
    #BP2_use="$(echo $BP2 | cut -b -$[$BP2_LEN-3])";
    #grep $chr1 all_sum_file.txt | grep $chr2 |grep $BP1_use | grep $BP2_use >> TP_fusion_list.txt
    #
    #done < "../../coverage_on_exons.txt_"$split_num"_"$span_num"_ref_bp.txt"
    python $fusion_simulator_path"QF_result_summary_on_simulation_not_considerID.py" "all_sum_file.txt" "../../coverage_on_exons.txt_"$split_num"_"$span_num"_ref_bp.txt" $split_num $span_num "../../QF_on_"$split_num"_"$span_num"_ref_bp_summary"
        
    #job5
    #group_num=`wc -l "../../coverage_on_exons.txt_"$split_num"_"$span_num".expression_groups" | cut -f1 -d" "`
    #for i in `seq 1 $group_num`;
    #do
    #    for j in `seq 1 $group_num`; 
    #    do
    #        if [ -s "TP_fusion_list_"$i"_"$j".txt" ] 
    #        then
    #            rm -rf "TP_fusion_list_"$i"_"$j".txt"
    #        fi
    #        while read myfile1;
    #        do
    #        chr1="$(echo $myfile1 | cut -f1 -d"_")"; 
    #        chr2="$(echo $myfile1 | cut -f4 -d"_")";
    #        ENS1="$(echo $myfile1 | cut -f1 -d"_")"; 
    #        ENS2="$(echo $myfile1 | cut -f4 -d"_")"; 
    #        BP1="$(echo $myfile1 | cut -f2 -d"_" | cut -f2 -d"p")";
    #        #get the length of the string
    #        BP1_LEN=${#BP1}
    #        BP2="$(echo $myfile1 | cut -f5 -d"_" | cut -f2 -d"p")"; 
    #        #get the length of the string
    #        BP2_LEN=${#BP2}
    #        BP1_use="$(echo $BP1 | cut -b -$[$BP1_LEN-2])";
    #        BP2_use="$(echo $BP2 | cut -b -$[$BP2_LEN-2])";
    #        #try to get all fusions found in each group as the TP+FP, but failed in this way
    #        #grep $ENS1 all_sum_file.txt | grep $ENS2 >> "all_fusion_found_list_"$i"_"$j".txt"
    #        grep $chr1 all_sum_file.txt | grep $chr2 |grep $BP1_use | grep $BP2_use >> "TP_fusion_list_"$i"_"$j"_bp.txt"
    #        grep $chr1 all_sum_file.txt | grep $chr2 |grep $BP1_use | grep $BP2_use | cut -f5,10 | sort -u >> "TP_fusion_list_"$i"_"$j".txt"
    #        done < "../../coverage_on_exons.txt_"$split_num"_"$span_num"_ref_bp_"$i"_"$j".txt"
    #    done
    #done

else
    echo "">all_sum_file.txt
    python $fusion_simulator_path"QF_result_summary_on_simulation_not_considerID.py" "all_sum_file.txt" "../../coverage_on_exons.txt_"$split_num"_"$span_num"_ref_bp.txt" $split_num $span_num "../../QF_on_"$split_num"_"$span_num"_ref_bp_summary"
    echo "all_sum_file.txt is empty, which means no fusions found."
    exit 1
fi

#echo `date`
#cut -f5,10 TP_fusion_list.txt | sort -u > TP_fusion_ID_pairlist.txt
#
#query_num=`wc -l "../../coverage_on_exons.txt_"$split_num"_"$span_num"_queryID_list.txt" | cut -f1 -d" "`
#query_found=`wc -l query_found.txt | cut -f1 -d" "`
#fusion_num=`wc -l "../../coverage_on_exons.txt_"$split_num"_"$span_num"_ref_bp.txt" | cut -f1 -d" "`
#fusion_all_found=`cut -f5,10 all_sum_file.txt |sort -u |wc -l | cut -f1 -d" "`
#fusion_TP=`wc -l TP_fusion_ID_pairlist.txt | cut -f1 -d" "`
#split_all_found=`awk '{print $(11)}' all_sum_file.txt | xargs |sed 's/ /+/g' | bc`
#split_TP=`awk '{print $(11)}' TP_fusion_list.txt | xargs |sed 's/ /+/g' | bc`
#span_all_found=`awk '{print $(12)}' all_sum_file.txt | xargs |sed 's/ /+/g' | bc`
#span_TP=`awk '{print $(12)}' TP_fusion_list.txt | xargs |sed 's/ /+/g' | bc`
#
#echo -e "$query_num\t$query_found\t$fusion_num\t$fusion_all_found\t$fusion_TP\t$split_all_found\t$split_TP\t$span_all_found\t$span_TP" >> QF_sim_summary_table.txt

#generate able for subgroups
#since I am not considering subgroup, this is not used anymore.
#for i in `seq 1 $group_num`;
#do
#    for j in `seq 1 $group_num`; 
#    do
#        ID="_"$i"_"$j
#        #all_found=`wc -l "all_fusion_found_list_"$i"_"$j".txt"`
#        TP=`wc -l "TP_fusion_list_"$i"_"$j".txt"| cut -f1 -d" "`
#        TF=`expr $fusion_num / $group_num / $group_num`
#        echo -e "$ID\t$TP\t$TF" >> QF_subgroup_summary_table.txt
#    done
#done

#Rscript $fusion_simulator_path"barplot_auto_QF_subgroup_gathering.R" file.in=QF_subgroup_summary_table.txt file.out=QF_subgroup_summary group.num=$group_num
