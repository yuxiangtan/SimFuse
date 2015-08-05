#!/bin/bash

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

#from the simulated fa file, grep to extract the simulated fusion breakpoints.


if [ $# -ne 2 ]
then
  echo ""
    echo "Usage: fusion_bp_extract.sh simulated_fa output_bp"
    echo ""
    echo "simulated_fa - either end fa file of the simulated fusion reads."
    echo "output_bp_prefix - output file name prefix for fusion breakpoints."
    echo "example: fusion_bp_extract.sh coverage_on_exons.txt_50_50_ref.fa1 coverage_on_exons.txt_50_50 "
    exit 1
fi

if [ ! -f $1 ] 
then
  echo ""
  echo "The file $1 does not exist"
  echo ""
  exit 1
fi

#name the parameters
sim_fa=$1
#For the breakpoint information is will be $2_ref_bp.txt"
bp_out=$2"_ref_bp.txt"
#For the queryID information is will be $2_queryID_list.txt
query_out=$2"_queryID_list.txt"
query_count=$2"_queryID_partner_count.txt"

#extract chr_bp_ID (left always forward, right side always reverse)
grep ">" $sim_fa |cut -f2,4,5,11,12,14 -d"_"| sort -u > $bp_out
cut -f3 -d"_" $bp_out | sort -u > $query_out

if [ ! -f $query_count ]
then
    rm -rf $query_count
fi

while read myfile
do
    partner_count=`grep $myfile $bp_out | wc -l` 
    echo -e "$myfile\t$partner_count" >> $query_count
done < $query_out
