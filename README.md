#SimFuse
<li>SimFuse is an RNA-Seq based fusion simulator. It uses real data to generate background reads. With one simple command, SimFuse can randomly generate multiple fusion samples with different numbers of fusion supporting reads. With enough sampling, users can minimize random effects and also check the performance of a fusion detection algorithm in more detail, such as calculating the differences of recall and precision rates between fusions with different numbers of supporting reads.If the result of a fusion detection algorithm can be converted to a specific format, SimFuse has a complementary summary function for all of the simulated data. </li>


##Publications


##RNA-Seq Datasets
<li>We use a ENCODE MCF-7 cell line dataset (SRR521521) to generate a simulated dataset by SimFuse as an example in the paper. This dataset is available at http://sra.dnanexus.com/ in SRA format. </li>



##Setup
<li>To run SimFuse, you need to have anaconda installed in your Linux system.</li>
<li>After that, you need to get SimFuse from binstar and create the running environment by:</li>
<li> conda create -c yuxiang simfuse -n $env_name #(Recommendation: use SF with version number as the env_name.) </li>
<li>Once it is setup correctly, it will tell you how to activate the environment. For example: source activate $env_name </li>




##Reference Dataset (available in: http://smonti.bumc.bu.edu/~montilab/zoho/QueryFuse/reference_files/) need to work on this section
<li>1.	Genome.fa and build the bowtie index. For bowtie website. [Note: the fa file must use "chr*" as chromosome rather than use just 1-22, X, Y and M.] The hg19.fa I used can be also downloaded in the reference folder. (You need to build bowtie index yourself based on the fa by using "bowtie2-build".)</li>
<li>2.	Gene annotation. They can be built in bio-mart. Columns to include and following this order. For example, the one for hg19 (hg19_whole_gene_list.bed in the reference folder). In order to have better extension, I generated a 5k bp extended version on both sides of each gene. (hg19_whole_gene_list_5k_expanded.bed in the reference folder). Because the chrM is shortly, need to pay attention that not to extent the boundary over chrM's length. However, generally, chrM is not considered in fusion detection and should be filtered out.)</li>



##Input data format</li>
<li>SimFuse now takes only paired end aligned bam file from aligners as input. To convert sam to bam, samtools command "samtools view" can be used.</li>



##How to run
<li>Before running SimFuse, you need to activate the SimFuse environment by:</li>
<li> source activate $env_name #(The one used in the setup step) </li>
<li>Running SimFuse is as simple as running a single command line by just typing "SimFuse.py ????????????". </li>
<li>[Details of all parameters: please see the help information in each SimFuse version.]</li>
<li>*Note: the output directory should be different from the directory containing the input bam file.*</li>
<li>??????add the summary function related info.......</li>


The following section need more works on, should also add the parameter section.

##Output
###Final output
<li>For a specific sample, all the most important outputs are in the "results" folder and each query gene will have its own subfolder. The "whole_fusion_sum_filtered.txt" file is the key one with score and ranked result after filtering. The "whole_fusion_sum_all.txt" file contains all the detect fusions before filtering and has scores of features but no ranks. If you want to have a look at the detail graph of ranked fusion events, you can find it in the "fusion_supporting_graph_final" folder by the fusion details.</li>
<li>There are three other folders in a run: "logs" folder with all running logs; "bams" folder with shared preprocess files for all queries; "intermedias" folder with query specific intermediate files for each query. </li>


###Identification
<li>PARTNER_GENE_NAME: the name of the partner gene, which can be ambiguous.</li>
<li>PARTNER_GENE_ID: the ENSEMBL ID of the partner gene, which is unique.</li>
<li>QUERY_GENE_NAME: the name of the query gene.</li>
<li>QUERY_GENE_ID: the ENSEMBL ID of the query gene.</li>

###Evidence 
<li>SPLIT_NUM: number of splitting reads supporting this fusion event. If it is fusion event with multiple alignments, the sum of all splitting reads is divided by the number of multiple locations.</li>
<li>SPAN_NUM: number of spanning reads supporting this fusion event.</li>
<li>SUPPORT_SUM_NUM: Total supporting read number of this fusion event. It is the sum of SPLIT_NUM and SPAN_NUM.</li>
<li>SPLIT_PVAL: the value to indicate how well the splitting reads spread around the breakpoint. (The bigger, the better.) It is the p-value from the KS-test.</li>
<li>SHIFT_RANGE: the length of shared sequence at the breakpoints of the pair of partners. Shifting the breakpoint in this region will not affect the fusion sequence. Biologically, this is locations that an enzyme can bind on and the double-strand DNA break or splicing event can happen at.</li>
<li>DINUCLEOTIDE_ENTROPY: entropy of 16 possible dinucleotides around the breakpoint.</li>
<li>MULTI-ALIGN_NUM: number of possible alignment locations for this fusion reference.</li>
<li>RANK_SCORE: The sum of ranks of all previous features.</li>


###Annotation
<li>CHR_PARTNER_GENE: the chromosome of the partner gene.</li>
<li>BREAKPOINT_PARTNER_GENE: the breakpoint location of the partner gene.</li>
<li>DIRECTION_PARTNER_GENE: the chromosomal connection direction of the partner gene. ("F" means forward direction, which means it is 3' end on plus strand or 5' end on the minus strand and the detected arm is on the left (side with smaller location number than the breakpoint) of a chromosome. "R" means reverse direction, which means it is 5' end on the plus strand or 3' end on the minus strand and the detected arm is on the right (side with bigger location number than the breakpoint) of a chromosome.</li>
<li>CHR_QUERY_GENE: the chromosome of the query gene</li>
<li>BREAKPOINT_QUERY_GENE: the breakpoint location of the query gene.</li>
<li>DIRECTION_QUERY_GENE: the chromosomal connection direction of the query gene.</li>
