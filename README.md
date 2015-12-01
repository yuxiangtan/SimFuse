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
<li>SimFuse now takes only paired end aligned bam file from aligners as input. To convert sam to bam, samtools command "samtools view" can be used. Other processing steps (i.e. PCR duplicate filtering, non-unique alignment filtering) are not required.</li>



##How to run
<li>Before running SimFuse, you need to activate the SimFuse environment by:</li>
<li> source activate $env_name #(The one used in the setup step) </li>
<li>Running SimFuse is as simple as running a single command line by just typing "SimFuse -i read_group_input -b bam_file -o working_dir -e exon_file -F SimFuse_path -G genome_ref -g log_file". </li>
<li>[Details of all parameters: please see the help information in each SimFuse version for latest updates.]</li>
<li>*Note: the output directory should be different from the directory containing the input bam file.*</li>
<li>??????add the summary function related info.......</li>

##Parameters
###The ones with no default values are required parameters, while the ones with default values are optional.
<li>-h help
<li>-i The read group file, which contains a desired split read number each row, complementary span read number will
       be calculated following read length distribution								        *[No default value, please see read_group_example.txt in the Example_files folder]</li>
<li>-b Targeted bam file location, from which the fusion background will be generated		                        *[No default value]</li>
<li>-o working/output directory [all folders should have / at the end]                                                  *[No default value]</li>
<li>-e The txt file with all exon annotations for a genome, can be generated from biomart		                *[No default value, please see HG.GRCh37.64.exons.txt in the Example_files folder]</li>
<li>-F simulator script path                                                                                            *[No default value]</li>
<li>-G tophat_genome_reference_fa - the path of the genome fa file (such as hg19.fa)                                    *[No default value]</li>
<li>-g LOG_folder, generally set to be within the working directory                                                     *[No default value]</li>
<li>-c bam filter script to get background reads with no potential fusion reads (only properly paired aligned reads)    [default value is clear_bg_filter in the SF folder]</li>
<li>-d minimum number of exons for each expression group to sample from 		                                [default value is 100*100(from -p)]</li>
<li>-E the length on each end of read you want to exclude when doing the simulation					[default value is 30 because default blat will not align read shorter than 30bp]</li>
<li>-l Read_length - length of reads		                                                                        [default value 99]</li>
<li>-M fragment size from library preparation                                                                           [default value is 0]</li>
<li>-m number of mutation allowed in each read (this is a in progress function but not supported yet)                   [default value is 0]</li>
<li>-p number of simulation gene pairs in each expression group.                                                        [default value is 100]</li>
<li>-r resume_status: check whether user want to skip finished step or start over                                       [default value 0, not resume]</li>
<li>-s standard deviation of fragment size                                                                              [default value is 0]</li>
<li>-t number of simulation rounds 											[default value is 100]</li>

##Output
###Final output stucture
<li>For each simulation round, a folder with the round number will generated under the working directory.</li>
<li>A "logs" folder is generated to record the simulation running log on both overview level and run-specific level (in folders with round number).</li>
<li>coverage_on_exons.txt shows the coverage of each exon in the input sample.</li>
<li>proper_pair_no_skip files (.bam, _1.fq and _2.fq) are the fusion-cleaned background files in different format.</li>
<li>stats.txt shows the read number of four different alignment categories of the input bam file.</li>



###In each simulation round, each split/span combination will have 14 files. Within them the two .fq files are the key files containing all simulated reads.
<li>coverage_on_exons.txt_#split_#span.bed: contains the exons randomly picked for simulation</li>
<li>coverage_on_exons.txt_#split_#span.expression_groups: shows under the input of -d, how many expression groups there are and the number of exons in each group. The header of this file is max_expression_level, group_ID, number_of_exons.</li>
<li>coverage_on_exons.txt_#split_#span.row_matrix_left and coverage_on_exons.txt_#split_#span.row_matrix_right: are just the row number of the picked exons</li>
<li>coverage_on_exons.txt_#split_#span_exon_pairs.bed: shows how the exons are paired for fusion simulation</li>
<li>coverage_on_exons.txt_#split_#span_queryID_list.txt: shows what genes can be used as query genes and with the set number of fusion partners in the simulation</li>
<li>coverage_on_exons.txt_#split_#span_queryID_partner_count.txt: shows the number of fusion partner for each query gene</li>
<li>coverage_on_exons.txt_#split_#span_ref.bed: bed file containing the reference sequence of each exon, from which fusion supporting reads are generated from.</li>
<li>coverage_on_exons.txt_#split_#span_ref.fa1 and coverage_on_exons.txt_#split_#span_ref.fa2: are the simulated reads in fa format.</li>
<li>coverage_on_exons.txt_#split_#span_ref.fq1 and coverage_on_exons.txt_#split_#span_ref.fq2: are the key outputs from SimFuse, which contains all the simulated reads. They are used to merge with back groun fq files (proper_pair_no_skip_1.fq and proper_pair_no_skip_2.fq) to generate simulated data. However, because the merged file is generally huge and causes space issues. Users need to merge them by "cat </li>
<li></li>
<li></li>
<li></li>
<li></li>
<li></li>
<li></li>
<li></li>
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
