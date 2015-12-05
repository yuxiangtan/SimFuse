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




##Reference Dataset (available in: http://smonti.bumc.bu.edu/~montilab/zoho/SimFuse/reference_files/)
<li>1.	Genome.fa. For bowtie website. [Note: the fa file must use "chr*" as chromosome rather than use just 1-22, X, Y and M.] The hg19.fa I used can be also downloaded in the reference folder. </li>
<li>2.	Exon annotation. They can be built in bio-mart. Columns to include and following this order. For example, the one for hg19 (HG.GRCh37.64.exons.txt in the reference folder). </li>



##Input Data Format</li>
<li>SimFuse now takes only paired end aligned bam file from aligners as input. To convert sam to bam, samtools command "samtools view" can be used. Other processing steps (i.e. PCR duplicate filtering, non-unique alignment filtering) are not required.</li>



##How to Run
<li>Before running SimFuse, you need to activate the SimFuse environment by:</li>
<li> source activate $env_name #(The one used in the setup step) </li>
<li>Running SimFuse is as simple as running a single command line by just typing "SimFuse -i read_group_input -b bam_file -o working_dir -e exon_file -F SimFuse_path -G genome_ref -g log_file". </li>
<li>[Details of all parameters: please see the help information in each SimFuse version for latest updates.]</li>
<li>*Note: the output directory should be different from the directory containing the input bam file.*</li>
<li>Scripts for summarization and comparison on results of fusion detection algorithms are in the SimFuse_summary_scripts folder. However,because different fusion detection algorithms have their own output format, there is no uniform adaptor to do this yet. SimFuse only supports format conversion for deFuse, TophatFusion and QueryFuse output. This section needs further improvement on building an automatic wrapper.</li>

##Parameters
###The ones with no default values are required parameters, while the ones with default values are optional.
<li>-h help
<li>-i The read group file, which contains a desired split read number each row, complementary span read number will
       be calculated following read length distribution								        *[No default value, please see read_group_example.txt in the Example_files folder]</li>
<li>-b Targeted bam file location, from which the fusion background will be generated		                        *[No default value]</li>
<li>-o working/output directory [all folders should have / at the end]                                                  *[No default value]</li>
<li>-e The txt file with all exon annotations for a genome, can be generated from biomart		                *[No default value, please see Reference Dataset 2 HG.GRCh37.64.exons.txt in the Example_files folder]</li>
<li>-F simulator script path (generally, it is ~/anaconda_envs/$env_name/lib/SimFuse/                                   *[No default value]</li>
<li>-G tophat_genome_reference_fa - the path of the genome fa file (such as hg19.fa)                                    *[No default value, please see Reference Dataset 1.]</li>
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
###Final Output Stucture
<li>For each simulation round, a folder with the round number will generated under the working directory.</li>
<li>A "logs" folder is generated to record the simulation running log on both overview level and run-specific level (in folders with round number).</li>
<li>coverage_on_exons.txt shows the coverage of each exon in the input sample.</li>
<li>proper_pair_no_skip files (.bam, _1.fq and _2.fq) are the fusion-cleaned background files in different format.</li>
<li>stats.txt shows the read number of four different alignment categories of the input bam file.</li>



###In each simulation round, each split/span combination will have 14 files. Within them the two .fq files are the key files containing all simulated reads.
<li>coverage_on_exons.txt_#split_#span_ref.fq1 and coverage_on_exons.txt_#split_#span_ref.fq2: are the key outputs from SimFuse, which contains all the simulated reads. They are used to merge with background fq files (proper_pair_no_skip_1.fq and proper_pair_no_skip_2.fq) to generate simulated data. However, because the merged file is generally huge and causes space issues. As a result, users need to merge them by: cat "coverage_on_exons.txt"$ID"_ref.fq1" $no_skip_bg_path"/proper_pair_no_skip_1.fq" > "coverage_on_exons.txt"$ID"_ref_merged_def.1.fastq" before running fusion detectors and then remove the merged files.</li>
###Other supplementary files that users may want to know.
<li>coverage_on_exons.txt_#split_#span.bed: contains the exons randomly picked for simulation</li>
<li>coverage_on_exons.txt_#split_#span.expression_groups: shows under the input of -d, how many expression groups there are and the number of exons in each group. The header of this file is max_expression_level, group_ID, number_of_exons.</li>
<li>coverage_on_exons.txt_#split_#span.row_matrix_left and coverage_on_exons.txt_#split_#span.row_matrix_right: are just the row number of the picked exons</li>
<li>coverage_on_exons.txt_#split_#span_exon_pairs.bed: shows how the exons are paired for fusion simulation</li>
<li>coverage_on_exons.txt_#split_#span_queryID_list.txt: shows what genes can be used as query genes and with the set number of fusion partners in the simulation</li>
<li>coverage_on_exons.txt_#split_#span_queryID_partner_count.txt: shows the number of fusion partner for each query gene</li>
<li>coverage_on_exons.txt_#split_#span_ref.bed: bed file containing the reference sequence of each exon, from which fusion supporting reads are generated from.</li>
<li>coverage_on_exons.txt_#split_#span_ref.fa1 and coverage_on_exons.txt_#split_#span_ref.fa2: are the simulated reads in fa format.</li>
<li>coverage_on_exons.txt_100_26_ref_bp.txt: shows the breakpoint locations of each pair of fusion partners.</li>


###Summary Function 
<li>By using the provided summary scripts in the SimFuse_summary_scripts folder, statistical summary tables (recall and precision) and complementary barplot figures (see 3method_detected_fusion_recall.png and 3method_detected_fusion_precision.png in Example_files) will be generated.</li>
<li>In the tables (see QF_summary_detected_fusion_recall.txt and QF_summary_detected_fusion_precision.txt in Example_files), as shown in the hearder, the first column indicates the split_span combination; the second column is the recall/pricision; the third column is its standard deviation; the fourth column indicates the 95% confidence interval, which is used to generate error bars in the plots.</li>


