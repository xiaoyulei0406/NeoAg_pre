# Neoantigen-Pred

## Introduction
Some codes summarizing data for neoanitgen prediciton

## Input files
Paired tumor-normal WES fastq files (optional: RNA-seq).

## Main analysis
Alignment

Fastq files were aligned to human reference genome by BWA-MEM. Picard and GATK toolkits were used for mark PCR duplicates, calculate base quality score relcalibration and local alignment. Bedtools and qualimaps were used for quality control include calculating coverage, mapping rate, duplication rate. Sample swapping was able to detected in this process.

Somatic mutation calling and annotation

Mutect2, Mutect, MuSE, Strelka and Varscan2 were used to call mutations. SNPs and Indels were kept if at least three of above tools were found.
Annovar and VEP were used to annotate amino acid changes, changed peptide sequence.
Mutation filtering
SNPs and Indels were combined if at least 3 variant reads in the tumor sample less than 2 variant reads in the normal samples. Variants with at least 10 (wild and variant type) reads in both normal and tumor samples were kept. Variants with allele frequency more than 5% were kept. Variants in the exonic region were kept. Variants annotated with synonymous SNV were removed.
Amino acid peptide sequences in wild and mutated type
Amino acid peptide sequences carrying the mutations were prepared with 25-mer length.

HLA typing 

Number of HLA-alleles with Genomic Sequences in this pipeline

|A    |B    |C    |DPA1 |DPB1 |DQA1 |DQB1 |DRB1 |DRB3 |DRB4 |DRB5 |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|3510 |4033 |3895 |128  |562  |178  |415  |308  |30   |22   |10   |


The HLA typing was identified with combination of Optitye, POLYSOLVER and HLA-VBSeq. If all these tools identified same allele in both tumor and normal samples, the alleles were taken with high confidence. If not, two of these tools identified same allele in normal samples, the alleles were taken. If all these tools not be able to have consistent results, all the potential results were listed as low confidence.

HLA-LOH event was tested based on McGranahan N, et al. Tumor purity, tumor ploidy, BAF and logR were used to calculate major and minor allele copy number. A copy number < 0.5, was identified as LOH.


Neoantigen prediction

IEDB methods for both class I and class II were used. 
Binding affinity <500 were used to filtered.

Tumor purity and clonal analysis

Sequenza was used to predict tumor purity and ploidy.

Pyclone was used to predict clonal cell population. Notes: Pyclone needs specical inputs which includes mutationID, reference read counts, mutated read counts, number of normal copy number, number of major copy number, number of minor copy number, genotyp, Gene Name, Sample ID, Copy number Abberation. Examples are like below:

|mutation_id	|ref_counts	  |var_counts	  |normal_cn	  |major_cn	    |minor_cn	    |genotype	    |Gene_Name	  |Sample_ID	                 |Copy_number_Abberation |
|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|----------------------------|-----------------------|
|17:7577121	  |30	          |150	        |2	          |2            |0	          |BB	          |TP53	        |TCGA-02-0003-01A-01D-1490-08|Copy Number Neutral    |
|3:178936082	|100	        |50           |2	          |1            |1	          |AB	          |PIK3CA	      |TCGA-02-0003-01A-01D-1490-08|Normal                 |
|7:150500865	|300	        |10           |2	          |2            |1	          |ABB          |TMEM176A	    |TCGA-02-0003-01A-01D-1490-08|Amplification          |
|7:1097861	  |217	        |26	          |2	          |1	          |0	          |B	          |GPR146	      |TCGA-02-0003-01A-01D-1490-08|Deletion               |




