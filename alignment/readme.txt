
##############################################scripts:
###########################################

1 run_coverage.py This script generates statistical summaries of QC table (total reads, mapping rate, duplication rate, mean coverage, 100X, 150X, 200X, 500X)

#Example command:
python run_coverage.py -i /data/GBM03052021/bam/ -o /data/GBM03052021/coverage.xls

ouput:
coverage.xls : quality control information


2 run_NGScheckmate.py This script performs sample swapping detection and then generate a summary table with mislabelled pairs
#Usage:
python run_NGScheckmate.py -i /data/GBM03052021/fastq/ -o /data/GBM03052021/

input:
 -i : fastq path
 -o : current work path

ouput:
1 *.pdf correlation plot of sample pairs
2 mislabel.xls (in the work path): if there is any mislabelled samples.

