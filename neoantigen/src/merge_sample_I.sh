ptID=$1
inputdir=$2

awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length8.combine.txt > ${inputdir}/${ptID}_length8_iedb_mhcI.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length9.combine.txt > ${inputdir}/${ptID}_length9_iedb_mhcI.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length10.combine.txt > ${inputdir}/${ptID}_length10_iedb_mhcI.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length11.combine.txt > ${inputdir}/${ptID}_length11_iedb_mhcI.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length12.combine.txt > ${inputdir}/${ptID}_length12_iedb_mhcI.txt

