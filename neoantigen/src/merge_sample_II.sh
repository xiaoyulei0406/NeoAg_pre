ptID=$1
inputdir=$2
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length15.combine.txt > ${inputdir}/${ptID}_length15_iedb_mhcII.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length16.combine.txt > ${inputdir}/${ptID}_length16_iedb_mhcII.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length17.combine.txt > ${inputdir}/${ptID}_length17_iedb_mhcII.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length18.combine.txt > ${inputdir}/${ptID}_length18_iedb_mhcII.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length19.combine.txt > ${inputdir}/${ptID}_length19_iedb_mhcII.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length20.combine.txt > ${inputdir}/${ptID}_length20_iedb_mhcII.txt
awk 'FNR>1 || NR==1' ${inputdir}/${ptID}*_length21.combine.txt > ${inputdir}/${ptID}_length21_iedb_mhcII.txt
