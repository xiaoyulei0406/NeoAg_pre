###Prediction binding of peptides and MHC molecules.

######################################################################################
###MHC-I epitope prediction in IEDB-I:
######################################################################################
##Tools
ANN 4.0
NetMHCpan BA 4.0
NetMHCpan EL 4.0
SMM
SMMPMBEC
PickPocket
CombLib_Sidney2008

## Input files
Peptides with length range from 8mer-12mer

##Cutoff value
IC50 < 500nM



######################################################################################
###MHC-II epitope predictionin IEDB-I:
######################################################################################
## Tools
NetMHCIIpan BA 4.0
NetMHCIIpan EL 4.0
NN-align 2.2
SMM-align
Combinatorial library
Sturniolo

## Input files
Peptides with length range from 15mer-25mer

##Cutoff value 
IC50 < 500nM




######################################################################################
###scripts:
######################################################################################
#Go to current work directionary and then run the scripts

##1 prepare input fasta files
vep2fasta.py

#Example command:
python vep2fasta.py

#output dir: ./data/neoantigen/pep_fasta/
#output files: 
 1.*fasta--- input fasta files
 2.*pep.loc--- input fasta location information

##2 run IEDB classI 

2.1. run IEDB class I
#Usage:

python run_iedbI.py \
-i /data/GBM03052021/ \
-p /data/GBM03052021/pairs.txt \
-s 3 -e 3

2.2. run_iedbI_summary.py 
#Usage:

python run_iedbI_summary.py \
-i /data/GBM03052021/ \
-p /data/GBM03052021/pairs.txt \
-s 3 -e 3


##3 run IEDB classII

3.1 run IEDB class II
#Usage:
python run_iedbII.py \
-i /data/GBM03052021/
-p /data/GBM03052021/pairs.txt \
-s 3 -e 3

3.2 run IEDB class II summary
#Usage:
python run_iedbII_summary.py \
-i /data/GBM03052021/
-p /data/GBM03052021/pairs.txt \
-s 3 -e 3
