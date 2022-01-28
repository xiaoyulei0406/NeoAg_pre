
1 HLA_allele_summary.py: summarize HLAtyping results for each of three tools(HLAVBSeq + Optitype + Polysolver) and get HLA allele input for class I neoantigen prediction
Usage:
python HLA_allele_summary.py \
-i /data/GBM03052021/ \
-p /data/GBM03052021/samples.txt -s 1 -e 1

2 HLA_tool_summary.py: get a consensus result from OptiType, HLA-VBSeq and PolySolver(class I)
Usage:
python HLA_tool_summary.py \
-i /data/GBM03052021/ \
-p /data/GBM03052021/samples.txt -s 1 -e 5




