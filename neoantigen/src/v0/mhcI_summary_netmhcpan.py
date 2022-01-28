'''

v0:change filter option to ic50 <500
v1 :change filter option to rank <2 or ic50 <500


DRB1 in sturniolo

'''
import sys
import matplotlib
import pandas as pd
import pandas.io.common
import numpy as np
import argparse
import subprocess
import re
import os

columns = ['allele', 'seq_num', 'start', 'end', 'length', 'core_peptide', 'peptide', 'ic50', 'percentile_rank',
		   'adjusted_rank']




def get_combine(ptID, input_dir, out_dir, length, mer_dir,flag):
	in_net_tocheck = input_dir + 'classI_binder_' + ptID + '_netMHCpanI_' + 'length' + str (length) + '.txt'
	colnames = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'ic50', 'rank']

	if os.path.isfile (in_net_tocheck):
		try:
			in_net = in_net_tocheck
			df_net = pd.read_csv (in_net, sep='\t')
			df_net_ic50 = df_net.loc[df_net['ic50'] < 500]
			df_net_ic50.columns = ['allele', 'start', 'length', 'peptide', 'net_Core', 'net_Icore', 'net_ic50',
								   'net_rank']
			df_net_ic50 = df_net_ic50[
				['allele', 'start', 'length', 'peptide', 'net_Core', 'net_Icore', 'net_ic50', 'net_rank']]
		#			print(df_net_ic50)
		except pandas.io.common.EmptyDataError:
			in_net = pd.DataFrame (
				columns=['allele', 'start', 'length', 'peptide', 'net_Core', 'net_Icore', 'net_ic50', 'net_rank'])
			df_net_ic50 = in_net
	else:
		in_net = pd.DataFrame (
			columns=['allele', 'start', 'length', 'peptide', 'net_Core', 'net_Icore', 'net_ic50', 'net_rank'])
		df_net_ic50 = in_net

	df_net_ic50 = df_net_ic50.drop_duplicates ()
	df_net_ic50_1 = df_net_ic50.loc[df_net_ic50['start'] < 14]
	df_net_ic50_uniq = df_net_ic50_1.drop_duplicates ()
	df_net_ic50_uniq['end'] = df_net_ic50_uniq['start'] + df_net_ic50_uniq['length']
	out_raw = out_dir + '/' + ptID + '.length'+str(length)+'.netmhcpan.classI.raw.txt'
	df_net_ic50_uniq.to_csv (out_raw, index=False, header=True, sep='\t')
	in_25mer = mer_dir + ptID + '.'+ flag + '.uniq.25mer.txt'
	get_checkmut(in_25mer,out_dir,out_raw,length,ptID)


# print(df_new)

def get_checkmut(in_25mer,out_dir,out_raw,length,ptID):
	df_25mer=open(in_25mer,'r').readlines()
	df_raw=open(out_raw,'r').readlines()
	out_mut=out_dir + '/' + ptID + '.length'+str(length)+'.netmhcpan.classI.txt'
	fw=open(out_mut,'w')
	fw.writelines('HLA-allele'+"\t"+'start'+"\t"+'length'+"\t"+ 'peptide'+"\t"+ "identity" + "\t" +
				'net_Core'+"\t"+'net_Icore'+"\t"+ 'net_ic50'+"\t"+ 'net_rank'+"\n")

	fd=[]
	for line in df_raw:
		ar=line.strip ().split ("\t")
		for line_mer in df_25mer:
			line_mer1 = line_mer.strip ().split ("\t")
			pepseq = line_mer1[20]
			loc_mut=line_mer1[30]
			geneName=line_mer1[1]
			if line_mer1[8]=="_":
				aaName=str(line_mer1[7]).split(":")[-1]
			elif line_mer1[8]!="_":
				aaName=line_mer1[8]
			if str (ar[3]) in line_mer:
				print(ar[3])
				print(loc_mut)
				if (int(loc_mut) >= int(ar[1])) and ((int(ar[8])-1) >= int(loc_mut)):
					fw.writelines(str(ar[0])+"\t"+str(ar[1])+"\t"+str(ar[2])+"\t"+str(ar[3])+"\t"+geneName+"_"+aaName+ "\t"+
					str(ar[4])+"\t"+str(ar[5])+"\t"+ str(ar[6])+"\t"+str(ar[7])+"\n")
	fw.close ()
	get_uniq (out_mut,out_dir,length,ptID)
def get_uniq(out_mut,out_dir,length,ptID):
    out_uniq = out_dir + '/' + ptID + '.length' + str (length) + '.netmhcpan.uniq.classI.txt'
    df = pd.read_csv (out_mut, sep='\t', na_filter=False)
    df_uniq = df.drop_duplicates ()
    df_uniq.to_csv (out_uniq, index=False, header=True, sep='\t')



def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Combine all mhcI results')
	parser.add_argument ('-i', dest='input_dir', help='The inputdir path')
	parser.add_argument ('-m', dest='mer_dir', help='The inputdir path')
	parser.add_argument ('-p', dest='pt_file', help='Input patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')
	parser.add_argument ('-o', dest='outdir', help='The output path')
	parser.add_argument ('-f', dest='flag', help='The input flag: SNP, MNP, INDEL')
	args = parser.parse_args ()
	return args


def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	mer_dir = args.mer_dir
	out_dir = args.outdir
	pt_file = args.pt_file
	start = args.start
	end = args.end
	flag = args.flag

	pts = open (pt_file)  # note there is a header
	lns = pts.readlines ()

	for i in range ((int) (start), (int) (end) + 1):
		print (lns[i])
		tmp = lns[i].strip ('\n').split (',')

		ptID = tmp[1]  # tumor patient ID

		for length in range (8, 13):
			get_combine(ptID, input_dir, out_dir, length, mer_dir,flag)

	pts.close ()


if __name__ == '__main__':
	main ()

'''

/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/neoantigen/mhcI_summary_netmhcpan.py \
-i /data/cyu/gbm/data/neoantigen/classI/ \
-m /data/cyu/gbm/data/VEP/25mer/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 1 -e 3 -o /data/cyu/gbm/data/neoantigen/classI_summary_netmhcpan/ \
-f SNP

/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/neoantigen/mhcI_summary_netmhcpan.py \
-i /data/cyu/gbm/data/neoantigen/classI/ \
-m /data/cyu/gbm/data/VEP/25mer/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 5 -e 5 -o /data/cyu/gbm/data/neoantigen/classI_summary_netmhcpan/ \
-f SNP

'''
