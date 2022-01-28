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


def get_netmhcpan(ptID, input_dir):
	for length in range (8, 13):

		in_net_tocheck = input_dir + 'classI_' + ptID + '_netMHCpanI_length' + str (length) + '.log'
		output = input_dir + 'classI_binder_' + ptID + '_netMHCpanI_length' + str (length) + '.txt'
		fw = open (output, 'w')
		fw.writelines (
			"allele" + "\t" + "pos" + "\t" + "length" + "\t" + "Peptide" + "\t" + "Core" + "\t" + "Icore" + "\t" + "ic50" + "\t" + "Rank" + "\n")
		if os.path.isfile (in_net_tocheck):
			in_net = in_net_tocheck
			df_net = open (in_net, 'r')
			for line in df_net.readlines ():
				# print(line)
				if line.startswith (" "):
					# print(line)
					line = re.sub (' +', '\t', line)
					line1 = line.strip ("\n").split ("\t")
					if str (line1[14]) != "%Rank" and float (line1[1]) < 14:
						if float (line1[14]) < 2 or float (line1[13]) < 500 :
							fw.writelines (str (line1[2]) + "\t" + str (line1[1]) + "\t" + str (length) + "\t" + str (
								line1[3]) + "\t" + str (line1[4]) + "\t" + str (line1[10]) + "\t" + str (
								line1[13]) + "\t" + str (line1[14]) + "\n")
		else:
			pass
		fw.close ()


def get_combine(ptID, input_dir, out_dir, length, mer_dir,flag):
	in_smm_tocheck = input_dir + 'classI_' + ptID + '_smm_' + 'length' + str (length) + '.log'
	in_smmpmbec_tocheck = input_dir + 'classI_' + ptID + '_smmpmbec_' + 'length' + str (length) + '.log'
	in_net_tocheck = input_dir + 'classI_binder_' + ptID + '_netMHCpanI_' + 'length' + str (length) + '.txt'
	in_comblib_tocheck = input_dir + 'classI_' + ptID + '_comblib_sidney2008_' + 'length' + str (length) + '.log'
	in_ann_tocheck = input_dir + 'classI_' + ptID + '_ann_' + 'length' + str (length) + '.log'
	in_pickpocket_tocheck = input_dir + 'classI_' + ptID + '_pickpocket_' + 'length' + str (length) + '.log'
	colnames = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'ic50', 'rank']
	if os.path.isfile (in_smm_tocheck):
		try:
			in_smm = in_smm_tocheck
			df_smm = pd.read_csv (in_smm, sep='\t', na_filter=False)
			#		df_mut_smm = df_smm[['allele','seq_num','start','end','length','core_peptide','peptide','ic50','percentile_rank','adjusted_rank']]
			df_smm_ic50 = df_smm.loc[(df_smm['rank'] < 2) | (df_smm['ic50'] < 500 )]
			df_smm_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'smm_ic50', 'smm_rank']
			df_smm_ic50 = df_smm_ic50[['allele', 'start', 'length', 'peptide', 'smm_ic50', 'smm_rank']]
		#			print(df_smm_ic50)
		except pandas.io.common.EmptyDataError:
			in_smm = pd.DataFrame (columns=['allele', 'start', 'length', 'peptide', 'smm_ic50', 'smm_rank'])
			df_smm_ic50 = in_smm
	else:
		in_smm = pd.DataFrame (columns=['allele', 'start', 'length', 'peptide', 'smm_ic50', 'smm_rank'])
		df_smm_ic50 = in_smm

	df_smm_ic50 = df_smm_ic50.drop_duplicates ()

	if os.path.isfile (in_ann_tocheck):
		try:
			in_ann = in_ann_tocheck
			df_ann = pd.read_csv (in_ann, sep='\t')
			df_ann_ic50 = df_ann.loc[(df_ann['rank'] < 2) | (df_ann['ic50'] < 500)]
			df_ann_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'ann_ic50', 'ann_rank']
			df_ann_ic50 = df_ann_ic50[['allele', 'start', 'length', 'peptide', 'ann_ic50', 'ann_rank']]
		#			print(df_ann_ic50)
		except pandas.io.common.EmptyDataError:
			in_ann = pd.DataFrame (columns=['allele', 'start', 'length', 'peptide', 'ann_ic50', 'ann_rank'])
			df_ann_ic50 = in_ann
	else:
		in_ann = pd.DataFrame (columns=['allele', 'start', 'length', 'peptide', 'ann_ic50', 'ann_rank'])
		df_ann_ic50 = in_ann

	df_ann_ic50 = df_ann_ic50.drop_duplicates ()

	if os.path.isfile (in_net_tocheck):
		try:
			in_net = in_net_tocheck
			df_net = pd.read_csv (in_net, sep='\t')
			df_net_ic50 = df_net.loc[(df_net['Rank'] < 2 ) | (df_net['Rank'] < 500)]
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

	if os.path.isfile (in_comblib_tocheck):
		try:
			in_comblib = in_comblib_tocheck
			df_comblib = pd.read_csv (in_comblib, sep='\t')
			df_comblib_ic50 = df_comblib.loc[(df_comblib['rank'] < 2) | (df_comblib['ic50'] < 500)]
			df_comblib_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'comblib_ic50',
									   'comblib_rank']
			df_comblib_ic50 = df_comblib_ic50[['allele', 'start', 'length', 'peptide', 'comblib_ic50', 'comblib_rank']]
		#			print(df_comblib_ic50)
		except pandas.io.common.EmptyDataError:
			in_comblib = pd.DataFrame (columns=['allele', 'start', 'length', 'peptide', 'comblib_ic50', 'comblib_rank'])
			df_comblib_ic50 = in_comblib
	else:
		in_comblib = pd.DataFrame (columns=['allele', 'start', 'length', 'peptide', 'comblib_ic50', 'comblib_rank'])
		df_comblib_ic50 = in_comblib

	df_comblib_ic50 = df_comblib_ic50.drop_duplicates ()

	if os.path.isfile ((in_smmpmbec_tocheck)):
		try:
			in_smmpmbec = in_smmpmbec_tocheck
			df_smmpmbec = pd.read_csv (in_smmpmbec, sep='\t')
			df_smmpmbec_ic50 = df_smmpmbec.loc[(df_smmpmbec['rank'] < 2) | (df_smmpmbec['ic50']<500)]
			df_smmpmbec_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'smmpmbec_ic50',
										'smmpmbec_rank']
			df_smmpmbec_ic50 = df_smmpmbec_ic50[
				['allele', 'start', 'length', 'peptide', 'smmpmbec_ic50', 'smmpmbec_rank']]
		#			print(df_smmpmbec_ic50)
		except pandas.io.common.EmptyDataError:
			in_smmpmbec = pd.DataFrame (
				columns=['allele', 'start', 'length', 'peptide', 'smmpmbec_ic50', 'smmpmbec_rank'])
			df_smmpmbec_ic50 = in_smmpmbec
	else:
		in_smmpmbec = pd.DataFrame (columns=['allele', 'start', 'length', 'peptide', 'smmpmbec_ic50', 'smmpmbec_rank'])
		df_smmpmbec_ic50 = in_smmpmbec

	df_smmpmbec_ic50 = df_smmpmbec_ic50.drop_duplicates ()
	if os.path.isfile ((in_pickpocket_tocheck)):
		try:
			in_pickpocket = in_pickpocket_tocheck
			df_pickpocket = pd.read_csv (in_pickpocket, sep='\t')
			df_pickpocket_ic50 = df_pickpocket.loc[(df_pickpocket['rank'] < 2) | (df_pickpocket['ic50'] < 500)]
			df_pickpocket_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'pickpocket_ic50',
										  'pickpocket_rank']
			df_pickpocket_ic50 = df_pickpocket_ic50[
				['allele', 'start', 'length', 'peptide', 'pickpocket_ic50', 'pickpocket_rank']]
		#			print(df_pickpocket_ic50)
		except pandas.io.common.EmptyDataError:
			in_pickpocket = pd.DataFrame (
				columns=['allele', 'start', 'length', 'peptide', 'pickpocket_ic50', 'pickpocket_rank'])
			df_pickpocket_ic50 = in_pickpocket
	else:
		in_pickpocket = pd.DataFrame (
			columns=['allele', 'start', 'length', 'peptide', 'pickpocket_ic50', 'pickpocket_rank'])
		df_pickpocket_ic50 = in_pickpocket

	df_pickpocket_ic50 = df_pickpocket_ic50.drop_duplicates ()
	############################################################

	if ((df_net_ic50.shape[0] == 0) and (df_comblib_ic50.shape[0] == 0) and
			(df_smm_ic50.shape[0] == 0) and (df_ann_ic50.shape[0] == 0) and
			(df_smmpmbec_ic50.shape[0] == 0) and (df_pickpocket_ic50.shape[0] == 0)):
		print ("all are null")
		pass
	else:
		mhci_join (df_net_ic50, df_comblib_ic50, df_smm_ic50, df_ann_ic50, df_smmpmbec_ic50, df_pickpocket_ic50,
				   out_dir, ptID, length, mer_dir,flag)


def mhci_join(df_net_ic50, df_comblib_ic50, df_smm_ic50, df_ann_ic50, df_smmpmbec_ic50, df_pickpocket_ic50, out_dir,
			  ptID, length, mer_dir, flag):
	df_join = pd.merge (df_net_ic50, df_comblib_ic50,
						on=['allele', 'start', 'length', 'peptide'],
						how='outer', indicator=False)
	df_join = pd.merge (df_join, df_smm_ic50,
						on=['allele', 'start', 'length', 'peptide'],
						how='outer', indicator=False)
	df_join = pd.merge (df_join, df_ann_ic50,
						on=['allele', 'start', 'length', 'peptide'],
						how='outer', indicator=False)
	df_join = pd.merge (df_join, df_smmpmbec_ic50,
						on=['allele', 'start', 'length', 'peptide'],
						how='outer', indicator=False)
	df_join = pd.merge (df_join, df_pickpocket_ic50,
						on=['allele', 'start', 'length', 'peptide'],
						how='outer', indicator=False)
	#	print (df_join.shape)
	cols = []
	for indx, row in df_join.iterrows ():
		temp = []
		if (df_net_ic50[
			((df_net_ic50['allele'] == row['allele']) &
			 (df_net_ic50['start'] == row['start']) &
			 (df_net_ic50['length'] == row['length']) & (df_net_ic50['peptide'] == row['peptide']))].shape[0] == 1):
			temp.append ('netMHCpanI')
		else:
			temp.append (np.nan)
		if (df_comblib_ic50[
			((df_comblib_ic50['allele'] == row['allele']) &
			 (df_comblib_ic50['start'] == row['start']) &
			 (df_comblib_ic50['length'] == row['length']) &
			 (df_comblib_ic50['peptide'] == row['peptide']))].shape[0] == 1):
			temp.append ('comblib_sidney2008')
		else:
			temp.append (np.nan)
		if (df_smm_ic50[
			((df_smm_ic50['allele'] == row['allele']) &
			 (df_smm_ic50['start'] == row['start']) &
			 (df_smm_ic50['length'] == row['length']) &
			 (df_smm_ic50['peptide'] == row['peptide']))].shape[0] == 1):
			temp.append ('smm')
		else:
			temp.append (np.nan)
		if (df_ann_ic50[((df_ann_ic50['allele'] == row['allele']) &
						 (df_ann_ic50['start'] == row['start']) &
						 (df_ann_ic50['length'] == row['length']) &
						 (df_ann_ic50['peptide'] == row['peptide']))].shape[0] == 1):
			temp.append ('ann')
		else:
			temp.append (np.nan)
		if (df_smmpmbec_ic50[
			((df_smmpmbec_ic50['allele'] == row['allele']) &
			 (df_smmpmbec_ic50['start'] == row['start']) &
			 (df_smmpmbec_ic50['length'] == row['length']) &
			 (df_smmpmbec_ic50['peptide'] == row['peptide']))].shape[0] == 1):
			temp.append ('smmpmbec')
		else:
			temp.append (np.nan)
		if (df_pickpocket_ic50[
			((df_pickpocket_ic50['allele'] == row['allele']) &
			 (df_pickpocket_ic50['start'] == row['start']) &
			 (df_pickpocket_ic50['length'] == row['length']) &
			 (df_pickpocket_ic50['peptide'] == row['peptide']))].shape[0] == 1):
			temp.append ('pickpocket')
		else:
			temp.append (np.nan)
		cols.append (temp)

	##column name
	sample_smm = 'smm'
	sample_ann = 'ann'
	sample_net = 'netMHCpanI'
	sample_smmpmbec = 'smmpmbec'
	sample_comblib = 'comblib_sidney2008'
	sample_pickpocket = 'pickpocket'

	df_join[sample_net] = [x[0] for x in cols]
	df_join[sample_comblib] = [x[1] for x in cols]
	df_join[sample_smm] = [x[2] for x in cols]
	df_join[sample_ann] = [x[3] for x in cols]
	df_join[sample_smmpmbec] = [x[4] for x in cols]
	df_join[sample_pickpocket] = [x[5] for x in cols]
	### add all callers into a new column
	df_join['Methods'] = df_join[['netMHCpanI', 'comblib_sidney2008', 'smm', 'ann', 'smmpmbec', 'pickpocket']].apply (
		lambda x: ';'.join (x.dropna ()), axis=1)
	### add number of callers into a new column
	df_join['Nmethods'] = df_join['Methods'].str.count (';') + 1
	df_join = df_join.fillna ('_')
	# Get a series containing maximum value of each row

	df_join = df_join[['allele', 'start', 'length', 'peptide', 'Methods', 'Nmethods',
					   'ann_ic50', 'ann_rank', 'comblib_ic50', 'comblib_rank',
					   'net_Core', 'net_Icore', 'net_ic50', 'net_rank',
					   'smm_ic50', 'smm_rank', 'smmpmbec_ic50', 'smmpmbec_rank',
					   'pickpocket_ic50', 'pickpocket_rank']]

	df_join1 = df_join.loc[df_join['start'] < 14]
	df_join_uniq = df_join1.drop_duplicates ()
	df_join_uniq['end'] = df_join_uniq['start'] + df_join_uniq['length']
	out_raw = out_dir + '/' + ptID + '.length'+str(length)+'.consensus.classI.raw.txt'
	df_join_uniq.to_csv (out_raw, index=False, header=True, sep='\t')
	in_25mer = mer_dir + ptID + '.'+ flag + '.uniq.25mer.txt'
	get_checkmut(in_25mer,out_dir,out_raw,length,ptID)


# print(df_new)

def get_checkmut(in_25mer,out_dir,out_raw,length,ptID):
	df_25mer=open(in_25mer,'r').readlines()
	df_raw=open(out_raw,'r').readlines()
	out_mut=out_dir + '/' + ptID + '.length'+str(length)+'.consensus.classI.txt'
	fw=open(out_mut,'w')
	fw.writelines('HLA-allele'+"\t"+'start'+"\t"+'length'+"\t"+ 'peptide'+"\t"+ 'Methods'+"\t"+ 'Nmethods'+"\t"+
					"identity"  "\t" +
					'ann_ic50'+"\t"+ 'ann_rank'+"\t"+ 'comblib_ic50'+"\t"+ 'comblib_rank'+"\t"+
					'net_Core'+"\t"+'net_Icore'+"\t"+ 'net_ic50'+"\t"+ 'net_rank'+"\t"+
					'smm_ic50'+"\t"+'smm_rank'+"\t"+ 'smmpmbec_ic50'+"\t"+'smmpmbec_rank'+"\t"+
					'pickpocket_ic50'+"\t"+'pickpocket_rank'+"\n")

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
				if (int(loc_mut) >= int(ar[1])) and ((int(ar[20])-1) >= int(loc_mut)):
					fw.writelines(str(ar[0])+"\t"+str(ar[1])+"\t"+str(ar[2])+"\t"+str(ar[3])+"\t"+
					str(ar[4])+"\t"+str(ar[5])+"\t"+geneName+"_"+aaName+ "\t"+ str(ar[6])+"\t"+str(ar[7])+"\t"+
					str(ar[8])+"\t"+str(ar[9])+"\t"+str(ar[10])+"\t"+str(ar[11])+"\t"+
					str(ar[12])+"\t"+str(ar[13])+"\t"+str(ar[14])+"\t"+str(ar[15])+"\t"+
					str(ar[16])+"\t"+str(ar[17])+"\t"+str(ar[18])+"\t"+str(ar[19])+"\n")
	fw.close ()
	get_uniq (out_mut,out_dir,length,ptID)
def get_uniq(out_mut,out_dir,length,ptID):
    out_uniq = out_dir + '/' + ptID + '.length' + str (length) + '.consensus.uniq.classI.txt'
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
		get_netmhcpan (ptID, input_dir)
		for length in range (8, 13):
			get_combine(ptID, input_dir, out_dir, length, mer_dir,flag)

	pts.close ()


if __name__ == '__main__':
	main ()

'''

/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/mhcI_summary_IEDB.py \
-i /data/cyu/gbm/data/neoantigen/classI/ \
-m /data/cyu/gbm/data/VEP/25mer/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 1 -e 3 -o /data/cyu/gbm/data/neoantigen/classI_summary_v1/ \
-f SNP


'''
