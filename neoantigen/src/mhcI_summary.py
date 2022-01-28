import sys
import matplotlib
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
import pandas.io.common

columns = ['allele', 'seq_num', 'start', 'end', 'length', 'ic50', 'rank']


def get_combine(ptID, input_dir, outDir, length, ale):

	in_ann_tocheck = input_dir + '/' + ptID + '_classI_'+ str (ale) + '_ann_length' + str (length) + '.log'
	in_smm_tocheck = input_dir + '/' + ptID+ '_classI_' + str (ale) + '_smm_length' + str (length) + '.log'
	in_smmpmbec_tocheck = input_dir + '/' + ptID + '_classI_' + str (ale) + '_smmpmbec_length' + str (length) + '.log'
	in_netBA_tocheck = input_dir + '/' + ptID+ '_classI_' + str (ale) + '_netmhcpan_ba_length' + str (length) + '.log'
	in_netEL_tocheck = input_dir + '/' + ptID + '_classI_' + str (ale) + '_netmhcpan_el_length' + str (length) + '.log'
	in_pickpocket_tocheck= input_dir + '/' + ptID + '_classI_' + str (ale) + '_pickpocket_length' + str (length) + '.log'
	in_comblib_tocheck = input_dir + '/' + ptID+ '_classI_' + str (ale) + '_comblib_sidney2008_length' + str (length) + '.log'
	if os.path.isfile (in_smm_tocheck):
		in_smm = in_smm_tocheck
		try:
			df_smm = pd.read_csv (in_smm, sep='\t')
		except pd.errors.EmptyDataError:
			#print ("The CSV file is empty")
			in_smm = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'smm_ic50',
					 'smm_rank'])
			df_smm_ic50 = in_smm
		else:
			df_smm_ic50 = df_smm.loc[df_smm['ic50'] < 500]
			df_smm_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length','peptide','smm_ic50','smm_rank']
	else:
		in_smm = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length',  'peptide', 'smm_ic50',
					 'smm_rank'])
		df_smm_ic50 = in_smm
	#print("df_smm_ic50 is :")
	#print (df_smm_ic50.shape[0])
	if os.path.isfile (in_ann_tocheck):
		in_ann = in_ann_tocheck
		try:
			df_ann = pd.read_csv (in_ann, sep='\t')
		except pd.errors.EmptyDataError:
			#print ("The CSV file is empty")
			in_ann = pd.DataFrame (
					columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'ann_ic50','ann_rank'])
			df_ann_ic50 = in_ann
		else:
			df_ann_ic50 = df_ann.loc[df_ann['ic50'] < 500]
			df_ann_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length',  'peptide','ann_ic50','ann_rank']
	else:
		in_ann = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length','peptide', 'ann_ic50','ann_rank'])
		df_ann_ic50 = in_ann
	#print("df_ann_ic50 is:")
	#print(df_ann_ic50.shape[0])
	if os.path.isfile (in_netBA_tocheck):
		in_netBA = in_netBA_tocheck

		try:
			df_netBA = pd.read_csv (in_netBA, sep='\t')
		except pd.errors.EmptyDataError:
			in_netBA = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide','netmhcpan_BA_core', 'netmhcpan_BA_icore',
						 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank'])
			df_netBA_ic50 = in_netBA
		else:
			df_netBA_ic50 = df_netBA.loc[df_netBA['ic50'] < 500]
			df_netBA_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide','netmhcpan_BA_core', 'netmhcpan_BA_icore',
						 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank']
	else:
		in_netBA = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide','netmhcpan_BA_core', 'netmhcpan_BA_icore',
						 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank'])
		df_netBA_ic50 = in_netBA
	#print ("df_netBA_ic50 is:")
	#print (df_netBA_ic50.shape[0])

	if os.path.isfile (in_comblib_tocheck):
		in_comblib = in_comblib_tocheck
		try:
			df_comblib = pd.read_csv (in_comblib, sep='\t')
		except pd.errors.EmptyDataError:
			in_comblib = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length',  'peptide','comblib_score','comblib_rank'])
			df_comblib_ic50 = in_comblib
		else:
			df_comblib_ic50 = df_comblib.loc[df_comblib['rank'] < 2]
			df_comblib_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide','comblib_score','comblib_rank']
	else:
		in_comblib = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'comblib_score','comblib_rank'])
		df_comblib_ic50 = in_comblib
	#print("df_comblib_ic50 is:")
	#print(df_comblib_ic50.shape[0])

	if os.path.isfile ((in_netEL_tocheck)):
		in_netEL = in_netEL_tocheck
		try:
			df_netEL = pd.read_csv (in_netEL, sep='\t')
		except pd.errors.EmptyDataError:
			in_netEL = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide','netmhcpan_el_core','netmhcpan_el_icore',
						 'netmhcpan_el_score', 'netmhcpan_el_rank'])
			df_netEL_ic50 = in_netEL
		else:
			df_netEL_ic50 = df_netEL.loc[df_netEL['rank'] <= 2]
			df_netEL_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide','netmhcpan_el_core','netmhcpan_el_icore',
						 'netmhcpan_el_score', 'netmhcpan_el_rank']
	else:
		in_netEL = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide','netmhcpan_el_core','netmhcpan_el_icore',
						 'netmhcpan_el_score', 'netmhcpan_el_rank'])
		df_netEL_ic50 = in_netEL
	#print("df_netEL_ic50 is:")
	#print(df_netEL_ic50.shape[0])
	#print(df_netEL_ic50)

	if os.path.isfile ((in_pickpocket_tocheck)):
		in_pickpocket = in_pickpocket_tocheck
		try:
			df_pickpocket = pd.read_csv (in_pickpocket, sep='\t')
		except pd.errors.EmptyDataError:
			in_pickpocket = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide','pickpocket_ic50', 'pickpocket_rank'])
			df_pickpocket_ic50 = in_pickpocket
		else:
			df_pickpocket_ic50 = df_pickpocket.loc[df_pickpocket['ic50'] <= 500]
			df_pickpocket_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide','pickpocket_ic50', 'pickpocket_rank']
	else:
		in_pickpocket = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide','pickpocket_ic50', 'pickpocket_rank'])
		df_pickpocket_ic50 = in_pickpocket

	if os.path.isfile ((in_smmpmbec_tocheck)):
		in_smmpmbec = in_smmpmbec_tocheck
		try:
			df_smmpmbec = pd.read_csv (in_smmpmbec, sep='\t')

		except pd.errors.EmptyDataError:
			in_smmpmbec = pd.DataFrame (columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide','smmpmbec_ic50', 'smmpmbec_rank'])
			df_smmpmbec_ic50 = in_smmpmbec
		else:
			df_smmpmbec_ic50 = df_smmpmbec.loc[df_smmpmbec['ic50'] <= 500]
			df_smmpmbec_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'peptide','smmpmbec_ic50', 'smmpmbec_rank']
	else:
		in_smmpmbec = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide','smmpmbec_ic50', 'smmpmbec_rank'])
		df_smmpmbec_ic50 = in_smmpmbec
	############################################################

	##column name
	sample_smm = 'SMM_align'
	sample_ann = 'ANN_align'
	sample_netBA = 'NetMHCpan_BA'
	sample_netEL = 'NetMHCpan_EL'
	#	sample_iedb = 'iedb'
	#	sample_consensus = 'consensus'
	sample_comblib = 'comblib'
	sample_pickpocket = 'pickpocket'
	sample_smmpmbec ="smmpmbec"

	if (df_netBA_ic50.shape[0]==0 and df_comblib_ic50.shape[0]==0):
		df_join =pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_score','comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore','netmhcpan_BA_ic50', 'netmhcpan_BA_rank'])
	else:
		df_join = pd.merge (df_netBA_ic50, df_comblib_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0]==0 and df_smm_ic50.shape[0]==0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_score','comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore','netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50','smm_rank'])
	else:
		df_join = pd.merge (df_join, df_smm_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0]==0 and df_ann_ic50.shape[0]==0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_score','comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore','netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50','smm_rank',
					 'ann_ic50','ann_rank'])
	else:
		df_join = pd.merge (df_join, df_ann_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0] == 0 and df_netEL_ic50.shape[0] == 0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_score','comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore','netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50','smm_rank',
					 'ann_ic50','ann_rank',
					 'netmhcpan_el_core','netmhcpan_el_icore','netmhcpan_el_score', 'netmhcpan_el_rank'])
	else:
		df_join = pd.merge (df_join, df_netEL_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0] == 0 and df_pickpocket_ic50.shape[0] == 0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_score', 'comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore', 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50', 'smm_rank',
					 'ann_ic50', 'ann_rank',
					 'netmhcpan_el_core', 'netmhcpan_el_icore', 'netmhcpan_el_score', 'netmhcpan_el_rank',
					 'pickpocket_ic50', 'pickpocket_rank'])
	else:
		df_join = pd.merge (df_join, df_pickpocket_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0] == 0 and df_smmpmbec_ic50.shape[0] == 0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_score', 'comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore', 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50', 'smm_rank',
					 'ann_ic50', 'ann_rank',
					 'netmhcpan_el_core', 'netmhcpan_el_icore', 'netmhcpan_el_score', 'netmhcpan_el_rank',
					 'pickpocket_ic50', 'pickpocket_rank',
					 'smmpmbec_ic50', 'smmpmbec_rank'])
	else:
		df_join = pd.merge (df_join, df_smmpmbec_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)

	#	print (df_join.shape)
	cols = []
	if (df_join.shape[0] != 0):
		for indx, row in df_join.iterrows ():
			temp = []
			if (df_netBA_ic50[
			((df_netBA_ic50['allele'] == row['allele']) & (df_netBA_ic50['seq_num'] == row['seq_num']) &
			 (df_netBA_ic50['start'] == row['start']) & (df_netBA_ic50['end'] == row['end']) &
			 (df_netBA_ic50['length'] == row['length']) & (df_netBA_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('NetMHCpan_BA')
			else:
				temp.append (np.nan)
			if (df_comblib_ic50[
			((df_comblib_ic50['allele'] == row['allele']) & (df_comblib_ic50['seq_num'] == row['seq_num']) &
			 (df_comblib_ic50['start'] == row['start']) & (df_comblib_ic50['end'] == row['end']) &
			 (df_comblib_ic50['length'] == row['length']) &
			 (df_comblib_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('comblib')
			else:
				temp.append (np.nan)
			if (df_smm_ic50[
			((df_smm_ic50['allele'] == row['allele']) & (df_smm_ic50['seq_num'] == row['seq_num']) &
			 (df_smm_ic50['start'] == row['start']) & (df_smm_ic50['end'] == row['end']) &
			 (df_smm_ic50['length'] == row['length']) &
			 (df_smm_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('SMM_align')
			else:
				temp.append (np.nan)
			if (df_ann_ic50[((df_ann_ic50['allele'] == row['allele']) & (df_ann_ic50['seq_num'] == row['seq_num']) &
						(df_ann_ic50['start'] == row['start']) & (df_ann_ic50['end'] == row['end']) &
						(df_ann_ic50['length'] == row['length']) &
						(df_ann_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('ANN_align')
			else:
				temp.append (np.nan)
			if (df_netEL_ic50[
			((df_netEL_ic50['allele'] == row['allele']) & (df_netEL_ic50['seq_num'] == row['seq_num']) &
			 (df_netEL_ic50['start'] == row['start']) & (df_netEL_ic50['end'] == row['end']) &
			 (df_netEL_ic50['length'] == row['length']) &
			 (df_netEL_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('NetMHCpan_EL')
			else:
				temp.append (np.nan)
			if (df_pickpocket_ic50[
			((df_pickpocket_ic50['allele'] == row['allele']) & (df_pickpocket_ic50['seq_num'] == row['seq_num']) &
			 (df_pickpocket_ic50['start'] == row['start']) & (df_pickpocket_ic50['end'] == row['end']) &
			 (df_pickpocket_ic50['length'] == row['length']) &
			 (df_pickpocket_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('pickpocket')
			else:
				temp.append (np.nan)
			if (df_smmpmbec_ic50[
			((df_smmpmbec_ic50['allele'] == row['allele']) & (df_smmpmbec_ic50['seq_num'] == row['seq_num']) &
			 (df_smmpmbec_ic50['start'] == row['start']) & (df_smmpmbec_ic50['end'] == row['end']) &
			 (df_smmpmbec_ic50['length'] == row['length']) &
			 (df_smmpmbec_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('smmpmbec')
			else:
				temp.append (np.nan)

			cols.append (temp)
		df_join[sample_netBA] = [x[0] for x in cols]
		df_join[sample_comblib] = [x[1] for x in cols]
		df_join[sample_smm] = [x[2] for x in cols]
		df_join[sample_ann] = [x[3] for x in cols]
		df_join[sample_netEL] = [x[4] for x in cols]
		df_join[sample_pickpocket] = [x[5] for x in cols]
		df_join[sample_smmpmbec] = [x[6] for x in cols]

	### add all callers into a new column
		df_join['Methods'] = df_join[['NetMHCpan_BA', 'comblib', 'SMM_align', 'ANN_align', 'NetMHCpan_EL','pickpocket','smmpmbec']].apply (
		lambda x: ';'.join (x.dropna ()), axis=1)
	### add number of callers into a new column
		df_join['Nmethods'] = df_join['Methods'].str.count (';') + 1
	# print(df_join)
		df_join = df_join.fillna ('_')
	# Get a series containing maximum value of each row
		df_join = df_join[['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'Methods', 'Nmethods', 'comblib_score', 'comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore', 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50', 'smm_rank',
					 'ann_ic50', 'ann_rank',
					 'netmhcpan_el_core', 'netmhcpan_el_icore', 'netmhcpan_el_score', 'netmhcpan_el_rank',
					 'pickpocket_ic50', 'pickpocket_rank',
					 'smmpmbec_ic50', 'smmpmbec_rank']]
		df_join_uniq = df_join.drop_duplicates ()
		df_join_uniq.to_csv (outDir + '/' + ptID + '_' + str (ale) + '_length' + str (length) + '.combine.txt', index=False, header=True, sep='\t')
		return (df_join_uniq)
	else:
		df_join_uniq=[]


def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Combine all the annotation results')
	parser.add_argument ('-i', dest='input_dir', help='The inputdir path')
	parser.add_argument ('-p', dest='pt_file', help='Input patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')
	parser.add_argument ('-o', dest='outdir', help='The output path')
	args = parser.parse_args ()
	return args


def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	out_dir = args.outdir
	pt_file = args.pt_file
	start = args.start
	end = args.end
	if os.path.isdir(out_dir):
		pass
	else:
		subprocess.call ('mkdir -p ' + out_dir, shell=True)
	pts = open (pt_file)  # note there is a header
	lns = pts.readlines ()
	input_classI = input_dir + '/data/neoantigen/IEDB_i/'

	for i in range ((int) (start), (int) (end)+1):
		tmp = lns[i].strip ('\n').split (',')
		ptID = tmp[1]   # tumor patient ID
		print (ptID)
		hla_file = input_dir + '/data/neoantigen/inputHLA/' + ptID + '.hla.IEDBI.txt'
		hlas = open (hla_file)
		hla_list1 = hlas.readlines ()
		#hla_list = str(hla_list1).split ('\n')[0]
		#print (len (hla_list))
		for line in hla_list1:
			hla_list=str(line).split ('\n')[0]
			print (hla_list)
			ale = hla_list.replace ('*', '')
			print(ale)
			get_combine (ptID, input_classI, out_dir, 8, ale)
			get_combine (ptID, input_classI, out_dir, 9, ale)
			get_combine (ptID, input_classI, out_dir, 10, ale)
			get_combine (ptID, input_classI, out_dir, 11, ale)
			get_combine (ptID, input_classI, out_dir, 12, ale)
	#		print(df)
	#		df_join = pd.DataFrame (np.concatenate ([df_SNP.values, df_INDEL.values,df_MNP_nn.values]),columns=df_INDEL.columns)
	#		df_join.to_csv (outdir + '/' + ptID + '.combine.all.txt', index=False, header=True, sep='\t')

	pts.close ()


if __name__ == '__main__':
	main ()

'''

/home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/mhcI_summary.py \
-i /data/cyu/topspot/05252021/ \
-p /data/cyu/topspot/05252021/input/sample.txt \
-s 1 -e 17 -o /data/cyu/topspot/05252021/output/summary/classI/sep/


'''
