import sys
import matplotlib
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
import pandas.io.common

def get_combine(ptID, input_dir, outDir, length, ale):
	in_smm_tocheck = input_dir + '/' + ptID+ '_IEDB_smm_' + str (ale) + '_length_' + str (length) + '.log'
	in_nn_tocheck = input_dir + '/' + ptID+ '_IEDB_nn_' + str (ale) + '_length_' + str (length) + '.log'
	in_net_el_tocheck = input_dir + '/' + ptID+ '_IEDB_NetMHCIIpan_el_' + str (ale) + '_length_' + str (length) + '.log'
	in_net_ba_tocheck = input_dir + '/' + ptID + '_IEDB_NetMHCIIpan_ba_' + str (ale) + '_length_' + str (length) + '.log'
	in_comblib_tocheck = input_dir + '/' + ptID+ '_IEDB_comblib_' + str (ale) + '_length_' + str (length) + '.log'
	in_sturniolo_tocheck = input_dir + '/' + ptID+ '_IEDB_sturniolo_' + str (ale) + '_length_' + str (length) + '.log'
	if os.path.isfile (in_smm_tocheck):
		in_smm = in_smm_tocheck
		try:
			df_smm = pd.read_csv (in_smm, sep='\t')
		except pd.errors.EmptyDataError:
			#print ("The CSV file is empty")
			in_smm = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'smm_core_peptide', 'peptide', 'smm_ic50',
					 'smm_percentile_rank', 'smm_adjusted_rank'])
			df_smm_ic50 = in_smm
		else:
			df_smm_ic50 = df_smm.loc[df_smm['ic50'] < 500]
			df_smm_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'smm_core_peptide', 'peptide',
								   'smm_ic50',
								   'smm_percentile_rank', 'smm_adjusted_rank']
	else:
		in_smm = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'smm_core_peptide', 'peptide', 'smm_ic50',
					 'smm_percentile_rank', 'smm_adjusted_rank'])
		df_smm_ic50 = in_smm
	#print("df_smm_ic50 is :")
	#print (df_smm_ic50.shape[0])
	if os.path.isfile (in_nn_tocheck):
		in_nn = in_nn_tocheck
		try:
			df_nn = pd.read_csv (in_nn, sep='\t')
		except pd.errors.EmptyDataError:
			#print ("The CSV file is empty")
			in_nn = pd.DataFrame (
					columns=['allele', 'seq_num', 'start', 'end', 'length', 'nn_core_peptide', 'peptide', 'nn_ic50',
						 'nn_percentile_rank', 'nn_adjusted_rank'])
			df_nn_ic50 = in_nn
		else:
			df_nn_ic50 = df_nn.loc[df_nn['ic50'] < 500]
			df_nn_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'nn_core_peptide', 'peptide',
								  'nn_ic50',
								  'nn_percentile_rank', 'nn_adjusted_rank']
	else:
		in_nn = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'nn_core_peptide', 'peptide', 'nn_ic50',
					 'nn_percentile_rank', 'nn_adjusted_rank'])
		df_nn_ic50 = in_nn
	#print("df_nn_ic50 is:")
	#print(df_nn_ic50.shape[0])
	if os.path.isfile (in_net_el_tocheck):
		in_net = in_net_el_tocheck

		try:
			df_net = pd.read_csv (in_net, sep='\t')
		except pd.errors.EmptyDataError:
			in_net = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'netmhciipan_el_core_peptide', 'peptide',
						 'netmhciipan_el_score', 'netmhciipan_el_percentile_rank', 'netmhciipan_el_adjusted_rank'])
			df_net_ic50 = in_net
		else:
			df_net_ic50 = df_net.loc[df_net['percentile_rank'] <=10]
			df_net_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'netmhciipan_el_core_peptide', 'peptide',
							   'netmhciipan_el_score', 'netmhciipan_el_percentile_rank', 'netmhciipan_el_adjusted_rank']
	else:
		in_net = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'netmhciipan_el_core_peptide', 'peptide',
					 'netmhciipan_el_score', 'netmhciipan_el_percentile_rank', 'netmhciipan_el_adjusted_rank'])
		df_net_ic50 = in_net
	#print ("df_net_ic50 is:")
	#print (df_net_ic50.shape[0])

	if os.path.isfile (in_net_ba_tocheck):
		in_net_ba = in_net_ba_tocheck

		try:
			df_net_ba = pd.read_csv (in_net_ba, sep='\t')
		except pd.errors.EmptyDataError:
			in_net_ba = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'netmhciipan_ba_core_peptide', 'peptide',
						 'netmhciipan_ba_ic50', 'netmhciipan_ba_percentile_rank', 'netmhciipan_ba_adjusted_rank'])
			df_net_ba_ic50 = in_net_ba
		else:
			df_net_ba_ic50 = df_net_ba.loc[df_net_ba['ic50'] < 500]
			df_net_ba_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'netmhciipan_ba_core_peptide', 'peptide',
							   'netmhciipan_ba_ic50', 'netmhciipan_ba_percentile_rank', 'netmhciipan_ba_adjusted_rank']
	else:
		in_net_ba = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'netmhciipan_ba_core_peptide', 'peptide',
					 'netmhciipan_ba_ic50', 'netmhciipan_ba_percentile_rank', 'netmhciipan_ba_adjusted_rank'])
		df_net_ba_ic50 = in_net_ba

	if os.path.isfile (in_comblib_tocheck):
		in_comblib = in_comblib_tocheck
		try:
			df_comblib = pd.read_csv (in_comblib, sep='\t')
		except pd.errors.EmptyDataError:
			in_comblib = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'comblib_core_peptide', 'peptide',
						 'comblib_ic50',
						 'comblib_percentile_rank', 'comblib_adjusted_rank'])
			df_comblib_ic50 = in_comblib
		else:
			df_comblib_ic50 = df_comblib.loc[df_comblib['ic50'] < 500]
			df_comblib_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'comblib_core_peptide', 'peptide',
								   'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank']
	else:
		in_comblib = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'comblib_core_peptide', 'peptide', 'comblib_ic50',
					 'comblib_percentile_rank', 'comblib_adjusted_rank'])
		df_comblib_ic50 = in_comblib
	#print("df_comblib_ic50 is:")
	#print(df_comblib_ic50.shape[0])

	if os.path.isfile ((in_sturniolo_tocheck)):
		in_sturniolo = in_sturniolo_tocheck
		try:
			df_sturniolo = pd.read_csv (in_sturniolo, sep='\t')

		except pd.errors.EmptyDataError:
			in_sturniolo = pd.DataFrame (
				columns=['allele', 'seq_num', 'start', 'end', 'length', 'sturniolo_core_peptide', 'peptide',
						 'sturniolo_score', 'sturniolo_percentile_rank', 'sturniolo_adjusted_rank'])
			df_sturniolo_ic50 = in_sturniolo
		else:
			df_sturniolo_ic50 = df_sturniolo.loc[df_sturniolo['percentile_rank'] <= 10]
			df_sturniolo_ic50.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'sturniolo_core_peptide', 'peptide',
									 'sturniolo_score', 'sturniolo_percentile_rank', 'sturniolo_adjusted_rank']
	else:
		in_sturniolo = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'sturniolo_core_peptide', 'peptide',
					 'sturniolo_score', 'sturniolo_percentile_rank', 'sturniolo_adjusted_rank'])
		df_sturniolo_ic50 = in_sturniolo
	#print("df_sturniolo_ic50 is:")
	#print(df_sturniolo_ic50.shape[0])
	#print(df_sturniolo_ic50)
	############################################################

	##column name
	sample_smm = 'SMM_align'
	sample_nn = 'NN_align'
	sample_net = 'NetMHCIIpan_el'
	sample_net_ba='NetMHCIIpan_ba'
	#	sample_iedb = 'iedb'
	#	sample_consensus = 'consensus'
	sample_comblib = 'comblib'
	sample_sturniolo = 'sturniolo'

	if (df_net_ic50.shape[0]==0 and df_comblib_ic50.shape[0]==0):
		df_join =pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
					 'netmhciipan_el_core_peptide','netmhciipan_el_score', 'netmhciipan_el_percentile_rank', 'netmhciipan_el_adjusted_rank'])
	else:
		df_join = pd.merge (df_net_ic50, df_comblib_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0]==0 and df_smm_ic50.shape[0]==0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
					 'netmhciipan_el_core_peptide','netmhciipan_el_score', 'netmhciipan_el_percentile_rank', 'netmhciipan_el_adjusted_rank',
					 'smm_core_peptide', 'smm_ic50', 'smm_percentile_rank','smm_adjusted_rank'])
	else:
		df_join = pd.merge (df_join, df_smm_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0]==0 and df_nn_ic50.shape[0]==0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
					 'netmhciipan_el_core_peptide','netmhciipan_el_score', 'netmhciipan_el_percentile_rank', 'netmhciipan_el_adjusted_rank',
					 'smm_core_peptide', 'smm_ic50', 'smm_percentile_rank','smm_adjusted_rank',
					 'nn_core_peptide', 'nn_ic50','nn_percentile_rank', 'nn_adjusted_rank'])
	else:
		df_join = pd.merge (df_join, df_nn_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0] == 0 and df_sturniolo_ic50.shape[0] == 0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
					 'netmhciipan_el_core_peptide', 'netmhciipan_el_score', 'netmhciipan_el_percentile_rank','netmhciipan_el_adjusted_rank',
					 'smm_core_peptide', 'smm_ic50', 'smm_percentile_rank','smm_adjusted_rank',
					 'nn_core_peptide', 'nn_ic50','nn_percentile_rank', 'nn_adjusted_rank',
					 'sturniolo_core_peptide', 'sturniolo_score', 'sturniolo_percentile_rank','sturniolo_adjusted_rank'])
	else:
		df_join = pd.merge (df_join, df_sturniolo_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	if (df_join.shape[0] == 0 and df_net_ba_ic50.shape[0] == 0):
		df_join = pd.DataFrame (
			columns=['allele', 'seq_num', 'start', 'end', 'length', 'peptide',
					 'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
					 'netmhciipan_el_core_peptide', 'netmhciipan_el_score', 'netmhciipan_el_percentile_rank','netmhciipan_el_adjusted_rank',
					 'smm_core_peptide', 'smm_ic50', 'smm_percentile_rank', 'smm_adjusted_rank',
					 'nn_core_peptide', 'nn_ic50', 'nn_percentile_rank', 'nn_adjusted_rank',
					 'sturniolo_core_peptide', 'sturniolo_score', 'sturniolo_percentile_rank','sturniolo_adjusted_rank',
					 'netmhciipan_ba_core_peptide','netmhciipan_ba_ic50', 'netmhciipan_ba_percentile_rank', 'netmhciipan_ba_adjusted_rank'])
	else:
		df_join = pd.merge (df_join, df_net_ba_ic50,
						on=['allele', 'seq_num', 'start', 'end', 'length', 'peptide'],
						how='outer', indicator=False)
	#	print (df_join.shape)
	cols = []
	if (df_join.shape[0] != 0):
		for indx, row in df_join.iterrows ():
			temp = []
			if (df_net_ic50[
			((df_net_ic50['allele'] == row['allele']) & (df_net_ic50['seq_num'] == row['seq_num']) &
			 (df_net_ic50['start'] == row['start']) & (df_net_ic50['end'] == row['end']) &
			 (df_net_ic50['length'] == row['length']) & (df_net_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('NetMHCIIpan_el')
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
			if (df_nn_ic50[((df_nn_ic50['allele'] == row['allele']) & (df_nn_ic50['seq_num'] == row['seq_num']) &
						(df_nn_ic50['start'] == row['start']) & (df_nn_ic50['end'] == row['end']) &
						(df_nn_ic50['length'] == row['length']) &
						(df_nn_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('NN_align')
			else:
				temp.append (np.nan)
			if (df_sturniolo_ic50[
			((df_sturniolo_ic50['allele'] == row['allele']) & (df_sturniolo_ic50['seq_num'] == row['seq_num']) &
			 (df_sturniolo_ic50['start'] == row['start']) & (df_sturniolo_ic50['end'] == row['end']) &
			 (df_sturniolo_ic50['length'] == row['length']) &
			 (df_sturniolo_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('sturniolo')
			else:
				temp.append (np.nan)
			if (df_net_ba_ic50[
			((df_net_ba_ic50['allele'] == row['allele']) & (df_net_ba_ic50['seq_num'] == row['seq_num']) &
			 (df_net_ba_ic50['start'] == row['start']) & (df_net_ba_ic50['end'] == row['end']) &
			 (df_net_ba_ic50['length'] == row['length']) &
			 (df_net_ba_ic50['peptide'] == row['peptide']))].shape[0] == 1):
				temp.append ('NetMHCIIpan_ba')
			else:
				temp.append (np.nan)

			cols.append (temp)
		df_join[sample_net] = [x[0] for x in cols]
		df_join[sample_comblib] = [x[1] for x in cols]
		df_join[sample_smm] = [x[2] for x in cols]
		df_join[sample_nn] = [x[3] for x in cols]
		df_join[sample_sturniolo] = [x[4] for x in cols]
		df_join[sample_net_ba] = [x[5] for x in cols]

	### add all callers into a new column
		df_join['Methods'] = df_join[['NetMHCIIpan_el', 'comblib', 'SMM_align', 'NN_align', 'sturniolo','NetMHCIIpan_ba']].apply (
		lambda x: ';'.join (x.dropna ()), axis=1)
	### add number of callers into a new column
		df_join['Nmethods'] = df_join['Methods'].str.count (';') + 1
	# print(df_join)
		df_join = df_join.fillna ('_')
	# Get a series containing maximum value of each row
		df_join = df_join[['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'Methods', 'Nmethods',
						   'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
						   'nn_core_peptide', 'nn_ic50','nn_percentile_rank', 'nn_adjusted_rank',
						   'smm_core_peptide', 'smm_ic50', 'smm_percentile_rank','smm_adjusted_rank',
						   'sturniolo_core_peptide', 'sturniolo_score', 'sturniolo_percentile_rank','sturniolo_adjusted_rank',
						   'netmhciipan_ba_core_peptide','netmhciipan_ba_ic50', 'netmhciipan_ba_percentile_rank', 'netmhciipan_ba_adjusted_rank',
						   'netmhciipan_el_core_peptide', 'netmhciipan_el_score', 'netmhciipan_el_percentile_rank','netmhciipan_el_adjusted_rank'
						   ]]
		df_join_uniq = df_join.drop_duplicates ()
		df_join_uniq.to_csv (outDir + '/' + ptID + '_' + str (ale) +  '_length' + str (length) + '.combine.txt', index=False, header=True, sep='\t')
		return (df_join_uniq)
	else:
		df_join_uniq=[]
def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Combine all the IEDB results')
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
	outdir = args.outdir
	pt_file = args.pt_file
	start = args.start
	end = args.end

	pts = open (pt_file)  # note there is a header
	lns = pts.readlines ()
	input_classII = input_dir + '/data/neoantigen/IEDB_ii/'

	for i in range ((int) (start), (int) (end)+1):
		#print (lns[i])
		tmp = lns[i].strip ('\n').split (',')
		ptID = tmp[1]  # tumor patient ID
		hla_file = input_dir + '/data/neoantigen/inputHLA/' + ptID + '.hla.IEDBII.txt'
		hlas = open (hla_file)
		hla_list1 = hlas.readlines ()
		hla_list = hla_list1[0].strip('\n').split (',')
		#print (len (hla_list))
		for i in range (0, len (hla_list)):
			#print (hla_list[i])
			if '/' in str (hla_list[i]):
				ale2 = hla_list[i].replace ('/', '')
				ale1 = ale2.replace ('*', '')
				ale = ale1.replace (':', '')
			else:
				ale1 = hla_list[i].replace ('*', '')
				ale = ale1.replace (':', '')
			get_combine (ptID, input_classII, outdir, 15, ale)
			get_combine (ptID, input_classII, outdir, 16, ale)
			get_combine (ptID, input_classII, outdir, 17, ale)
			get_combine (ptID, input_classII, outdir, 18, ale)
			get_combine (ptID, input_classII, outdir, 19, ale)
			get_combine (ptID, input_classII, outdir, 20, ale)
			get_combine (ptID, input_classII, outdir, 21, ale)
			get_combine (ptID, input_classII, outdir, 22, ale)
			get_combine (ptID, input_classII, outdir, 23, ale)
			get_combine (ptID, input_classII, outdir, 24, ale)
			get_combine (ptID, input_classII, outdir, 25, ale)
	#		print(df)
	#		df_join = pd.DataFrame (np.concatenate ([df_SNP.values, df_INDEL.values,df_MNP_nn.values]),columns=df_INDEL.columns)
	#		df_join.to_csv (outdir + '/' + ptID + '.combine.all.txt', index=False, header=True, sep='\t')

	pts.close ()


if __name__ == '__main__':
	main ()

'''

/home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/mhcII_summary.py \
-i /data/cyu/topspot/05252021/ \
-p /data/cyu/topspot/05252021/input/sample.txt \
-s 1 -e 17 -o /data/cyu/topspot/05252021/output/summary/classII/


'''
