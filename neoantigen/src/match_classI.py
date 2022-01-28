import argparse
import glob
import os
import pandas as pd
import subprocess
def get_filter(out_dir,out_tmpdir):
	for filepath in glob.glob (os.path.join (out_tmpdir, '*iedb_mhcI_summary_tmp.txt')):
		filename_tmp = filepath.split ("/")[-1]
		filename = filename_tmp.split (".")[0] + "_summary.txt"
		out_dir1=out_dir +"/filter/"
		with open ('{}/{}'.format (out_dir1, filename), 'w') as fw:
			pep = open (filepath)
			header = pep.readline ()
			fw.writelines (header)
			for l in pep.readlines ():
				line = l.strip ().split ("\t")
				if "KRAS_G12" in l:
					if (int(line[2])+int(line[3])) > 12 and int(line[2]) <13 and int(line[3]) >11:
						fw.writelines(l)
					else:
						pass
				else:
					if (int(line[2]) + int(line[3])) > 13 and int(line[2]) < 14 and int(line[3]) >12:
						fw.writelines (l)
					else:
						pass
		fw.close()
def get_exl(out_dir):
	for filepath in glob.glob (os.path.join (out_dir, '*_tmp_summary.txt')):
		filename_tmp = filepath.split ("/")[-1]
		print(filepath)
		out_name = filename_tmp.split ("_tmp")[0] + ".xlsx"
		out_file=out_dir + "/" + out_name
		df= pd.read_table (filepath, sep='\t',header=0,index_col=False)
		df_uniq=df.drop_duplicates ()
		df_uniq = df_uniq.fillna ('_')
		df_uniq.to_excel (out_file,index=False)

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser(
		description='Map to original 25mer file)')
	parser.add_argument('-o', dest='out_dir', help='The output path')
	parser.add_argument ('-r', dest='topspot_file', help='The raw topspot file')
	args = parser.parse_args()
	return args
def main():
	args = parse_arguments()
	out_dir = args.out_dir
	topspot_file = args.topspot_file
	out_tmpdir= out_dir + '/tmp/'
	if os.path.isdir (out_tmpdir):
		pass
	else:
		subprocess.call('mkdir -p ' + out_tmpdir , shell=True)

	out_filterdir = out_dir + '/filter/'
	if os.path.isdir(out_filterdir):
		pass
	else:
		subprocess.call('mkdir -p '+ out_filterdir, shell =True)
	for filepath in glob.glob (os.path.join (out_dir, '*iedb_mhcI.txt')):
		filename = filepath.split ("/")[-1]
		filename = filename.split (".txt")[0]
		df_pep = pd.read_table (topspot_file, sep='\t')
		df_mhci = pd.read_table (filepath, sep='\t')
		df_mhci = df_mhci.merge (df_pep, left_on = df_mhci.peptide.str.extract('(\d+)', expand = False), right_on = df_pep.peptide.str.extract('(\d+)', expand = False), how = 'inner').rename(columns = {'peptide_y': 'Right_peptide'})
		df_uniq = df_mhci.drop_duplicates ()
		##select  columns
		df_uniq = df_uniq[['allele', 'seq_num', 'start', 'end', 'length', 'identity','peptide_x', 'Methods', 'Nmethods',
					   'comblib_score', 'comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore', 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50', 'smm_rank',
					 'ann_ic50', 'ann_rank',
					 'netmhcpan_el_core', 'netmhcpan_el_icore', 'netmhcpan_el_score', 'netmhcpan_el_rank',
					 'pickpocket_ic50', 'pickpocket_rank',
					 'smmpmbec_ic50', 'smmpmbec_rank','Right_peptide']]
		df_uniq.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'identity','peptide', 'Methods', 'Nmethods',
					   'comblib_score', 'comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore', 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50', 'smm_rank',
					 'ann_ic50', 'ann_rank',
					 'netmhcpan_el_core', 'netmhcpan_el_icore', 'netmhcpan_el_score', 'netmhcpan_el_rank',
					 'pickpocket_ic50', 'pickpocket_rank',
					 'smmpmbec_ic50', 'smmpmbec_rank','MT_25mer_AA']
		###filter 25mer column  include peptide column
		df_uniq = df_uniq[df_uniq.apply (lambda x: x.peptide in x.MT_25mer_AA, axis=1)]
		##select new columns:
		df_uniq = df_uniq[['allele', 'seq_num', 'start', 'end', 'length', 'identity','peptide', 'Methods', 'Nmethods',
					   'comblib_score', 'comblib_rank',
					 'netmhcpan_BA_core', 'netmhcpan_BA_icore', 'netmhcpan_BA_ic50', 'netmhcpan_BA_rank',
					 'smm_ic50', 'smm_rank',
					 'ann_ic50', 'ann_rank',
					 'netmhcpan_el_core', 'netmhcpan_el_icore', 'netmhcpan_el_score', 'netmhcpan_el_rank',
					 'pickpocket_ic50', 'pickpocket_rank',
					 'smmpmbec_ic50', 'smmpmbec_rank']]
		df_uniq = df_uniq.drop_duplicates ()

		out_summary= out_tmpdir + filename + "_summary_tmp.txt"
		df_uniq.to_csv (out_summary, index=False, header=True, sep="\t")
	get_filter(out_dir,out_tmpdir)
	get_exl(out_filterdir)

if __name__ == '__main__':
	main()

'''
/home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/match_classI.py \
-o /data/cyu/topspot/05252021/output/summary/classI/ -r /data/cyu/topspot/05252021/input/pep_fasta/pep_loc.txt

'''
