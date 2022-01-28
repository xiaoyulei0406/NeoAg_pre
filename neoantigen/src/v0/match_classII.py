import argparse

import pandas as pd

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser(
		description='Map to original 25mer file)')
	parser.add_argument('-i', dest='hlabind_file', help='The inputdir path')
	parser.add_argument('-o', dest='out_dir', help='The output path')
	parser.add_argument ('-r', dest='topspot_file', help='The raw topspot file')
	args = parser.parse_args()
	return args
def main():
	args = parse_arguments()
	hlabind_file=args.hlabind_file
	filename = hlabind_file.split ("/")[-1]
	filename = filename.split (".txt")[0]
	sample_ID = filename.split ("_")[0]
	outfile = args.out_dir + filename + "_summary.xlsx"
	out_dir=args.out_dir
	topspot_file=args.topspot_file
	df_pep = pd.read_table (topspot_file)
	df_mhci = pd.read_table (hlabind_file)
	df_mhci = df_mhci.merge (df_pep, left_on = df_mhci.peptide.str.extract('(\d+)', expand = False), right_on = df_pep.peptide.str.extract('(\d+)', expand = False), how = 'inner').rename(columns = {'peptide_y': 'Right_peptide'})
	df_uniq = df_mhci.drop_duplicates ()
  ###filter 25mer column  include peptide column
	##select  columns
	df_uniq = df_uniq[['allele', 'seq_num', 'start', 'end', 'length', 'identity','peptide_x', 'Methods', 'Nmethods',
					   'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
					   'nn_core_peptide', 'nn_ic50','nn_percentile_rank', 'nn_adjusted_rank',
					   'smm_core_peptide', 'smm_ic50', 'smm_percentile_rank','smm_adjusted_rank',
					   'sturniolo_core_peptide', 'sturniolo_score', 'sturniolo_percentile_rank','sturniolo_adjusted_rank',
					   'netmhciipan_core_peptide', 'netmhciipan_ic50', 'netmhciipan_percentile_rank','netmhciipan_adjusted_rank','Right_peptide']]
	df_uniq.columns = ['allele', 'seq_num', 'start', 'end', 'length', 'identity','peptide', 'Methods', 'Nmethods',
					   'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
					   'nn_core_peptide', 'nn_ic50','nn_percentile_rank', 'nn_adjusted_rank',
					   'smm_core_peptide', 'smm_ic50', 'smm_percentile_rank','smm_adjusted_rank',
					   'sturniolo_core_peptide', 'sturniolo_score', 'sturniolo_percentile_rank','sturniolo_adjusted_rank',
					   'netmhciipan_core_peptide', 'netmhciipan_ic50', 'netmhciipan_percentile_rank','netmhciipan_adjusted_rank','MT_25mer_AA']
	###filter 25mer column  include peptide column
	df_uniq = df_uniq[df_uniq.apply (lambda x: x.peptide in x.MT_25mer_AA, axis=1)]
	df_uniq = df_uniq.drop_duplicates ()
	df_uniq = df_uniq[['allele', 'seq_num', 'start', 'end', 'length', 'identity','peptide', 'Methods', 'Nmethods',
					   'comblib_core_peptide', 'comblib_ic50', 'comblib_percentile_rank', 'comblib_adjusted_rank',
					   'nn_core_peptide', 'nn_ic50','nn_percentile_rank', 'nn_adjusted_rank',
					   'smm_core_peptide', 'smm_ic50', 'smm_percentile_rank','smm_adjusted_rank',
					   'sturniolo_core_peptide', 'sturniolo_score', 'sturniolo_percentile_rank','sturniolo_adjusted_rank',
					   'netmhciipan_core_peptide', 'netmhciipan_ic50', 'netmhciipan_percentile_rank','netmhciipan_adjusted_rank']]
	out_summary= out_dir + '/' + filename + "_summary.txt"
	df_uniq.to_csv (out_summary, index=False,header=True, sep="\t")
	df_uniq.to_excel (outfile, index=False)

if __name__ == '__main__':
	main()

'''
for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do /home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/match_classII.py -i /data/cyu/topspot/05252021/output/summary/classII/${i}_length15_iedb_mhcII.txt -r /data/cyu/topspot/05252021/input/pep_fasta/pep_loc.txt -o /data/cyu/topspot/05252021/output/summary/classII/
done

for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do /home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/match_classII.py -i /data/cyu/topspot/05252021/output/summary/classII/${i}_length16_iedb_mhcII.txt -r /data/cyu/topspot/05252021/input/pep_fasta/pep_loc.txt -o /data/cyu/topspot/05252021/output/summary/classII/
done

for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do /home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/match_classII.py -i /data/cyu/topspot/05252021/output/summary/classII/${i}_length17_iedb_mhcII.txt -r /data/cyu/topspot/05252021/input/pep_fasta/pep_loc.txt -o /data/cyu/topspot/05252021/output/summary/classII/
done

for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do /home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/match_classII.py -i /data/cyu/topspot/05252021/output/summary/classII/${i}_length18_iedb_mhcII.txt -r /data/cyu/topspot/05252021/input/pep_fasta/pep_loc.txt -o /data/cyu/topspot/05252021/output/summary/classII/
done

for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do /home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/match_classII.py -i /data/cyu/topspot/05252021/output/summary/classII/${i}_length19_iedb_mhcII.txt -r /data/cyu/topspot/05252021/input/pep_fasta/pep_loc.txt -o /data/cyu/topspot/05252021/output/summary/classII/
done

for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do /home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/match_classII.py -i /data/cyu/topspot/05252021/output/summary/classII/${i}_length20_iedb_mhcII.txt -r /data/cyu/topspot/05252021/input/pep_fasta/pep_loc.txt -o /data/cyu/topspot/05252021/output/summary/classII/
done

for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do /home/ec2-user/anaconda3/bin/python /data/cyu/topspot/05252021/scripts/match_classII.py -i /data/cyu/topspot/05252021/output/summary/classII/${i}_length21_iedb_mhcII.txt -r /data/cyu/topspot/05252021/input/pep_fasta/pep_loc.txt -o /data/cyu/topspot/05252021/output/summary/classII/
done
'''

