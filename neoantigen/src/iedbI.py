import os, sys, subprocess
import argparse
import glob

netMHCpan='/home/ec2-user/tools/netMHCpan-4.0/netMHCpan'
IEDB='/home/cyu/tools/IEDB/mhc_i/src/predict_binding.py'
methods=['ann','comblib_sidney2008','smm','smmpmbec','pickpocket','netmhcpan_ba','netmhcpan_el']
def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Sequenza process')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	parser.add_argument ('-a', dest='hla_file', help='hla_file')
	parser.add_argument ('-f', dest='hla_flag', help='hla flag')
	parser.add_argument ('-p', dest='pep_fasta', help='Input peptide fasta')
	args = parser.parse_args ()
	return args

def main():
	args = parse_arguments ()
	out_dir = args.out_dir
	HLA_file = args.hla_file
	hla_flag = args.hla_flag

	sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed

	if os.path.isdir (sh_dir) and os.path.isdir (log_dir):
		pass
	else:
		subprocess.call ('mkdir -p ' + out_dir + '/{sh,log}/', shell=True)
	fout = sh_dir + 'mhcI_iedb_'+ hla_flag + '.sh'
	sh = open (fout, 'w')
	log = log_dir + 'mhcI_iedb_'+ hla_flag + '.log'
	subprocess.call ('mkdir -p ' + out_dir + '/', shell=True)

	cmd = ''
	hlas1 = open (HLA_file)
	cmd += 'echo starting running IEDB methods at `date` \n\n'
	pep_fasta = args.pep_fasta
	out_dir1 = out_dir + '/'
	for line in hlas1.readlines ():
		HLA_list1 = line.split ('\n')[0]
		HLA_tmp=str(HLA_list1).split("*")
		print(str(HLA_tmp[0])+str(HLA_tmp[1]))
		for method in methods:
			cmd += IEDB + ' ' + method + ' ' + str (HLA_list1) + ' ' + ' 8 ' + pep_fasta
			cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + str(HLA_tmp[0])+str(HLA_tmp[1]) + '_' + method + '_length8.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLA_list1) + ' ' + ' 9 ' + pep_fasta
			cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + str(HLA_tmp[0])+str(HLA_tmp[1]) + '_' + method + '_length9.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLA_list1) + ' ' + ' 10 ' + pep_fasta
			cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + str(HLA_tmp[0])+str(HLA_tmp[1]) + '_' + method + '_length10.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLA_list1) + ' ' + ' 11 ' + pep_fasta
			cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + str(HLA_tmp[0])+str(HLA_tmp[1]) + '_' + method + '_length11.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLA_list1) + ' ' + ' 12 ' + pep_fasta
			cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + str(HLA_tmp[0])+str(HLA_tmp[1]) + '_' + method + '_length12.log' + '\n'

	sh.write (cmd)

	sh.close ()

	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)

if __name__ == '__main__':
	main()

