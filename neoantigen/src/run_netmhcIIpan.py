import os, sys, subprocess
import argparse
import glob

netMHCIIpan='/home/ec2-user/tools/netMHCIIpan3.1/netMHCIIpan '
IEDB='/home/cyu/tools/IEDB/mhc_i/src/predict_binding.py'
methods=['ann','comblib_sidney2008','consensus','IEDB_recommended','smm','smmpmbec','pickpocket']
def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Sequenza process')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	parser.add_argument ('-a', dest='hla_file', help='hla_file')
	parser.add_argument ('-f', dest='sample_flag', help='hla flag')
	parser.add_argument ('-p', dest='pep_fasta', help='Input peptide fasta')
	args = parser.parse_args ()
	return args

def main():
	args = parse_arguments ()
	out_dir = args.out_dir

	sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed

	fout = sh_dir + args.sample_flag + '_' + 'netMHCIIpan.sh'
	sh = open (fout, 'w')
	log = log_dir + args.sample_flag + '_' + 'netMHCIIpan.log'
	subprocess.call ('mkdir -p ' + out_dir + '/classII/'+args.sample_flag +'/', shell=True)
	out_dir1=out_dir + '/classII/'+args.sample_flag +'/'
	cmd = ''
	ale = open (args.hla_file)
	cmd += 'echo starting running netMHCIIpan methods at `date` \n\n'
	pep_fasta = args.pep_fasta
	for line in ale.readlines ():
		HLA_list1 = line.split ('\n')[0]
		cmd += netMHCIIpan + ' -a ' + str (HLA_list1) + ' -f ' + pep_fasta + ' -length 15 -s'
		cmd += ' > ' + out_dir1 + args.sample_flag + '_' + str (HLA_list1) + '_netMHCIIpan_length15.log' + '\n'
		cmd += netMHCIIpan + ' -a ' + str (HLA_list1) + ' -f ' + pep_fasta + ' -length 16 -s'
		cmd += ' > ' + out_dir1 + args.sample_flag + '_' + str (HLA_list1) + '_netMHCIIpan_length16.log' + '\n'
		cmd += netMHCIIpan + ' -a ' + str (HLA_list1) + ' -f ' + pep_fasta + ' -length 17 -s'
		cmd += ' > ' + out_dir1 + args.sample_flag + '_' + str (HLA_list1) + '_netMHCIIpan_length17.log' + '\n'
		cmd += netMHCIIpan + ' -a ' + str (HLA_list1) + ' -f ' + pep_fasta + ' -length 18 -s'
		cmd += ' > ' + out_dir1 + args.sample_flag + '_' + str (HLA_list1) + '_netMHCIIpan_length18.log' + '\n'
		cmd += netMHCIIpan + ' -a ' + str (HLA_list1) + ' -f ' + pep_fasta + ' -length 19 -s'
		cmd += ' > ' + out_dir1 + args.sample_flag + '_' + str (HLA_list1) + '_netMHCIIpan_length19.log' + '\n'
		cmd += netMHCIIpan + ' -a ' + str (HLA_list1) + ' -f ' + pep_fasta + ' -length 20 -s'
		cmd += ' > ' + out_dir1 + args.sample_flag + '_' + str (HLA_list1) + '_netMHCIIpan_length20.log' + '\n'
		cmd += netMHCIIpan + ' -a ' + str (HLA_list1) + ' -f ' + pep_fasta + ' -length 21 -s'
		cmd += ' > ' + out_dir1 + args.sample_flag + '_' + str (HLA_list1) + '_netMHCIIpan_length21.log' + '\n'
		cmd += netMHCIIpan + ' -a ' + str (HLA_list1) + ' -f ' + pep_fasta + ' -length 22 -s'
		cmd += ' > ' + out_dir1 + args.sample_flag + '_' + str (HLA_list1) + '_netMHCIIpan_length22.log' + '\n'
		sh.write (cmd)

	sh.close ()

	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)

if __name__ == '__main__':
	main()

'''
for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do python /data/cyu/topspot/05252021/scripts/run_netmhcIIpan.py -o /data/cyu/topspot/05252021/output/ -a /data/cyu/topspot/05252021/input/hlalist/classII/netmhcpanII/${i}.netmhcIIpan.txt -f ${i} -p /data/cyu/topspot/05252021/input/pep_fasta/6topspot.fasta
done


'''
