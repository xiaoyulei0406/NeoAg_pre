import os, sys, subprocess
import argparse
import glob

netMHCpan='/home/ec2-user/tools/netMHCpan-4.0/netMHCpan'
IEDB='/home/cyu/tools/IEDB/mhc_i/src/predict_binding.py'
methods=['ann','comblib_sidney2008','consensus','IEDB_recommended','smm','smmpmbec','pickpocket']
def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Sequenza process')
	parser.add_argument ('-a', dest='hla_file', help='hla_file')
	parser.add_argument ('-f', dest='hla_flag', help='hla flag')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	parser.add_argument ('-p', dest='pep_fasta', help='Input peptide fasta')
	args = parser.parse_args ()
	return args

def main():
	args = parse_arguments ()
	out_dir = args.out_dir
	hla_file = args.hla_file
	hla_flag = args.hla_flag

	sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed

	fout = sh_dir + 'mhcI_netMHCpan_'+ hla_flag + '.sh'
	sh = open (fout, 'w')
	log = log_dir + 'mhcI_netMHCpan_'+ hla_flag + '.log'
	subprocess.call ('mkdir -p ' + out_dir + '/'+ hla_flag +'/', shell=True)
	cmd = ''
	HLA_file = hla_file
	hlas = open (HLA_file)
	pep_fasta=args.pep_fasta
	cmd += 'echo starting running netMHCpan at $date \n\n'
	out_dir1=out_dir + '/'+ hla_flag +'/'
	for line in hlas.readlines ():
		HLA_list = line.split ('\n')[0]
		#pep_fasta = input_dir + '/pep_fasta/' + 'topSPOTs.fasta'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 8 -s '
		cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + HLA_list + '_NetMHCpan_length8.log' + '\n'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 9 -s '
		cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + HLA_list + '_NetMHCpan_length9.log' + '\n'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 10 -s '
		cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + HLA_list + '_NetMHCpan_length10.log' + '\n'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 11 -s '
		cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + HLA_list + '_NetMHCpan_length11.log' + '\n'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 12 -s '
		cmd += ' > ' + out_dir1 + '/' + hla_flag + '_classI_' + HLA_list + '_NetMHCpan_length12.log' + '\n'
		cmd += '\n\n'
	sh.write (cmd)
	sh.close ()

	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)

if __name__ == '__main__':
	main()

'''
python /data/cyu/topspot/05252021/scripts/run_netmhcpan.py \
-o /data/cyu/topspot/05252021/output/classI/ \
-a /data/cyu/topspot/05252021/input/hlalist/classI/BCD11.hla.mhcpanI.txt \
-f BCD11 \
-p /data/cyu/topspot/05252021/input/pep_fasta/6topspot.fasta

for i in {BCD4,BCD5,BCD6,BCD8,BCD9,BCD11,BCD17,BCD19,BCD20,BCD26,BCD27,BCD28,BCD30,BCD33,BCD35,BCD36,BCD37}
do python /data/cyu/topspot/05252021/scripts/run_netmhcpan.py -o /data/cyu/topspot/05252021/output/classI/ -a /data/cyu/topspot/05252021/input/hlalist/classI/${i}.hla.mhcpanI.txt -f ${i} -p /data/cyu/topspot/05252021/input/pep_fasta/6topspot.fasta
done

'''
