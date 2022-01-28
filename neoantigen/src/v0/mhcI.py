'''
seperate run hlas

'''
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
	parser.add_argument ('-i', dest='input_dir', help='The inputdir path')
	parser.add_argument ('-p', dest='pt_file', help='Input patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	args = parser.parse_args ()
	return args

def get_HLAfile(tmp,hla_path):
	for file in glob.iglob (hla_path):  # generator, search immediate subdirectories
		filename_tmp=file.split("/")[-1]
		sample_name=filename_tmp.split(".")[0]
		if tmp == sample_name:
			return file
		else:
			pass

def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	out_dir = args.out_dir
	pt_file = args.pt_file
	start = args.start
	end = args.end

	sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed

	fout = sh_dir + 'mhcI_v1_' + start + '_' + end + '.sh'
	sh = open (fout, 'w')
	pts = open (pt_file)  # note there is a header
	log = log_dir + 'mhcI_v1_' + start + '_' + end + '.log'
	lns = pts.readlines ()
	subprocess.call ('mkdir -p ' + out_dir + '/data/neoantigen/classI_v1/', shell=True)
	cmd = ''

	for i in range ((int) (start), (int) (end) + 1):
		print (lns[i])
		tmp = lns[i].strip('\n').split(',')
		sample_id = tmp[1]  # tumor patient ID
		print(sample_id)
		HLA_file= input_dir + '/inputHLA/' + sample_id + '.hla.mhcpanI.txt'
		hlas = open (HLA_file)  # note there is a header
		line = hlas.readlines ()
		HLA_list = line[0].split ('\n')[0]
		HLA_file1 = input_dir + '/inputHLA/' + sample_id + '.hla.IEDBI.txt'
		hlas1 = open (HLA_file1)  # note there is a header
		line1 = hlas1.readlines ()
		HLA_list1 = line1[0].split ('\n')[0]

		HLAale=HLA_list1.split(",")
		print(HLAale[0])

		pep_fasta = input_dir + '/pep_fasta/' + sample_id + '.classI.fasta'
		cmd += '\n\n'

		cmd += 'echo starting running netMHCpan at $(date) \n\n'
		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 8 -s '
		cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_netMHCpanI_length8.log' + '\n'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 9 -s '
		cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_netMHCpanI_length9.log' + '\n'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 10 -s '
		cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_netMHCpanI_length10.log' + '\n'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 11 -s '
		cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_netMHCpanI_length11.log' + '\n'
		cmd += '\n\n'

		cmd += netMHCpan + ' -a ' + HLA_list + ' -f ' + pep_fasta + ' -BA -l 12 -s '
		cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_netMHCpanI_length12.log' + '\n'
		cmd += '\n\n'

		cmd += 'echo starting running IEDB methods at `date` \n\n'
		for method in methods:
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[0]) + ' ' + ' 8 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[0])+ '_length8.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[1]) + ' ' + ' 8 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[1])+ '_length8.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[2]) + ' ' + ' 8 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[2])+ '_length8.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[3]) + ' ' + ' 8 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[3])+ '_length8.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[4]) + ' ' + ' 8 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[4])+ '_length8.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[5]) + ' ' + ' 8 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[5])+ '_length8.log' + '\n'

			cmd += IEDB + ' ' + method + ' ' + str (HLAale[0]) + ' ' + ' 9 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[0])+ '_length9.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[1]) + ' ' + ' 9 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[1])+ '_length9.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[2]) + ' ' + ' 9 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[2])+ '_length9.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[3]) + ' ' + ' 9 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[3])+ '_length9.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[4]) + ' ' + ' 9 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[4])+ '_length9.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[5]) + ' ' + ' 9 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[5])+ '_length9.log' + '\n'

			cmd += IEDB + ' ' + method + ' ' + str (HLAale[0]) + ' ' + ' 10 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[0])+ '_length10.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[1]) + ' ' + ' 10 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[1])+ '_length10.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[2]) + ' ' + ' 10 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[2])+ '_length10.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[3]) + ' ' + ' 10 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[3])+ '_length10.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[4]) + ' ' + ' 10 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[4])+ '_length10.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[5]) + ' ' + ' 10 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[5])+ '_length10.log' + '\n'

			cmd += IEDB + ' ' + method + ' ' + str (HLAale[0]) + ' ' + ' 11 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[0])+ '_length11.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[1]) + ' ' + ' 11 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[1])+ '_length11.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[2]) + ' ' + ' 11 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[2])+ '_length11.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[3]) + ' ' + ' 11 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[3])+ '_length11.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[4]) + ' ' + ' 11 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[4])+ '_length11.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[5]) + ' ' + ' 11 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[5])+ '_length11.log' + '\n'

			cmd += IEDB + ' ' + method + ' ' + str (HLAale[0]) + ' ' + ' 12 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[0])+ '_length12.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[1]) + ' ' + ' 12 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[1])+ '_length12.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[2]) + ' ' + ' 12 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[2])+ '_length12.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[3]) + ' ' + ' 12 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[3])+ '_length12.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[4]) + ' ' + ' 12 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[4])+ '_length12.log' + '\n'
			cmd += IEDB + ' ' + method + ' ' + str (HLAale[5]) + ' ' + ' 12 ' + pep_fasta
			cmd += ' > ' + out_dir + '/data/neoantigen/classI/classI_' + sample_id + '_' + method + '_'+ str (HLAale[5])+ '_length12.log' + '\n'

		sh.write (cmd)
	sh.close ()
	pts.close ()

	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)

if __name__ == '__main__':
	main()

'''

python /data/cyu/scripts/wes/neoantigen/mhcI_v1.py \
-i /data/cyu/gbm/data/neoantigen/ \
-o /data/cyu/gbm/ \
-s 2 -e 2 \
-p /data/cyu/gbm/sample_pairs.txt

python /data/cyu/scripts/wes/neoantigen/mhcI_v1.py \
-i /data/cyu/gbm/data/neoantigen/ \
-o /data/cyu/gbm/ \
-s 5 -e 5 \
-p /data/cyu/gbm/sample_pairs.txt

python /data/cyu/scripts/wes/neoantigen/mhcI_v1.py \
-i /data/cyu/gbm/data/neoantigen/ \
-o /data/cyu/gbm/ \
-s 3 -e 3 \
-p /data/cyu/gbm/sample_pairs.txt

python /data/cyu/scripts/wes/neoantigen/mhcI_v1.py \
-i /data/cyu/gbm/data/neoantigen/ \
-o /data/cyu/gbm/ \
-s 1 -e 1 \
-p /data/cyu/gbm/sample_pairs.txt


'''

