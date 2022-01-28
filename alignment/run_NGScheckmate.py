import argparse, os, sys, subprocess


def parse_arguments():
	description=("Running NGScheckmate for sample labelling."
	)
	parser = argparse.ArgumentParser (description=description)
	parser.add_argument ('-i', dest='input_dir', help='Input fastq directory')
	#parser.add_argument ('-p', dest='pt_file', help='Input patient file')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	args = parser.parse_args ()
	return args


def main():
	args = parse_arguments ()
	out_dir = args.out_dir
	input_dir = args.input_dir


	sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed
	subprocess.call ('mkdir -p ' + './' + out_dir + "/{sh,log}/", shell=True)
	subprocess.call ('mkdir -p ' + './' + out_dir + '/data/NGScheckmate/', shell=True)

	fout = sh_dir + 'NGScheckmate.sh'
	sh = open (fout, 'w')
	log = log_dir + 'NGScheckmate.log'
	out_dir1 = out_dir +'/'

	cmd = 'echo \'###########################################\' >>' + log + ' 2>&1 & \n'
	cmd += 'echo Running NGScheckMate for sample labelling >>' + log + ' 2>&1 & \n'
	cmd += '\n'
	cmd += 'echo Below are the tools information >>' + log + ' 2>&1 & \n'
	cmd += 'echo NGSCheckMate version: v1>>' + log + ' 2>&1 & \n'
	cmd += '\n'
	cmd += 'echo \'###########################################\' >>' + log + ' 2>&1 & \n'
	cmd += '\n\n\n'

	cmd += 'echo Get NGSCheckMate input format at `date` \n'
	cmd += 'echo \'###########################################\' \n'

	cmd += 'python3 /data/cyu/scripts/wes/1_alignment/src/get_NGSCheckMate_input.py'
	cmd += ' -i ' + input_dir
	cmd += ' -o ' + out_dir1 + ' & wait; \n'

	cmd += 'echo starting NGSCheckMate at `date` \n\n'
	cmd += 'echo 1.Running NGSCheckMate ' +  '...\n\n'

	pt_file = out_dir1 + '/fastqlist.txt'

	cmd += 'python /home/cyu/tools/NGSCheckMate/ncm_fastq_cyu.py -p 4 '
	cmd += '-l ' + pt_file + ' -pt /home/cyu/tools/NGSCheckMate/SNP/SNP.pt '
	cmd += '-O ' + out_dir1 + '/data/NGScheckmate/' + ' -nz & wait; \n'

	cmd += 'echo NGSCheckMate is Done \n\n'

	cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/1_alignment/src/run_checkpairs.py '
	cmd += out_dir1 + '/data/NGScheckmate/' + "output_matched.txt" + " " + out_dir1 + 'mislabel.xls'
	cmd += ' & wait; \n'

	sh.write (cmd)
	sh.close ()
	cmd = 'nohup bash ' + fout + '>>' + log + ' 2>&1 &'
	print ('Submitting NGSCheckMate jobs to process ' + '\n\n')
	subprocess.call (cmd, shell=True)

if __name__ == '__main__':
	main()

