import sys
import matplotlib
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Mutation calling')
	parser.add_argument ('-i', dest='input_dir', help='The input directionary path')
	parser.add_argument ('-p', dest='pt_file', help='Input patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')
	parser.add_argument ('-o', dest='outdir', help='The output path for filtered annovar annotation')
	args = parser.parse_args ()
	return args

def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	outdir = args.input_dir
	pt_file = args.pt_file
	start = args.start
	end = args.end

	sh_dir = outdir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = outdir + "/log/"  # The script log/alignment_start_end.sh will be executed
	if os.path.isdir(outdir + "/sh/") and os.path.isdir(outdir + "/log/"):
		pass
	else:
		subprocess.call ('mkdir -p ' + outdir + "/{sh,log}/", shell=True)
	if os.path.isdir (outdir + '/data/neoantigen/'):
		pass
	else:
		subprocess.call ('mkdir -p ' + outdir + '/data/neoantigen/', shell=True)
	if os.path.isdir(outdir + '/data/neoantigen/IEDB_ii/'):
		pass
	else:
		subprocess.call ('mkdir -p ' + outdir + '/data/neoantigen/IEDB_ii/', shell=True)
	fout = sh_dir + 'IEDB_ii_summary' + '_' + start + '_' + end + '.sh'
	sh = open (fout, 'w')
	pts = open (pt_file)  # note there is a header
	log = log_dir + 'IEDB_ii_summary' + '_' + start + '_' + end + '.log'
	lns = pts.readlines ()
	iedb_inputdir = outdir + '/data/neoantigen/inputHLA/'
	pep_inputdir = input_dir + '/neoantigen/pep_fasta/'
	outdir_iedbii = outdir + '/data/neoantigen/IEDB_ii/'
	outdir_classII_summary = outdir + '/data/neoantigen/summary/IEDB_ii/'

	if os.path.isdir(iedb_inputdir):
		pass
	else:
		subprocess.call ('mkdir -p ' + iedb_inputdir, shell=True)
	if os.path.isdir(outdir_iedbii + '/sh/') and os.path.isdir (outdir_iedbii + '/log/'):
		pass
	else:
		subprocess.call ('mkdir -p ' + iedb_inputdir + '/{sh,log}/', shell=True)

	if os.path.isdir(outdir_classII_summary ):
		pass
	else:
		subprocess.call ('mkdir -p ' + outdir_classII_summary, shell=True)

	cmd = ''

	for i in range ((int) (start), (int) (end) + 1):
		print (lns[i])
		tmp = lns[i].strip ('\n').split (',')
		ptID = tmp[1]  # tumor patient ID
		pnID = tmp[0]  # normal patient ID
		pepfasta = pep_inputdir + ptID + '.classII.fasta'
		pep_loc = input_dir + '/data/neoantigen/pep_fasta/' + ptID + '.pep.loc'


		cmd += 'echo summary of IEDB classII \n\n'
		cmd += 'python /data/cyu/scripts/wes/4_neoantigen/src/mhcII_summary.py'
		cmd += ' -i ' + input_dir
		cmd += ' -o ' + outdir_classII_summary
		cmd += ' -p ' + pt_file
		cmd += ' -s ' + start + ' -e ' + end
		cmd += ' & wait; \n'

		cmd += 'sh /data/cyu/scripts/wes/4_neoantigen/src/merge_sample_II.sh '
		cmd += ptID + ' ' + outdir_classII_summary + ' & wait; \n\n'

		if os.path.isdir(outdir_classII_summary + '/' + ptID + '_length25_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length25_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass

		if os.path.isdir (outdir_classII_summary + '/' + ptID + '_length24_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length24_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass

		if os.path.isdir (outdir_classII_summary + '/' + ptID + '_length23_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length23_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass
		if os.path.isdir (outdir_classII_summary + '/' + ptID + '_length22_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length22_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass
		if os.path.isdir (outdir_classII_summary + '/' + ptID + '_length21_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length21_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass

		if os.path.isdir (outdir_classII_summary + '/' + ptID + '_length20_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length20_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass

		if os.path.isdir (outdir_classII_summary + '/' + ptID + '_length19_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length19_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass

		if os.path.isdir (outdir_classII_summary + '/' + ptID + '_length18_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length18_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass

		if os.path.isdir (outdir_classII_summary + '/' + ptID + '_length17_iedb_mhcII.txt'):
			cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
			cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length17_iedb_mhcII.txt'
			cmd += ' -r ' + pep_loc
			cmd += ' -o ' + outdir_classII_summary
			cmd += ' & wait; \n'
		else:
			pass
		cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
		cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length16_iedb_mhcII.txt'
		cmd += ' -r ' + pep_loc
		cmd += ' -o ' + outdir_classII_summary
		cmd += ' & wait; \n'

		cmd += '/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/src/match_classII.py '
		cmd += ' -i ' + outdir_classII_summary + '/' + ptID + '_length15_iedb_mhcII.txt'
		cmd += ' -r ' + pep_loc
		cmd += ' -o ' + outdir_classII_summary
		cmd += ' & wait; \n'
		sh.write (cmd)
	sh.close ()
	pts.close ()

	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	#	print ('Submitting mutation calling jobs to process ' + ptID + '\n\n')
	subprocess.call (cmd, shell=True)


if __name__ == '__main__':
	main ()

'''
python /data/cyu/scripts/wes/4_neoantigen/run_iedbII_summary.py \
-i /data/gbm/dna/second/ \
-p /data/gbm/dna/second/pairs.txt \
-s 3 -e 3
'''
