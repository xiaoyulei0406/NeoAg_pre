import sys
import matplotlib
import pandas as pd
import numpy as np
import argparse


def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Get input for IEDB classI')
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

	pts = open (pt_file)  # note there is a header
	lns = pts.readlines ()

	for i in range ((int) (start)-1, (int) (end)):
		tmp = lns[i].split ('\n')
		ptID = tmp[0]  # tumor patient ID
		hla_file = input_dir + '/' + ptID +'.hla.consensus.txt'
		hlas = open (hla_file)
		hla_list1 = hlas.readlines ()
		outfile=out_dir + '/' + ptID + '.hla.IEDBI.txt'
		fw=open(outfile, 'w')
		hlaalleles = []
		A1a = ''
		A1b = ''
		B1a = ''
		B1b = ''
		C1a = ''
		C1b = ''

		for line in hla_list1:
			if line.startswith ('HLA-A'):
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].upper()
				line[1] = line[1].replace ("HLA_A_", 'HLA-A*')
				line[1] = line[1].replace ("_", ':')
				line[2] = line[2].replace ("hla_a_", 'HLA-A*')
				line[2] = line[2].replace ("_", ':')
				A1a = line[1]
				A1b = line[2]
			elif line.startswith ('HLA-B'):
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].upper ()
				line[1] = line[1].replace ("HLA_B_", 'HLA-B*')
				line[1] = line[1].replace ("_", ':')
				line[2] = line[2].replace ("hla_b_", 'HLA-B*')
				line[2] = line[2].replace ("_", ':')
				B1a = line[1]
				B1b = line[2]
			elif line.startswith ('HLA-C'):
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].upper ()
				line[1] = line[1].replace ("HLA_C_", 'HLA-C*')
				line[1] = line[1].replace ("_", ':')
				line[2] = line[2].replace ("hla_c_", 'HLA-C*')
				line[2] = line[2].replace ("_", ':')
				C1a = line[1]
				C1b = line[2]
			else :
				continue
		hlaalleles.append (A1a)
		hlaalleles.append (A1b)
		hlaalleles.append (B1a)
		hlaalleles.append (B1b)
		hlaalleles.append (C1a)
		hlaalleles.append (C1b)

		hlaalleles = list (set (hlaalleles))  # Remove duplicate alleles if there are any
		hlaalleles = sorted(hlaalleles) # Sort
		hlastring = '\n'.join (hlaalleles)
		#hlastring=hlastring.replace(',,',',')
		fw.writelines(hlastring)

	pts.close ()


if __name__ == '__main__':
	main ()

'''

/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/4_neoantigen/get_IEDB_input.py \
-i /data/cyu/gbm/data/HLAtyping/summary/ \
-p /data/cyu/gbm/samples.txt \
-s 1 -e 10 -o /data/cyu/gbm/data/neoantigen/inputHLA/

'''

