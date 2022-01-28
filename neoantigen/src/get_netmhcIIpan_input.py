import sys
import matplotlib
import pandas as pd
import numpy as np
import argparse


def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Get input for netmhcIIpan')
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
		print("Running sample ID :")
		print (ptID)
		hla_file = input_dir + '/' + ptID + '/' + ptID +'_report.d4.txt'
		hlas = open (hla_file)
		hla_list1 = hlas.readlines ()
		outfile=out_dir + '/' + ptID + '.hla.mhcpanII.txt'
		fw=open(outfile, 'w')
		hlaalleles = []
		flag_drb1 = "Not exist"
		flag_drb3 = "Not exist"
		flag_drb4 = "Not exist"
		flag_drb5 = "Not exist"
		flag=len(hla_list1)
		i=4
		DQA1aB1a = ''
		DQA1aB1b = ''
		DQA1bB1a = ''
		DQA1bB1b = ''
		DRB1a = ''
		DRB1b = ''
		DRB3a = ''
		DRB3b = ''
		DRB4a = ''
		DRB4b = ''
		DRB5a = ''
		DRB5b = ''
		for line in hla_list1:
			if line.startswith ('A') or line.startswith ('B') or line.startswith ('C'):
				continue
			if i>flag:
				break
				# may lose DRB1,DRB3,DRB4,DRB5
			elif line.startswith('DPA1'):
				line = line.strip ("\n").split ("\t")
				line[1]=line[1].replace("*",'')
				line[1] = line[1].replace (":", '')
				line[2] = line[2].replace ("*", '')
				line[2] = line[2].replace (":", '')
				DPA1a = "HLA-" + line[1]
				DPA1b = "HLA-" + line[2]
			elif line.startswith ('DPB1'):
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].replace ("*", '')
				line[1] = line[1].replace (":", '')
				line[2] = line[2].replace ("*", '')
				line[2] = line[2].replace (":", '')
				DPB1a = line[1]
				DPB1b = line[2]

				DPA1aB1a = DPA1a + '-' + DPB1a +','
				DPA1aB1b = DPA1a + '-' + DPB1b +','
				DPA1bB1a = DPA1b + '-' + DPB1a +','
				DPA1bB1b = DPA1b + '-' + DPB1b +','
			elif line.startswith ('DQA1'):
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].replace ("*", '')
				line[1] = line[1].replace (":", '')
				line[2] = line[2].replace ("*", '')
				line[2] = line[2].replace (":", '')
				DQA1a = "HLA-" + line[1]
				DQA1b = "HLA-" + line[2]
			elif line.startswith ('DQB1'):
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].replace ("*", '')
				line[1] = line[1].replace (":", '')
				line[2] = line[2].replace ("*", '')
				line[2] = line[2].replace (":", '')
				DQB1a = line[1]
				DQB1b = line[2]
				DQA1aB1a = DQA1a + '-' + DQB1a +','
				DQA1aB1b = DQA1a + '-' + DQB1b +','
				DQA1bB1a = DQA1b + '-' + DQB1a +','
				DQA1bB1b = DQA1b + '-' + DQB1b +','
			elif line.startswith ('DRB1') and flag_drb1 =="Not exist":
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].replace ("*", '_')
				line[1] = line[1].replace (":", '')
				line[2] = line[2].replace ("*", '_')
				line[2] = line[2].replace (":", '')
				DRB1a = line[1]
				DRB1b = line[2]
				flag_drb1 = "exist"
			elif line.startswith ('DRB3') and flag_drb1 =="Not exist":
				DRB1a = ''
				DRB1b = ''
			elif line.startswith ('DRB3') and flag_drb1 == "exist":
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].replace ("*", '_')
				line[1] = line[1].replace (":", '')
				line[2] = line[2].replace ("*", '_')
				line[2] = line[2].replace (":", '')
				if line[1] != '':
					DRB3a = line[1] +','
				elif line[1] == '':
					DRB3a = ''
				if line[2] != '':
					DRB3b = line[2] +','
				elif line[2] == '':
					DRB3b = ''
				flag_drb3 = "exist"
			elif line.startswith ('DRB4') and flag_drb3 =="Not exist":
				DRB3a = ''
				DRB3b = ''
			elif line.startswith ('DRB4') and flag_drb3 == "exist":
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].replace ("*", '_')
				line[1] = line[1].replace (":", '')
				line[2] = line[2].replace ("*", '_')
				line[2] = line[2].replace (":", '')
				if line[1] != '':
					DRB4a = line[1] +','
				elif line[1] == '':
					DRB4a = ''
				if line[2] != '':
					DRB4b = line[2] +','
				elif line[2] == '':
					DRB4b = ''
				flag_drb4 = "exist"
			elif line.startswith ('DRB5') and flag_drb4 =="Not exist":
				DRB4a = ''
				DRB4b = ''
			elif line.startswith ('DRB5') and flag_drb4 == "exist":
				line = line.strip ("\n").split ("\t")
				line[1] = line[1].replace ("*", '_')
				line[1] = line[1].replace (":", '')
				line[2] = line[2].replace ("*", '_')
				line[2] = line[2].replace (":", '')
				if line[1]!='':
					DRB5a = line[1] +','
				elif line[1] == '':
					DRB5a = ''
				if line[2]!='':
					DRB5b = line[2] +','
				elif line[2] == '':
					DRB5b = ''
				flag_drb5 = "exist"
			elif 'DRB5' not in line and flag_drb5 == "Not exist":
				DRB5a = ''
				DRB5b = ''
			i=i+1
		hlaalleles.append (DPA1aB1a)
		hlaalleles.append (DPA1aB1b)
		hlaalleles.append (DPA1bB1a)
		hlaalleles.append (DPA1bB1b)
		hlaalleles.append (DQA1aB1a)
		hlaalleles.append (DQA1aB1b)
		hlaalleles.append (DQA1bB1a)
		hlaalleles.append (DQA1bB1b)
		hlaalleles.append (DRB1a)
		hlaalleles.append (DRB1b)
		hlaalleles.append (DRB3a)
		hlaalleles.append (DRB3b)
		hlaalleles.append (DRB4a)
		hlaalleles.append (DRB4b)
		hlaalleles.append (DRB5a)
		hlaalleles.append (DRB5b)
		hlaalleles = list (set (hlaalleles))  # Remove duplicate alleles if there are any
		hlaalleles = sorted(hlaalleles) # Sort
		hlastring = '\n'.join (hlaalleles)
		#hlastring=hlastring.replace(',,',',')
		fw.writelines(hlastring)

	pts.close ()


if __name__ == '__main__':
	main ()


