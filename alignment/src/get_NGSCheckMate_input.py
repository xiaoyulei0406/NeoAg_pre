import os, sys, subprocess
import argparse
import glob


def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
	description='Run input file for NGScheckMate')
	parser.add_argument ('-i', dest='input_dir', help='Input fastq directory')
	parser.add_argument ('-o', dest='output_dir', help='Output fastq list for NGScheckMate input format')
	args = parser.parse_args ()
	return args

def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	output_dir = args.output_dir
	outfile = output_dir + '/fastqlist.txt'
	fw=open(outfile,'w')

	allfiles_arr=[]
	pairs=0
	#get full path of all the fastq files
	for file in glob.glob(input_dir+'/*.gz'):
		allfiles_arr.append(file)
	#allfiles_arr = allfiles_arr.sort()
	#allfiles = '\n'.join(allfiles_arr)
	#print(allfiles)
	print(allfiles_arr)
	for line in allfiles_arr :
		print (line)
		filename=str(line).split("/")[-1]
		p_id = filename.split ('_R')[0]
		if pairs==0 :
			fw.writelines(str(line)+'\t')
			pairs=+1
		elif pairs==1 :
			sample_id = str (line).split("_")[0]
			fw.writelines (str (line) + '\t' + p_id + '\n')
			pairs=0
	fw.close()

if __name__ == '__main__':
	main ()
