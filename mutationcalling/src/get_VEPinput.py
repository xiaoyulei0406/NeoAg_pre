import glob
import os
import argparse
import pandas as pd

def get_VEP_SNP(input_dir,out_dir):
	for filepath in glob.glob (os.path.join (input_dir, '*.SNP.combine.txt')):
		filename_tmp = filepath.split ("/")[-1]
		filename = filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1] + ".vep.txt"
		with open ('{}/{}'.format (out_dir, filename), 'w') as fw:
			with open (filepath, 'r') as fp:
				# header skip
				fw.writelines(fp.readline ())
				for readline in fp:
					ar = readline.strip ().split ('\t')
					if (str (ar[9]) != "." and str (ar[7]) != "UNKNOWN"):
						refaa = str (ar[9]).split (":")[3]
						refP = str (refaa).split (".")[1]
						if str (ar[3]) == str (refP[0]):
							fw.writelines (str (ar[0]) + "\t" + str (ar[1]) + "\t" + str (ar[2]) + "\t" + str (ar[3]) + "/" + str (ar[4]) + "\t" + "+" + "\n")
						if str (ar[3]) != str (refP[0]):
							if len (str (ar[3])) == 1:
								fw.writelines (str (ar[0]) + "\t" + str (ar[1]) + "\t" + str (ar[2]) + "\t" + str (refP[0]) + "/" + str (refP[-1]) + "\t" + "-" + "\n")
							if len (str (ar[3])) == 2:
								refPP = str (refaa).split (".")[-2:-1]
								if (ar[4]) == str (refPP):
									fw.writelines (str (ar[0]) + "\t" + str (ar[1]) + "\t" + str (ar[2]) + "\t" + str (ar[3]) + "/" + str (ar[4]) + "\t" + "+" + "\n")
								if (ar[4]) != str (refPP):
									fw.writelines (str (ar[0]) + "\t" + str (ar[1]) + "\t" + str (ar[2]) + "\t" + str (ar[3]) + "/" + str (ar[4]) + "\t" + "-" + "\n")

		infile = out_dir + "/" + filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1] + ".vep.txt"
		df = pd.read_csv (infile, sep='\t')
		df = df.sort_values (by=['Chromosome', 'Start_Position'])
		df_1 = df.iloc[:, [0, 1, 2, 3, 4]]
		df_1.to_csv (out_dir + "/" + filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1]  +".sorted.vep.txt", index=False, header=False, sep='\t')


def get_VEP_INDEL(input_dir,out_dir):
	for filepath in glob.glob (os.path.join (input_dir, '*.INDEL.combine.txt')):
		filename_tmp = filepath.split ("/")[-1]
		filename = filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1] + ".vep.txt"
		with open ('{}/{}'.format (out_dir, filename), 'w') as fw:
			with open (filepath, 'r') as fp:
				# header skip
				fw.writelines(fp.readline ())
				for readline in fp:
					ar = readline.strip ().split ('\t')
					if ("wholegene" in str (ar[9])):
						fw.writelines (str (ar[0]) + "\t" + str (ar[1]) + "\t" + str (int(ar[2])) + "\t" + str (ar[3]) + "/" + str (ar[4]) + "\t" + "+" + "\n")
					else:
						refaa = str (ar[9]).split (":")[3]
						refP = str (refaa).split (".")[1]
						print (refP)
						if (str (ar[3]) != "-"):
							fw.writelines (str (ar[0]) + "\t" + str (int(ar[1])) + "\t" + str (int(ar[2])) + "\t" + str (ar[3]) + "/" + str (ar[4]) + "\t" + "+" + "\n")
						if (str (ar[4]) != "-"):
							fw.writelines (str (ar[0]) + "\t" + str (int(ar[1])+1) + "\t" + str (int(ar[2])) + "\t" + "-/" + str (ar[4]) + "\t" + "+" + "\n")
	#			appendFile.close()
		infile = out_dir + "/" + filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1] + ".vep.txt"
		print(infile)
		df = pd.read_csv (infile, sep='\t')
		df = df.sort_values (by=['Chromosome', 'Start_Position'])
		df_1=df.iloc[:,[0,1,2,3,4]]
		df_1.to_csv (out_dir + "/" + filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1]  +".sorted.vep.txt", index=False, header=False, sep='\t')


def get_VEP_MNP(input_dir,out_dir):
	for filepath in glob.glob (os.path.join (input_dir, '*.MNP.combine.txt')):
		print (filepath)
		filename_tmp = filepath.split ("/")[-1]
		filename = filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1] + ".vep.txt"
		with open ('{}/{}'.format (out_dir, filename), 'w') as fw:
			with open (filepath, 'r') as fp:
				# header skip
				fw.writelines(fp.readline ())
				for readline in fp:
					#				print(readline)
					ar = readline.strip ().split ('\t')
					if ("wholegene" in str (ar[9])):
						fw.writelines (str (ar[0]) + "\t" + str (ar[1]) + "\t" + str (ar[2]) + "\t" + str (ar[3]) + "/" + str (
								ar[4]) + "\t" + "+" + "\n")
					else:
						refaa = str (ar[9]).split (":")[3]
						refP = str (refaa).split (".")[1]
						print (refP)
						if (str (ar[4]) not in str (ar[9]).split (":")[3]):
							fw.writelines (
								str (ar[0]) + "\t" + str (ar[1]) + "\t" + str (ar[2]) + "\t" + str (ar[3]) + "/" + str (
									ar[4]) + "\t" + "-" + "\n")
						if (str (ar[4]) in str (ar[9]).split (":")[3]):
							fw.writelines (
								str (ar[0]) + "\t" + str (ar[1]) + "\t" + str (ar[2]) + "\t" + str (ar[3]) + "/" + str (
									ar[4]) + "\t" + "+" + "\n")

		infile = out_dir + "/" + filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1] + ".vep.txt"
		df = pd.read_csv (infile, sep='\t')
		df = df.sort_values (by=['Chromosome', 'Start_Position'])
		df_1=df.iloc[:,[0,1,2,3,4]]
		df_1.to_csv (out_dir + "/" + filename_tmp.split (".")[0] + "." + filename_tmp.split (".")[1]  +".sorted.vep.txt", index=False, header=False, sep='\t')


def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Combine all the annotation results')
	parser.add_argument ('-i', dest='input_dir', help='The inputdir path')
	parser.add_argument ('-o', dest='outdir', help='The output path')
	args = parser.parse_args ()
	return args


def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	out_dir = args.outdir

	get_VEP_SNP(input_dir,out_dir)
	get_VEP_INDEL(input_dir,out_dir)
	get_VEP_MNP(input_dir,out_dir)



if __name__ == '__main__':
	main()



"""
Get VEP input from annovar SNP and INDEL files
Usage:
python get_VEPinput.py <input_dir> <out_dir>
python /Users/cyu/Documents/work/mc3data/2020JUL14/somaticCalling/get_VEPinput.py -i /Users/cyu/Documents/work/mislabel_ICON/mafcombine/ \
-o /Users/cyu/Documents/work/mislabel_ICON/VEP/input/

/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/get_VEPinput.py -i /data/cyu/gbm/data/annovar/mafcombine/ \
-o /data/cyu/gbm/data/VEP/input


"""
