import argparse
import os
import pandas as pd
import glob
import openpyxl
def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
		description='Check if there are any files mislabelled')
	parser.add_argument ('input', help='Input file of NGSCheckMate output_matched.txt', type=str)
	parser.add_argument ('output', help='Output for mislabelled files', type=str)
	args = parser.parse_args ()
	return args


def main():
	args = parse_arguments ()
	input = args.input
	output = args.output
	misArr=[]

	with open (input, 'r') as f:
		for line in f.readlines():
			if "matched" in line:
				l1_tmp = line.split ("\t")[0]
				l1 = str(l1_tmp).split ("_")[0]
				print (l1)
				l2_tmp = line.split ("\t")[2]
				l2 = str (l2_tmp).split ("_")[0]
				#print("l1: "+l1)
				#print("l2: " +l2)
				if (l1!= l2):
					misLabel = "{0},{1},{2},{3},{4}".format (line.split ("\t")[0], line.split ("\t")[1], line.split ("\t")[2], line.split ("\t")[3],
																 line.split ("\t")[4])
					misLabels = list (misLabel.split (","))
					misArr.append (misLabels)
			else:
				continue

	y = pd.DataFrame (misArr)
	y.to_excel (output, index=False, float_format="%.2f")

if __name__ == '__main__':
	main ()



