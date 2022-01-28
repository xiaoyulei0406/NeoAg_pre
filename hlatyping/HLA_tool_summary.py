#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import glob

def get_opti(input_opti,output):
	list_hla = []
	fmt_converter = lambda x: x.lower ().replace ('-', '_').replace ('*', '_').replace (':', '_')
	fw=open(output,'w')
	fw.write ("sample id" + '\t' + 'HLA-A1' + '\t' + 'HLA-A2' + '\t' + 'HLA-B1' + '\t' + 'HLA-B2' + '\t' + 'HLA-C1' + '\t' + 'HLA-C2' + '\n')
	for file in glob.glob (os.path.join (input_opti, '*.tsv')):  # generator, search immediate subdirectories
		list_hla=[]
		sample_id_tmp1 = str(file).split ("/HLAtyping/optitype/")[0] ##note for sampleID location
		sample_id_tmp = sample_id_tmp1.split ("/2021")[0]  ##note for sampleID location
		sample_id = sample_id_tmp.split ("/")[-1]  ##note for sampleID location
		line=sample_id+"\t"
		handle=open (file,'r')
		lines = handle.readlines ()
		list_hla.extend (['hla_' + fmt_converter (x) for x in lines[1].split ('\t')[1:7]])
		fw.write (line+ '\t'.join(list_hla) + '\n')
	fw.close()

def get_HLAVBSeq(input_vb,outputI,outputII):
	list_hla = []
	fmt_converter = lambda x: x.lower ().replace ('-', '_').replace ('*', '_').replace (':', '_')
	fw = open (outputI, 'w')
	fw.write ("sample id" + '\t' + 'HLA-A1' + '\t' + 'HLA-A2' + '\t' + 'HLA-B1' + '\t' + 'HLA-B2' + '\t' + 'HLA-C1' + '\t' + 'HLA-C2' + '\n')

	for file in glob.glob (os.path.join (input_vb, '*_report.d4.txt')):  # generator, search immediate subdirectories
		list_hla = []
		sample_id_tmp = file.split ("/")[-1]  ##note for sampleID location
		sample_id = sample_id_tmp.split ("_report.d4.txt")[0]  ##note for sampleID location
		line = sample_id + "\t"
		handle = open (file, 'r')
		lines = handle.readlines ()
		for i in range(1,4):
			lines[i] = lines[i].replace ('\n', '')
			list_hla.extend (['hla_' + fmt_converter (x) for x in lines[i].split ('\t')[1:7]])
		fw.write (line + '\t'.join (list_hla) + '\n')
	fw.close ()


	fmt_converter1 = lambda x: x.upper ().replace ('-', '_').replace ('*', '_').replace (':', '_')
	fw1 = open (outputII, 'w')
	for file in glob.glob (os.path.join (input_vb, '*_report.d4.txt')):  # generator, search immediate subdirectories

		sample_id_tmp = file.split ("/")[-1]  ##note for sampleID location
		sample_id = sample_id_tmp.split ("_report.d4.txt")[0]  ##note for sampleID location
		sid = sample_id + "\t"
		handle = open (file, 'r')
		lines = handle.readlines ()
		fw1.write (sid)
		list_hla = []
		flag_drb3=0
		flag_drb4=0
		for line in lines:
			if line.startswith('A') or line.startswith ('B') or line.startswith ('C'):
				continue
			elif line.startswith('DPA1'):
				#line = line.replace ('\n', '')
				#list_hla.extend ([fmt_converter1 (x) for x in line.split ('\t')[1:7]])
				line=line.strip("\n").split("\t")
				list_hla.extend ("HLA-"+line[1]+"\t"+ "HLA-"+line[2]+"\t")
			elif line.startswith('DPB1'):
				#line = line.replace ('\n', '')
				#list_hla.extend ([fmt_converter1 (x) for x in line.split ('\t')[1:7]])
				line = line.strip ("\n").split ("\t")
				list_hla.extend ("HLA-" + line[1] + "\t" + "HLA-" + line[2]+"\t")
			elif line.startswith('DQA1'):
				line = line.replace ('\n', '')
				#list_hla.extend ([fmt_converter1 (x) for x in line.split ('\t')[1:7]])
				line = line.strip ("\n").split ("\t")
				list_hla.extend ("HLA-" + line[1] + "\t" + "HLA-" + line[2]+"\t")
			elif line.startswith('DQB1'):
				line = line.replace ('\n', '')
				#list_hla.extend ([fmt_converter1 (x) for x in line.split ('\t')[1:7]])
				line = line.strip ("\n").split ("\t")
				list_hla.extend ("HLA-" + line[1] + "\t" + "HLA-" + line[2]+"\t")
			elif line.startswith('DRB1'):
				line = line.replace ('\n', '')
				#list_hla.extend ([fmt_converter1 (x) for x in line.split ('\t')[1:7]])
				line = line.strip ("\n").split ("\t")
				list_hla.extend ("HLA-" + line[1] + "\t" + "HLA-" + line[2]+"\t")
			elif line.startswith('DRB3'):
				line = line.replace ('\n', '')
				#list_hla.extend ([fmt_converter1 (x) for x in line.split ('\t')[1:7]])
				line = line.strip ("\n").split ("\t")
				list_hla.extend ("HLA-" + line[1] + "\t" + "HLA-" + line[2]+"\t")
				flag_drb3=1
			elif line.startswith('DRB4'):
				line = line.replace ('\n', '')
				#list_hla.extend ([fmt_converter1 (x) for x in line.split ('\t')[1:7]])
				line = line.strip ("\n").split ("\t")
				list_hla.extend ("HLA-" + line[1] + "\t" + "HLA-" + line[2]+"\t")
				flag_drb4 = 1
			elif 'DRB3' not in line:
				if flag_drb3 == 1 and flag_drb4 == 1:
					pass
				elif flag_drb3 == 1 and flag_drb4 == 0:
					#list_hla.extend (["DRB4 not found" + '\t' + "DRB4 not found"])
					list_hla.extend (["_" + '\t' + "_"+"\t"])
				elif flag_drb3 == 0 and flag_drb4 == 1:
					#list_hla.extend (["DRB3 not found" + '\t' + "DRB3 not found"])
					list_hla.extend (["_" + '\t' + "_"+"\t"])
				else:
					pass

		#fw1.writelines ('\t'.join (list_hla) + '\n')
		fw1.writelines (''.join (list_hla) + '\n')

	fw1.close ()


def get_polysolver(input_poly,output):
	fmt_converter = lambda x: x.lower ().replace ('-', '_').replace ('*', '_').replace (':', '_')
	fw=open(output,'w')
	fw.write ("sample id" + '\t' + 'HLA-A1' + '\t' + 'HLA-A2' + '\t' + 'HLA-B1' + '\t' + 'HLA-B2' + '\t' + 'HLA-C1' + '\t' + 'HLA-C2' + '\n')
	for file in glob.glob (os.path.join (input_poly, 'winners.hla.txt')):  # generator, search immediate subdirectories
		list_hla=[]
		sample_id = file.split ("/")[-2] ##note for sampleID location
		sid=sample_id+"\t"
		fw.write(sid)
		handle=open (file,'r')
		lines = handle.readlines ()
		for line in lines:
			line=line.replace("\n",'')
			list_hla.extend ([fmt_converter (x) for x in line.split ('\t')[1:7]])
		fw.write ('\t'.join(list_hla) + '\n')
	fw.close()



def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
	description='Obtain consensus hla alleles based on HLA-VBSeq, OptiType, Polysolver')
	parser.add_argument ('-i', dest='input_dir', help='Input directory')

	args = parser.parse_args ()
	return args
def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	out_dir = args.input_dir

	input_opti = input_dir + '/optitype/*/2*/'
	input_vb = input_dir + '/HLAVBSeq/*/'
	input_poly = input_dir + '/PolySolver/*/'
	output_opti = out_dir + '/' + 'optitype.summary.classI.txt'
	output_vbI=out_dir + '/' + 'HLAVBSeq.summary.classI.txt'
	output_vbII = out_dir + '/' + 'HLAVBSeq.summary.classII.txt'
	output_poly = out_dir + '/' + 'PolySolver.summary.classI.txt'
	list_hla_vb = get_HLAVBSeq(input_vb,output_vbI,output_vbII)
	list_hla_op = get_opti (input_opti,output_opti)
	list_hla_poly=get_polysolver(input_poly,output_poly)




if __name__ == '__main__':
	main ()

'''
python /data/cyu/scripts/wes/3_hlatyping/HLA_tool_summary.py \
-i /data/cyu/gbm/data/HLAtyping/ \
-o /data/cyu/gbm/data/HLAtyping/

python /data/cyu/scripts/wes/3_hlatyping/HLA_tool_summary.py \
-i /data/cyu/gbm/data/HLAtyping/

'''
