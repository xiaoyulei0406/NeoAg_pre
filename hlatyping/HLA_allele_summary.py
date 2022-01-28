#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import argparse
import glob

def csv_parser(input_opti): #optitype
	#print(input_opti)
	list_hla = []
	fmt_converter = lambda x: x.lower ().replace ('-', '_').replace ('*', '_').replace (':', '_')
	for file in glob.iglob (input_opti):  # generator, search immediate subdirectories
		with open (file) as handle:
			lines = handle.readlines ()
		if len (lines) == 2:
			list_hla.extend (['hla_' + fmt_converter (x) for x in lines[1].split ('\t')[1:7]])
	return list_hla

def txt_parser(input_vb): #HLAVBSeq
	list_hla = []
	fmt_converter = lambda x: x.lower ().replace ('-', '_').replace ('*', '_').replace (':', '_')
	with open (input_vb) as handle:
		lines = handle.readlines ()
		for i in range(1,4):
			lines[i]=lines[i].replace ('\n', '')
			list_hla.extend (['hla_' + fmt_converter (x) for x in lines[i].split ('\t')[1:7]])
	return list_hla

def write_output(list_hla, output):
	hla_genes = ['hla_a', 'hla_b', 'hla_c']
	with open (output, 'w') as handle:
		if not list_hla:  # empty list_hla
			print ("[WARNING] Nothing was written to {}!".format (output))
			return None
		for gene in hla_genes:
			try:  # if missing gene,then continue
				allele1, allele2 = [x for x in list_hla if x.startswith (gene)]
			except ValueError:
				continue
			if allele1 == allele2:  # only one uniq allele
				print ('[INFO] {} is HOM, allele={}'.format (gene, allele1))
				alleles = allele1 + '\t' + allele2
				line = gene.replace ('_', '-').upper () + '\t' + alleles + '\n'
				handle.write (line)
			else:
				alleles = allele1 + '\t' + allele2
				line = gene.replace ('_', '-').upper () + '\t' + alleles + '\n'
				handle.write (line)

def write_mhcpan_input(list_hla, output1):
	hla_genes = ['hla_a', 'hla_b', 'hla_c']
	with open (output1, 'w') as handle:
		if not list_hla:  # empty list_hla
			print ("[WARNING] Nothing was written to {}!".format (output1))
			return None
		for gene in hla_genes:
			try:  # if missing gene,then continue
				allele1, allele2 = [x for x in list_hla if x.startswith (gene)]
			except ValueError:
				continue
			if allele1 == allele2:  # only one uniq allele
				print ('[INFO] {} is HOM, allele={}'.format (gene, allele1))
				ale_tmp = allele1.split("_")
				ale1=ale_tmp[0].replace ('_', '-').upper ()+ale_tmp[2]+":"+ale_tmp[3]
				ale_tmp1 = allele2.split("_")
				ale2=ale_tmp1[0].replace ('_', '-').upper ()+ale_tmp1[2]+":"+ale_tmp1[3]
#				alleles = ale1 + ',' + ale2
				line = ale1 + ',' + ale2
				handle.write (line)
			else:
				ale_tmp = allele1.split ("_")
				ale1 = ale_tmp[0].replace ('_', '-').upper () + "-"+ ale_tmp[1].upper () + ale_tmp[2] + ":" + ale_tmp[3]
				ale_tmp1 = allele2.split ("_")
				ale2 = ale_tmp1[0].replace ('_', '-').upper () + "-"+ ale_tmp1[1].upper () + ale_tmp1[2] + ":" + ale_tmp1[3]
#				alleles = ale1 + ',' + ale2
				line = ale1 + ',' + ale2
				handle.write (line)
			if allele1.startswith("hla_c"):
				handle.write ("\n")
			else:
				handle.write (",")

def write_IEDB_input (list_hla, output2):
	hla_genes = ['hla_a', 'hla_b', 'hla_c']
	with open (output2, 'w') as handle:
		if not list_hla:  # empty list_hla
			print ("[WARNING] Nothing was written to {}!".format (output2))
			return None
		for gene in hla_genes:
			try:  # if missing gene,then continue
				allele1, allele2 = [x for x in list_hla if x.startswith (gene)]
			except ValueError:
				continue
			if allele1 == allele2:  # only one uniq allele
				print ('[INFO] {} is HOM, allele={}'.format (gene, allele1))
				ale_tmp = allele1.split ("_")
				ale1 = ale_tmp[0].replace ('_', '-').upper () + "*" + ale_tmp[2] + ":" + ale_tmp[3]
				ale_tmp1 = allele2.split ("_")
				ale2 = ale_tmp1[0].replace ('_', '-').upper () + "*" + ale_tmp1[2] + ":" + ale_tmp1[3]
				#				alleles = ale1 + ',' + ale2
				line = ale1 + ',' + ale2
				handle.write (line)
			else:
				ale_tmp = allele1.split ("_")
				ale1 = ale_tmp[0].replace ('_', '-').upper () + "-" + ale_tmp[1].upper () +  "*" +ale_tmp[2] + ":" + ale_tmp[3]
				ale_tmp1 = allele2.split ("_")
				ale2 = ale_tmp1[0].replace ('_', '-').upper () + "-" + ale_tmp1[1].upper () +  "*" +ale_tmp1[2] + ":" + ale_tmp1[3]
				#				alleles = ale1 + ',' + ale2
				line = ale1 + ',' + ale2
				handle.write (line)
			if allele1.startswith ("hla_c"):
				handle.write ("\n")
			else:
				handle.write (",")
def counts_parser(counts):
	zip_allele_scores = []
	with open (counts) as handle:
		for line in handle:
			allele, score = line.strip ().split ('\t')
			score = float (score)
			zip_allele_scores.append ((allele, score))
	return sorted (zip_allele_scores, key=lambda x: x[1], reverse=True)


def get_rank_score(zip_allele_scores_gene, allele):
	# zip_allele_scores_gene = [x for x in zip_allele_scores if x[0].startswith(gene)]
	rank = len (zip_allele_scores_gene) + 1
	score = None
	allele_poly = allele
	for i, item in enumerate (zip_allele_scores_gene):
		allele_poly, score_poly = item
		score_poly = float (score_poly)
		if i == 0:
			max_score = score_poly
		if allele_poly.startswith (allele):
			rank = i + 1
			score_percent = score_poly / max_score
			break
	return (allele_poly, rank, score_percent)

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
	description='Obtain consensus hla alleles based on HLA-VBSeq, OptiType, Polysolver')
	parser.add_argument ('-i', dest='input_dir', help='Input directory')
	parser.add_argument ('-p', dest='pt_file', help='Input patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')

	args = parser.parse_args ()
	return args
def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	out_dir = args.input_dir
	start = args.start
	end = args.end
	pt_file = args.pt_file
	pts = open (pt_file)  # there is a header
	lns = pts.readlines ()
	for i in range ((int) (start)-1, (int) (end)):
		tmp = lns[i].strip ('\n').split (',')
		sample_id = tmp[0]  # patient ID
		input_opti = input_dir + '/optitype/' + sample_id + '/2*/*.tsv'
		input_vb = input_dir + '/HLAVBSeq/' + sample_id + '/'+ sample_id + '_report.d4.txt'
		input_poly = input_dir + '/PolySolver/' + sample_id + '/'
		output = out_dir + '/' + sample_id + '.hla.consensus.txt'
		output1 = out_dir + '/' + sample_id + '.hla.mhcpanI.txt'
		output2 = out_dir + '/' + sample_id + '.hla.IEDBI.txt'
		list_hla_vb = txt_parser (input_vb)
		list_hla = csv_parser (input_opti)
		counts1 = os.path.join (input_poly, 'counts1.R0k6')
		counts2 = os.path.join (input_poly, 'counts2.R0k6')
		if not os.path.exists (counts1) or not os.path.exists (counts2):
			print ("[WARNING] Missing PolySolver counts file: {} and {}".format (counts1, counts2))
			write_output (list_hla, output)
		else:
			list_hla_out = []
			zip_allele1_scores = counts_parser (counts1)
			zip_allele2_scores = counts_parser (counts2)
			hla_genes = ['hla_a', 'hla_b', 'hla_c']
			for gene in hla_genes:
				zip_allele1_scores_gene = [x for x in zip_allele1_scores if x[0].startswith (gene)]
				zip_allele2_scores_gene = [x for x in zip_allele2_scores if x[0].startswith (gene)]
				allele1, allele2 = [x for x in list_hla if gene in x]
				#print(gene)
				#print(str(list_hla_vb))
				if gene not in str(list_hla_vb):
					allele1_vb, allele2_vb = [x for x in list_hla if gene in x]
				else:
					if 'hla_' in str(list_hla_vb):
						allele1_vb, allele2_vb = [x for x in list_hla if gene in x]
					else:
						allele1_vb, allele2_vb = [x for x in list_hla_vb if gene in x]
				#print(allele1_vb, allele2_vb)
#                print(allele1_vb, allele2_vb)
				#print(allele1, allele2)
				#print('''!!!''')
				poly_allele1 = zip_allele1_scores_gene[0][0]
				#print(poly_allele1)
				#print ('''!!!''')

				if (allele1 == allele2) and (allele1_vb==allele2_vb):
					list_hla_out.extend ([allele1, allele2])
					#print(list_hla_out)
					print ("[INFO] {} is HOM, allele={}, skipping this allele".format (gene, allele1))
					continue

				# if OptiType's allele1 and allele 2 dose not match PolySolver's allele1 either, then compare HLAVBSeq's allele1 and allele2
				if not poly_allele1.startswith (allele1) and not poly_allele1.startswith (allele2):
					print (
						"[INFO] OptiType two alleles({}, {}) do not matches PolySolver allele1({}). Try to get test HLAVBSeq allele in {}..." \
						.format (allele1, allele2, poly_allele1, counts1))

					if poly_allele1.startswith (allele1_vb):
						pass
					elif poly_allele1.startswith (allele2_vb):
						allele1_vb, allele2_vb = allele2_vb, allele1_vb  # swap allele1 and allele2:
					if not poly_allele1.startswith (allele1_vb) and not poly_allele1.startswith (allele2_vb):
						if allele1.startswith (allele1_vb):
							pass
						elif allele1.startswith (allele2_vb):
							allele1_vb, allele2_vb = allele2_vb, allele1_vb  # swap allele1 and allele2:
				if poly_allele1.startswith (allele1):
					pass
				elif poly_allele1.startswith (allele2):
					allele1, allele2 = allele2, allele1  # swap allele1 and allele2:

				poly_allele2, rank, score_percent = get_rank_score (zip_allele2_scores_gene, allele2)
				#print (poly_allele2)
				# High confidence allele2
				if rank < 100 and score_percent > 0.5:
					poly_allele2_tmp=zip_allele2_scores_gene[0][0]
					#print(poly_allele2_tmp)
					#print (allele2_vb)
					if not poly_allele2_tmp.startswith (allele1_vb) and not poly_allele2_tmp.startswith (allele2_vb):
						list_hla_out.extend ([allele1, allele2])
					elif poly_allele2_tmp.startswith (allele2_vb):
						list_hla_out.extend ([allele1, allele2_vb])
					elif poly_allele2_tmp.startswith (allele1_vb):
						list_hla_out.extend ([allele1, allele1_vb])

				# Low confidance allele2
				else:
					print (
						"[WARNING] OptiType allele2({}) not be supported by PolySolver ({}, rank={:.2f}, score_percent={:.2f}) !".format (
							allele2, poly_allele2, rank, score_percent))
					poly_allele2, rank, score_percent = get_rank_score (zip_allele2_scores_gene, allele2_vb)
					if rank < 100 and score_percent > 0.5:
						list_hla_out.extend ([allele1_vb, allele2_vb])
					else:
						list_hla_out.extend ([allele1, allele2])
			#print(list_hla_out)
			write_output (list_hla_out, output)
			#write_mhcpan_input (list_hla_out, output1)
			#write_IEDB_input (list_hla_out, output2)



if __name__ == '__main__':
	main ()
