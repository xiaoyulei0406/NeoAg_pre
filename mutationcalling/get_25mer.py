# Import modules needed for program
import argparse
import Bio
from Bio.Seq import Seq
import sys
import os
import re
import numpy as np
import pandas as pd
from xlsxwriter.workbook import Workbook
import xlrd

import subprocess


def get_pep_SNP(filename, outfile): ##extrace 25mer peptide from SNP annotations
	with open (outfile, 'w') as fw:
		with open (filename, 'r') as fp:
			for readline in fp:
				if '##' in readline:
					continue
				if '#Uploaded_variation' in readline:
					readline = readline.replace ('\n', '')
					ar = readline.strip ().split ('\t')
					ar.insert (len (readline), "25merpeptide")
					fw.write (''.join (readline) + "\t" + "WT_25mer_AA" + "\t" + "MT_25mer_AA" + "\t" + "GeneName" + "\t" + "Mutation"
							  + "\t" + "REF" + "\t" + "loc_mut" + "\t" + "len_mutation" + "\t" + "Start_Position"+ "\t" + "chrom"+ "\t" + "cDNA_Change"+ "\t" + "nm_Feature"
							  + "\t" + "nm_Featuretype" + "\t" + "EXON"+ "\t" + "gnomAD_AF" + "\t" +"1000genome_AF" + "\t" +"ESP_AA_AF" + "\t" +"ESP_EA_AF"+"\t" + "dbSNP_ID" + "\t" + "COSMIC_ID" + "\n")
					continue
				if 'WildtypeProtein' not in readline:
					continue
				if 'intron_variant' in readline:
					continue
				if 'downstream_gene_variant' in readline:
					continue
				if 'regulatory_region_variant' in readline:
					continue
				if 'intergenic_variant' in readline:
					continue
				if 'stop_gained' in readline:
					continue
				if 'stop_lost' in readline:
					continue
				if 'TF_binding_site_variant' in readline:
					continue
				if 'coding_sequence_variant' in readline:
					continue
				if 'upstream_gene_variant' in readline:
					continue
				if 'non_coding_transcript_exon_variant' in readline:
					continue
				if '5_prime_UTR_variant' in readline:
					continue
				if '3_prime_UTR_variant' in readline:
					continue
				if 'splice_acceptor_variant' in readline:
					continue
				if 'splice_donor_variant' in readline:
					continue
				else:
					mypep_seq = readline.split ("WildtypeProtein=", 1)[-1]
					mypep_seq = mypep_seq.replace ('\n', '')
					readline = readline.replace ('\n', '')
					ar = readline.strip ().split ('\t')
					(seq_raw, loc) = get_25mers_SNP_rawseq (ar[9], mypep_seq)

					if 'missense_variant' in (ar[6]):
						##non-synonymous_variant
						aaChange_seq = str (ar[10])[-1]
						(seq_mut, len_mutation) = get_25mers_mut (ar[9], aaChange_seq, mypep_seq)
						##print gene_symbol
						if "SYMBOL" in readline:
							start = ";SYMBOL="
							end = ";SYMBOL_SOURCE="
							symbol = readline[readline.find (start) + len (start):readline.rfind (end)]
						else:
							symbol = "_"
						start_pos_tmp = str(ar[1]).split(":")[1]
						start_pos=start_pos_tmp.split("-")[0]
						chrom=str(ar[1]).split(":")[0]
						cDNA_change = "c." + ar[0].split ("_")[-1].split ("/")[0] + str (ar[8]) + ar[0].split ("_")[-1].split ("/")[1]
						nm_Feature= str(ar[4])
						nm_Featuretype = str (ar[5])
						exon=readline[readline.find (";EXON=") + len (";HGVSc"):readline.rfind (";HGVSc")].split("/")[0]
						if ";gnomAD_AF" in readline:
							start1 = ";gnomAD_AF="
							end1 = ";gnomAD_AFR_AF="
							gnomAD_AF_tmp = readline[readline.find (start1) + len (start):readline.rfind (end1)]
							if ";" in gnomAD_AF_tmp:
								gnomAD_AF = gnomAD_AF_tmp.split (";")[0]
							else:
								gnomAD_AF = gnomAD_AF_tmp
						else:
							gnomAD_AF = "_"
						if ";AF" in readline:
							start1 = ";AF="
							end1 = ";AFR_AF="
							onegenome_AF_tmp = readline[readline.find (start1) + 1:readline.rfind (end1)]
							if "_" in onegenome_AF_tmp:
								onegenome_AF = "_"
							else:
								onegenome_AF = onegenome_AF_tmp
						else:
							onegenome_AF = "_"
						if ";AA_AF=" in readline:
							start1 = ";AA_AF="
							end1 = ";EA_AF="
							ESP_AA_AF = readline[readline.find (start1) + 1:readline.rfind (end1)]
						else:
							ESP_AA_AF = "_"
						if ";EA_AF=" in readline:
							start1 = ";EA_AF="
							end1 = ";gnomAD_AF="
							ESP_EA_AF_tmp = readline[readline.find (start1) + 1:readline.rfind (end1)]
							if ";" in ESP_EA_AF_tmp:
								ESP_EA_AF = ESP_EA_AF_tmp.split (";")[0]
							else:
								ESP_EA_AF = ESP_EA_AF_tmp
						else:
							ESP_EA_AF = "_"
						## seperate dbsnp and cosmic annotation
						anno=get_dbsnp_cosmic(ar[12])
						fw.write ("".join (readline) + "\t" + str (seq_raw) + "\t" + str (seq_mut) + "\t" + str (symbol) +
								  "\t"+ "p." + str (ar[10])[0] + str (ar[9]) + str (ar[10])[-1] + "\t" + ar[0].split ("_")[-1].split ("/")[0] + "\t"
								  + str(loc)+ "\t" + str(len_mutation) + "\t" + str(start_pos)+ "\t"+ str(chrom)+ "\t" + str(cDNA_change)+ "\t" + str(nm_Feature)+"\t" + str(nm_Featuretype)+
								  "\t" + 'exon'+str(exon)+ "\t" + str(gnomAD_AF)+"\t" + str(onegenome_AF)+"\t" + str(ESP_AA_AF)+"\t" + str(ESP_EA_AF)+ "\t"+anno+"\n")
					else:
#						fw.write(''.join(readline) + "\t" + "NULL"+ "\n")
						continue

	fw.close ()


def get_25mers_mut(loc, aaChange_seq, seq_pep):
# This function gets 25mer sequence for frameshift variants
	int_loc = int (loc)
	len_mut = len (aaChange_seq)
	if (len (seq_pep) - int_loc < 12):
		pos_pep = seq_pep[int_loc:len (seq_pep)]
	else:
		pos_pep = seq_pep[int_loc:int_loc + 12]
	if int_loc < 13:
		pre_pep = seq_pep[0:int_loc - 1]
		newloc = int_loc
	else:
		pre_pep = seq_pep[int_loc - 13:int_loc - 1]
		newloc = 13
	aa = Seq (aaChange_seq)
	seq = pre_pep + aa + pos_pep
	return (seq, len_mut)


def get_25mers_SNP_rawseq(loc, seq_pep):
	if "-" not in loc:
		int_loc = int (loc)
	if "-" in loc:
		loc1 = loc.split ("-")[0]
		int_loc = int (loc1)
	if (len (seq_pep) - int_loc < 12):
		pos_pep = seq_pep[int_loc:len (seq_pep)]
	else:
		pos_pep = seq_pep[int_loc:int_loc + 12]
	if int_loc < 13:
		pre_pep = seq_pep[0:int_loc]
		newloc = int_loc
	else:
		pre_pep = seq_pep[int_loc - 13:int_loc]
		newloc = 13

	seq = pre_pep + pos_pep
	return (seq, newloc)

def get_dbsnp_cosmic(readline):
	exit_var = readline.split (",")
	if ("rs" in readline) and ("COSM" in readline):
		dbSNP = []
		cosM = []
		for item in exit_var:
			if item.startswith ("rs"):
				dbSNP.append (item)
			if item.startswith ("COSM"):
				cosM.append (item)
		anno=",".join (str (item) for item in dbSNP) + "\t" + ",".join (str (item) for item in cosM)
	if ("rs" in readline) and ("COSM" not in readline):
		dbSNP = []
		for item in exit_var:
			if item.startswith ("rs"):
				dbSNP.append (item)
		anno= ",".join (str (item) for item in dbSNP) + "\t" + "_"
	if ("rs" not in readline) and ("COSM" in readline):
		cosM = []
		for item in exit_var:
			if item.startswith ("COSM"):
				cosM.append (item)
		anno="_" + "\t" + ",".join (str (item) for item in cosM)
	if ("rs" not in readline) and ("COSM" not in readline):
		anno="_" + "\t" + "_"
	return(anno)

def get_map(tmpfile,annovar_file,outfile,type,outdir,tumor_ID,VEP_c):

	##Read vep outputfile into a data frame:
	vepOutfile = pd.read_table (tmpfile, sep='\t')

	##Read ANNOVAR MAF file into a data frame:
	annovar_FILE = pd.read_table (annovar_file, sep='\t')
	##Obtain several columns from raw vep output file
	vepOutfile_1 = vepOutfile[
		['#Uploaded_variation', 'Location', 'chrom', 'Start_Position', 'REF', 'Allele', 'GeneName', 'Mutation',
		 'Consequence', 'WT_25mer_AA','MT_25mer_AA', 'loc_mut', 'len_mutation', 'cDNA_Change', 'Feature', 'Feature_type',
		 'EXON', 'gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF", 'dbSNP_ID','COSMIC_ID']]


	##Rename column
	vepOutfile_1.columns = ['Variant_Key', 'ID', 'Chr', 'Start', 'Ref', 'Alt', 'Hugo_Symbol', 'AAchange', 'Consequence',
							'WT_25mer_AA', 'MT_25mer_AA', 'loc_mut', 'len_mutation','cDNA_Change', 'Feature', 'Feature_type', 'EXON',
							'gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",'dbSNP_ID','COSMIC_ID']
	vepOutfile_uniq = vepOutfile_1.drop_duplicates ()
	vepOutfile_uniq['Chr'] = vepOutfile_uniq['Chr'].astype (str)
	annovar_FILE = annovar_FILE[['Chromosome', "Start_Position", "End_Position", 'Reference_Allele', 'Tumor_Seq_Allele2',
								 'ExonicFunc.refGene', 'Callers', 'NCallers', 'Tumor_Sample_Barcode','tumor_vaf', 'tumor_ref_count',
								 'tumor_alt_count', 'normal_vaf', 'normal_ref_count','normal_alt_count']]
	#    print(vepOutfile_uniq)
	annovar_FILE.columns = ['Chr', "Start", "End", 'Ref_raw', 'Alt_raw', 'ExonicFunc.refGene', 'Callers', 'NCallers',
							'Tumor_Sample_Barcode','tumor_vaf', 'tumor_ref_count',
							'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count']
	annovar_FILE['Chr'] = annovar_FILE['Chr'].astype (str)
	df_new = pd.merge (vepOutfile_uniq, annovar_FILE, how='inner', on=['Chr', 'Start'])
	#    print(df_new)
	df_final = df_new[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'Feature', 'ExonicFunc.refGene', 'EXON', 'Consequence',
		 'cDNA_Change', 'AAchange', 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA',
		 'loc_mut', 'len_mutation','gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",'dbSNP_ID','COSMIC_ID',
		 'Tumor_Sample_Barcode', 'tumor_vaf', 'tumor_ref_count','tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count']]
	df_final['TYPE']=type
	df_final_uniq = df_final.drop_duplicates ()
	#    print(df_final_uniq)
	df_final_uniq = df_final_uniq[df_final_uniq['Consequence'] != "synonymous_variant"]

	df_selcols = df_final_uniq[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'Feature', 'ExonicFunc.refGene', 'EXON', 'Consequence',
		 'cDNA_Change', 'AAchange', 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA',
		 'loc_mut', 'len_mutation','gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",'dbSNP_ID','COSMIC_ID',
		 'Tumor_Sample_Barcode', 'TYPE','tumor_vaf', 'tumor_ref_count','tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count']]
	df_selcols = df_selcols.drop_duplicates ()

	df_selcols['RefGene.feature'] = df_selcols[['Feature', 'EXON', 'cDNA_Change', 'AAchange']].apply (
		lambda x: ':'.join (x.astype (str)), axis=1)

	df_merge = df_selcols[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'RefGene.feature', 'ExonicFunc.refGene', 'Consequence',
		 'Callers', 'NCallers', 'gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",'dbSNP_ID','COSMIC_ID',
		 'WT_25mer_AA', 'MT_25mer_AA', 'loc_mut','len_mutation', 'Tumor_Sample_Barcode', 'TYPE','tumor_vaf', 'tumor_ref_count',
		 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count']]
	# merge AAchange in gene that have same 25mer sequences
	df_counts_merge = df_merge.groupby (['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'ExonicFunc.refGene', 'Consequence',
		'Callers', 'NCallers','gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF", 'dbSNP_ID','COSMIC_ID','WT_25mer_AA', 'MT_25mer_AA',
		 'loc_mut','len_mutation', 'Tumor_Sample_Barcode', 'TYPE','tumor_vaf', 'tumor_ref_count','tumor_alt_count', 'normal_vaf', 'normal_ref_count',
		'normal_alt_count'], as_index=False).agg ({'RefGene.feature': ','.join})
	print (df_counts_merge.shape[0])
	##change column order
	'''
	df_counts_merge = df_counts_merge[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt','RefGene.feature','Consequence',
		'ExonicFunc.refGene', 'tumor_vaf', 'tumor_ref_count','tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
		'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA','COSMIC_ID','dbSNP_ID','gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",
		'Tumor_Sample_Barcode', 'TYPE','loc_mut','len_mutation']]
	'''
	df_counts_merge = df_counts_merge.drop_duplicates ()

	##get neighbor base within 36
	df_counts_merge['neighbor_distance'] = df_counts_merge.groupby ('Chr')['Start'].apply (lambda x: x - x.shift (1))
	#print(df_counts_merge)
	df_counts_merge['neighbor_distance'] = df_counts_merge['neighbor_distance'].replace (np.nan, "0")
	df_counts_merge['neighbor_distance'] = df_counts_merge['neighbor_distance'].astype(int)
	#print(df_counts_merge)
	df_counts_merge.loc[(df_counts_merge['neighbor_distance'] > 36),'neighbor_distance'] = "0"
	#print(df_counts_merge)
	df_counts_merge['neighbor_distance'] = df_counts_merge['neighbor_distance'].astype(int)
	df_counts_merge.loc[(df_counts_merge['neighbor_distance'] < 0), 'neighbor_distance'] = "0"
	df_counts_merge.loc[(df_counts_merge['neighbor_distance'] == 0), 'neighbor_distance'] = "0"
	df_counts_merge['neighbor_distance'] = df_counts_merge['neighbor_distance'].astype(int)
	df_counts_merge.loc[(df_counts_merge['neighbor_distance'] == 0), 'neighbor_distance'] = "_"
	#print(df_counts_merge)

	for i in range(0, len(df_counts_merge)):
		if (df_counts_merge.loc[i, 'neighbor_distance'] != '_'):
			df_counts_merge.loc[i, 'neighbor_variants'] = df_counts_merge.loc[i - 1, 'ID']
		elif (df_counts_merge.loc[i, 'neighbor_distance'] == '_'):
			df_counts_merge.loc[i, 'neighbor_variants'] = "_"
	df_counts_merge=df_counts_merge.replace (np.nan, "_")

		##change column order
	df_counts_merge = df_counts_merge[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'RefGene.feature', 'Consequence',
			 'ExonicFunc.refGene', 'tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
			 'normal_alt_count',
			 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA', 'COSMIC_ID', 'dbSNP_ID', 'gnomAD_AF', "1000genome_AF",
			 "ESP_AA_AF", "ESP_EA_AF",
			 'Tumor_Sample_Barcode', 'TYPE', 'neighbor_variants','loc_mut', 'len_mutation']]
	df_counts_merge.to_csv (outfile, index=False, sep="\t")
	##add cannonical information
	df_1=get_canonical(df_counts_merge,VEP_c)
	df_1=df_1.drop_duplicates()
	##change column order
	df_1 = df_1[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'RefGene.feature', 'carnonical_info', 'Consequence',
		 'ExonicFunc.refGene', 'tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
		 'normal_alt_count',
		 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA', 'COSMIC_ID', 'dbSNP_ID', 'gnomAD_AF', "1000genome_AF",
		 "ESP_AA_AF", "ESP_EA_AF",
		 'Tumor_Sample_Barcode', 'TYPE', 'neighbor_variants', 'loc_mut', 'len_mutation']]
	#print(df_1)
	out_txt = outdir + '/' + tumor_ID + "." + type + ".uniq.25mer.txt"
	df_1.to_csv (out_txt, index=False, sep="\t")

	out_highlight=outdir+'/'+tumor_ID + "."+ type +".highlight.xlsx"
#	out1=get_highlight(df_counts_merge,out_highlight,type)
	df_Var = df_1[["Chr", "Start", "End", 'Ref', 'Alt', 'Hugo_Symbol','Consequence','tumor_vaf', 'normal_vaf', 'tumor_ref_count',
				  'tumor_alt_count',  'normal_ref_count','normal_alt_count','Callers', 'NCallers','COSMIC_ID',
				  'dbSNP_ID', 'gnomAD_AF', "1000genome_AF","ESP_AA_AF", "ESP_EA_AF",'Tumor_Sample_Barcode', 'TYPE']]

	##get special columns
	df_uniqVar = df_Var.drop_duplicates()

	out2= outdir + '/' + tumor_ID + "." + type + ".vcf4pyclone.txt"
	df_uniqVar.to_csv(out2,index=False, sep="\t")
	out1=get_highlight(df_1,out_highlight,type)

#	out2=outdir+'/'+tumor_ID + ".raw1.xlsx"
#	df_counts_merge.to_excel (out2, index=None, header=True)

def get_map_MNP(tmpfile,annovar_file,outfile,type,outdir,tumor_ID,VEP_c):

	##Read vep outputfile into a data frame:
	vepOutfile = pd.read_table (tmpfile, sep='\t')

	##Read ANNOVAR MAF file into a data frame:
	annovar_FILE = pd.read_table (annovar_file, sep='\t')
	##Obtain several columns from raw vep output file
	vepOutfile_1 = vepOutfile[
		['#Uploaded_variation', 'Location', 'chrom', 'Start_Position', 'REF', 'Allele', 'GeneName', 'Mutation',
		 'Consequence', 'WT_25mer_AA','MT_25mer_AA', 'loc_mut', 'len_mutation', 'cDNA_Change', 'Feature', 'Feature_type',
		 'EXON', 'gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF", 'dbSNP_ID','COSMIC_ID']]


	##Rename column
	vepOutfile_1.columns = ['Variant_Key', 'ID', 'Chr', 'Start', 'Ref', 'Alt', 'Hugo_Symbol', 'AAchange', 'Consequence',
							'WT_25mer_AA', 'MT_25mer_AA', 'loc_mut', 'len_mutation','cDNA_Change', 'Feature', 'Feature_type', 'EXON',
							'gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",'dbSNP_ID','COSMIC_ID']
	vepOutfile_uniq = vepOutfile_1.drop_duplicates ()
	vepOutfile_uniq['Chr'] = vepOutfile_uniq['Chr'].astype (str)
	annovar_FILE = annovar_FILE[['Chromosome', "Start_Position", "End_Position", 'Reference_Allele', 'Tumor_Seq_Allele2',
								 'ExonicFunc.refGene', 'Callers', 'NCallers', 'Tumor_Sample_Barcode','tumor_vaf', 'tumor_ref_count',
								 'tumor_alt_count', 'normal_vaf', 'normal_ref_count','normal_alt_count']]
	#    print(vepOutfile_uniq)
	annovar_FILE.columns = ['Chr', "Start", "End", 'Ref_raw', 'Alt_raw', 'ExonicFunc.refGene', 'Callers', 'NCallers',
							'Tumor_Sample_Barcode','tumor_vaf', 'tumor_ref_count',
							'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count']
	annovar_FILE['Chr'] = annovar_FILE['Chr'].astype (str)
	df_new = pd.merge (vepOutfile_uniq, annovar_FILE, how='inner', on=['Chr', 'Start'])
	#    print(df_new)
	df_final = df_new[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'Feature', 'ExonicFunc.refGene', 'EXON', 'Consequence',
		 'cDNA_Change', 'AAchange', 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA',
		 'loc_mut', 'len_mutation','gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",'dbSNP_ID','COSMIC_ID',
		 'Tumor_Sample_Barcode', 'tumor_vaf', 'tumor_ref_count','tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count']]
	df_final['TYPE']=type
	df_final_uniq = df_final.drop_duplicates ()
	#    print(df_final_uniq)
	df_final_uniq = df_final_uniq[df_final_uniq['Consequence'] != "synonymous_variant"]

	df_selcols = df_final_uniq[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'Feature', 'ExonicFunc.refGene', 'EXON', 'Consequence',
		 'cDNA_Change', 'AAchange', 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA',
		 'loc_mut', 'len_mutation','gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",'dbSNP_ID','COSMIC_ID',
		 'Tumor_Sample_Barcode', 'TYPE','tumor_vaf', 'tumor_ref_count','tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count']]
	df_selcols = df_selcols.drop_duplicates ()

	df_selcols['RefGene.feature'] = df_selcols[['Feature', 'EXON', 'cDNA_Change', 'AAchange']].apply (
		lambda x: ':'.join (x.astype (str)), axis=1)

	df_merge = df_selcols[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'RefGene.feature', 'ExonicFunc.refGene', 'Consequence',
		 'Callers', 'NCallers', 'gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",'dbSNP_ID','COSMIC_ID',
		 'WT_25mer_AA', 'MT_25mer_AA', 'loc_mut','len_mutation', 'Tumor_Sample_Barcode', 'TYPE','tumor_vaf', 'tumor_ref_count',
		 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count']]
	# merge AAchange in gene that have same 25mer sequences
	df_counts_merge = df_merge.groupby (['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'ExonicFunc.refGene', 'Consequence',
		'Callers', 'NCallers','gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF", 'dbSNP_ID','COSMIC_ID','WT_25mer_AA', 'MT_25mer_AA',
		 'loc_mut','len_mutation', 'Tumor_Sample_Barcode', 'TYPE','tumor_vaf', 'tumor_ref_count','tumor_alt_count', 'normal_vaf', 'normal_ref_count',
		'normal_alt_count'], as_index=False).agg ({'RefGene.feature': ','.join})
	#print (df_counts_merge.shape[0])
	##change column order
	'''
	df_counts_merge = df_counts_merge[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt','RefGene.feature','Consequence',
		'ExonicFunc.refGene', 'tumor_vaf', 'tumor_ref_count','tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
		'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA','COSMIC_ID','dbSNP_ID','gnomAD_AF', "1000genome_AF", "ESP_AA_AF", "ESP_EA_AF",
		'Tumor_Sample_Barcode', 'TYPE','loc_mut','len_mutation']]
	'''
	df_counts_merge = df_counts_merge.drop_duplicates ()

	##get neighbor base within 36
	if (df_counts_merge.shape[0]==1):
		df_counts_merge['neighbor_variants'] = "_"
	elif (df_counts_merge.shape[0] > 1):
		df_counts_merge['neighbor_distance'] = df_counts_merge.groupby ('Chr')['Start'].apply (lambda x: x - x.shift (1))
		# print(dfCon)
#		print(df_counts_merge)
		df_counts_merge['neighbor_distance'] = df_counts_merge['neighbor_distance'].replace (np.nan, "0")
		df_counts_merge['neighbor_distance'] = df_counts_merge['neighbor_distance'].astype(int)
		#print(df_counts_merge)
		df_counts_merge.loc[(df_counts_merge['neighbor_distance'] > 36),'neighbor_distance'] = "0"
		#print(df_counts_merge)
		df_counts_merge['neighbor_distance'] = df_counts_merge['neighbor_distance'].astype(int)
		df_counts_merge.loc[(df_counts_merge['neighbor_distance'] < 0), 'neighbor_distance'] = "0"
		df_counts_merge.loc[(df_counts_merge['neighbor_distance'] == 0), 'neighbor_distance'] = "0"
		df_counts_merge['neighbor_distance'] = df_counts_merge['neighbor_distance'].astype(int)
		df_counts_merge.loc[(df_counts_merge['neighbor_distance'] == 0), 'neighbor_distance'] = "_"

		for i in range(1, len(df_counts_merge)):
			if (df_counts_merge.loc[i, 'neighbor_distance'] != '_'):
				df_counts_merge.loc[i, 'neighbor_variants'] = df_counts_merge.loc[i - 1, 'ID']
			else:
				df_counts_merge.loc[i, 'neighbor_variants'] = "_"
	df_counts_merge=df_counts_merge.replace (np.nan, "_")
		##change column order
	df_counts_merge = df_counts_merge[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'RefGene.feature', 'Consequence',
			 'ExonicFunc.refGene', 'tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
			 'normal_alt_count',
			 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA', 'COSMIC_ID', 'dbSNP_ID', 'gnomAD_AF', "1000genome_AF",
			 "ESP_AA_AF", "ESP_EA_AF",
			 'Tumor_Sample_Barcode', 'TYPE', 'neighbor_variants','loc_mut', 'len_mutation']]
	df_counts_merge.to_csv (outfile, index=False, sep="\t")
	##add cannonical information
	df_1=get_canonical(df_counts_merge,VEP_c)
	df_1=df_1.drop_duplicates()
	##change column order
	df_1 = df_1[['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'RefGene.feature', 'carnonical_info', 'Consequence',
		 'ExonicFunc.refGene', 'tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
		 'normal_alt_count',
		 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA', 'COSMIC_ID', 'dbSNP_ID', 'gnomAD_AF', "1000genome_AF",
		 "ESP_AA_AF", "ESP_EA_AF",
		 'Tumor_Sample_Barcode', 'TYPE', 'neighbor_variants', 'loc_mut', 'len_mutation']]
	#print(df_1)
	out_txt = outdir + '/' + tumor_ID + "." + type + ".uniq.25mer.txt"
	df_1.to_csv (out_txt, index=False, sep="\t")

	out_highlight=outdir+'/'+tumor_ID + "."+ type +".highlight.xlsx"
#	out1=get_highlight(df_counts_merge,out_highlight,type)
	df_Var = df_1[["Chr", "Start", "End", 'Ref', 'Alt', 'Hugo_Symbol','Consequence','tumor_vaf', 'normal_vaf', 'tumor_ref_count',
				  'tumor_alt_count',  'normal_ref_count','normal_alt_count','Callers', 'NCallers','COSMIC_ID',
				  'dbSNP_ID', 'gnomAD_AF', "1000genome_AF","ESP_AA_AF", "ESP_EA_AF",'Tumor_Sample_Barcode', 'TYPE']]

	##get special columns
	df_uniqVar = df_Var.drop_duplicates()

	out2= outdir + '/' + tumor_ID + "." + type + ".vcf4pyclone.txt"
	df_uniqVar.to_csv(out2,index=False, sep="\t")
	out1=get_highlight(df_1,out_highlight,type)

def get_canonical(df_counts_merge,VEP_c):
	'''

	:param df_counts_merge:
	:param VEP_c:
	:param out_file_c:
	:param type:
	:return:
	'''
	df=pd.read_table (VEP_c, sep='\t',skiprows=38)
	new_header=['Uploaded_variation','Location','Allele','Gene','Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','Extra','CANONICAL']
	df.columns = new_header
	df= df.drop_duplicates (subset='Feature', keep="first")
	df[['refAA','refP']] = df.Amino_acids.str.split ("/",expand=True)
	df['aaName']='p.'
	df['carnonical_info']= df[['aaName','refAA','Protein_position','refP']].apply (lambda x: ''.join (x.astype (str)), axis=1)
	#select out new columns
	df_filtered = df[df['CANONICAL'] == 'YES']
	df_filtered = df_filtered[['Location','carnonical_info']]
	df_filtered.columns = ['ID','carnonical_info']
	df_filtered = df_filtered[df_filtered['carnonical_info'] != 'p.--None']
	df_filtered=df_filtered.drop_duplicates(subset='ID', keep="first")

	df_merge= pd.merge (df_counts_merge, df_filtered, how='left', on=['ID'])
	df_merge = df_merge.replace (np.nan, "_")
	df_merge= df_merge.drop_duplicates ()
	return(df_merge)



def get_highlight(df_counts_merge,outfile,type):
	workbook = Workbook (outfile)
	worksheet = workbook.add_worksheet ()
	red = workbook.add_format ({'color': 'red', 'bold': True, 'font_size': '12'})
	blue = workbook.add_format ({'color': 'blue', 'bold': True, 'font_size': '12'})
	for col in range (50):
		worksheet.set_column (col, col, 10)
	df1 = df_counts_merge
	# Add a bold format to use to highlight cells.
	bold = workbook.add_format ({'bold': True})
	# Write data headers.

	headers=['ID', 'Hugo_Symbol', "Chr", "Start", "End", 'Ref', 'Alt', 'RefGene.feature','carnonical_info', 'Consequence',
			 'ExonicFunc.refGene', 'tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
			 'normal_alt_count',
			 'Callers', 'NCallers', 'WT_25mer_AA', 'MT_25mer_AA', 'COSMIC_ID', 'dbSNP_ID', 'gnomAD_AF', "1000genome_AF",
			 "ESP_AA_AF", "ESP_EA_AF",
			 'Tumor_Sample_Barcode', 'TYPE', 'neighbor_variants']
	worksheet.write_row('A1', headers, bold)


	wildSeqs = df1['WT_25mer_AA']
	mutSeqs = df1['MT_25mer_AA']
	locMuts = df1['loc_mut']
	len_mutation = df1['len_mutation']
	workbook.use_zip64 ()

	for rowNum, data in enumerate (df1['ID']):
		worksheet.write (rowNum+1, 0, data)
	for rowNum, data in enumerate (df1['Hugo_Symbol']):
		worksheet.write (rowNum+1, 1, data)
	for rowNum, data in enumerate (df1['Chr']):
		worksheet.write (rowNum+1, 2, data)
	for rowNum, data in enumerate (df1['Start']):
		worksheet.write (rowNum+1, 3, data)
	for rowNum, data in enumerate (df1['End']):
		worksheet.write (rowNum+1, 4, data)
	for rowNum, data in enumerate (df1['Ref']):
		worksheet.write (rowNum+1, 5, data)
	for rowNum, data in enumerate (df1['Alt']):
		worksheet.write (rowNum+1, 6, data)
	for rowNum, data in enumerate (df1['RefGene.feature']):
		worksheet.write (rowNum+1, 7, data)
	for rowNum, data in enumerate (df1['carnonical_info']):
		worksheet.write (rowNum + 1, 8, data)
	for rowNum, data in enumerate (df1['Consequence']):
		worksheet.write (rowNum+1, 9, data)
	for rowNum, data in enumerate (df1['ExonicFunc.refGene']):
		worksheet.write (rowNum+1, 10, data)
	for rowNum, data in enumerate (df1['tumor_vaf']):
		worksheet.write (rowNum+1, 11, data)
	for rowNum, data in enumerate (df1['tumor_ref_count']):
		worksheet.write (rowNum+1, 12, data)
	for rowNum, data in enumerate (df1['tumor_alt_count']):
		worksheet.write (rowNum+1, 13, data)
	for rowNum, data in enumerate (df1['normal_vaf']):
		worksheet.write (rowNum+1, 14, data)
	for rowNum, data in enumerate (df1['normal_ref_count']):
		worksheet.write (rowNum+1, 15, data)
	for rowNum, data in enumerate (df1['normal_alt_count']):
		worksheet.write (rowNum+1, 16, data)
	for rowNum, data in enumerate (df1['Callers']):
		worksheet.write (rowNum+1, 17, data)
	for rowNum, data in enumerate (df1['NCallers']):
		worksheet.write (rowNum+1, 18, data)
	for rowNum, sequence in enumerate (wildSeqs):
		highwildSeq = []
		i = 0
		for base in sequence:
			if (i == int (locMuts[rowNum]) - 1):
				highwildSeq.extend ((blue, base))
			else:
				highwildSeq.extend (base)
			i = i + 1
		worksheet.write_rich_string (rowNum+1, 19, *highwildSeq)

	for rowNum, sequence in enumerate (mutSeqs):
		highwildSeq = []
		j = 0
		flag=0
		for base in sequence:
			if ( type =="INDEL" or type =="SNP"):
#			if len_mutation[rowNum]==1:
				if (j == int (locMuts[rowNum]) - 1):
					highwildSeq.extend ((red, base))
				else:
					highwildSeq.extend (base)
				j = j + 1
#			elif len_mutation[rowNum]>1:
			elif (type =="MNP"):
				if len_mutation[rowNum] == 1:
					if (j == int (locMuts[rowNum]) - 1):
						highwildSeq.extend ((red, base))
					else:
						highwildSeq.extend (base)
					j = j + 1
				elif len_mutation[rowNum] > 1:
					if flag < len_mutation[rowNum]:
						if (j >= int (locMuts[rowNum]) - 1):
							highwildSeq.extend ((red, base))
							flag=flag+1
							#print("flag:" + str(flag))
						else:
							highwildSeq.extend (base)
				else:
					highwildSeq.extend (base)
					j = j + 1
			else:
				pass
		flag=0
		worksheet.write_rich_string (rowNum+1, 20, *highwildSeq)
	for rowNum, data in enumerate (df1['COSMIC_ID']):
		worksheet.write (rowNum+1, 21, data)
	for rowNum, data in enumerate (df1['dbSNP_ID']):
		worksheet.write (rowNum+1, 22, data)
	for rowNum, data in enumerate (df1['gnomAD_AF']):
		worksheet.write (rowNum+1, 23, data)
	for rowNum, data in enumerate (df1['1000genome_AF']):
		worksheet.write (rowNum+1, 24, data)
	for rowNum, data in enumerate (df1['ESP_AA_AF']):
		worksheet.write (rowNum+1, 25, data)
	for rowNum, data in enumerate (df1['ESP_EA_AF']):
		worksheet.write (rowNum+1, 26, data)
	for rowNum, data in enumerate (df1['Tumor_Sample_Barcode']):
		worksheet.write (rowNum+1, 27, data)
	for rowNum, data in enumerate (df1['TYPE']):
		worksheet.write (rowNum+1, 28, data)
	#print(1)
	for rowNum, data in enumerate (df1['neighbor_variants']):
		worksheet.write (rowNum+1, 29, data)


	workbook.close ()
	return(outfile)
'''
def map_highlight(out1,hightmp,outdir,tumor_ID,outfile):

	pd.set_option ('display.float_format', '{:.2g}'.format)
	highlight_cols = pd.read_excel (hightmp, index =0, usecols="A:B", sheet_name="Sheet1")
	df = pd.read_excel (out1, index =0, sheet_name="Sheet1")
	highlight_cols.columns = ['WT_hl', 'MT_hl']

	df['WT_hl'] = highlight_cols[['WT_hl']]
	df['MT_hl'] = highlight_cols[['MT_hl']]
	df.to_excel (outfile, sheet_name='report_1')
'''

def get_pep_INDEL(filename, outfile): ##extrace 25mer peptide from INDEL annotations
	with open (outfile, 'w') as fw:
		with open (filename, 'r') as fp:
			for readline in fp:
				if '##' in readline:
					continue
				if '#Uploaded_variation' in readline:
					readline = readline.replace ('\n', '')
					ar = readline.strip ().split ('\t')
					ar.insert (len (readline), "25merpeptide")
					fw.write (''.join (readline) + "\t" + "WT_25mer_AA" + "\t" + "MT_25mer_AA" + "\t" + "GeneName" + "\t" + "Mutation"
							  + "\t" + "REF" + "\t" + "loc_mut" + "\t" + "len_mutation" + "\t" + "Start_Position"+ "\t" + "chrom"+ "\t" + "cDNA_Change"+ "\t" + "nm_Feature"
							  + "\t" + "nm_Featuretype" + "\t" + "EXON"+ "\t" + "gnomAD_AF" +"\t" +"1000genome_AF" + "\t" +"ESP_AA_AF" + "\t" +"ESP_EA_AF"
							  +"\t" + "dbSNP_ID" + "\t" + "COSMIC_ID" + "\n")
					continue
				if 'WildtypeProtein' not in readline:
					continue
				if 'intron_variant' in readline:
					continue
				if 'downstream_gene_variant' in readline:
					continue
				if 'regulatory_region_variant' in readline:
					continue
				if 'intergenic_variant' in readline:
					continue
				if 'splice_region_variant' in readline:
					continue
				if 'stop_gained' in readline:
					continue
				if 'stop_lost' in readline:
					continue
				if 'TF_binding_site_variant' in readline:
					continue
				if 'coding_sequence_variant' in readline:
					continue
				if 'upstream_gene_variant' in readline:
					continue
				if 'non_coding_transcript_exon_variant' in readline:
					continue
				if '5_prime_UTR_variant' in readline:
					continue
				if '3_prime_UTR_variant' in readline:
					continue
#				if 'splice_acceptor_variant':
#					continue
#				if 'splice_donor_variant':
#					continue
#				if 'start_lost':
#					continue
				else:
					mypep_seq = readline.split ("WildtypeProtein=", 1)[-1]
					mypep_seq = mypep_seq.replace ('\n', '')
					readline = readline.replace ('\n', '')
					ar = readline.strip ().split ('\t')
					(seq_raw,loc) = get_25mers_indel_rawseq(ar[9], ar[10], mypep_seq)
					if 'frameshift_variant'  in (ar[6]):
						if "frameshift_variant,splice_region_variant" in (ar[6]):
							continue
						aaChange_seq = str (ar[10]).split("/")[1]
						aaRaw_seq = str (ar[10]).split ("/")[0]
						downstream_seq_tmp = readline.split ("DownstreamProtein=")[-1]
						downstream_seq = downstream_seq_tmp.split (";ProteinLengthChange")[0]
						(seq_mut, len_mutation) = get_frameshift(ar[9], aaRaw_seq,aaChange_seq, mypep_seq, downstream_seq)
					if 'inframe_insertion' in (ar[6]):
						aaChange_seq = str (ar[10]).split ("/")[1]
						aaRaw_seq = str (ar[10]).split ("/")[0]
						(seq_mut, len_mutation) = get_25mers_inframe_insertion (ar[9], aaRaw_seq,aaChange_seq, mypep_seq)
					#	print (aaRaw_seq)
					if 'inframe_deletion' in (ar[6]):
						aaChange_seq = str (ar[10]).split ("/")[1]
						aaRaw_seq = str (ar[10]).split ("/")[0]
						(seq_mut, len_mutation) = get_25mers_inframe_deletion (ar[9], aaRaw_seq, aaChange_seq, mypep_seq)
					if "SYMBOL" in readline:
						start = ";SYMBOL="
						end = ";SYMBOL_SOURCE="
						symbol = readline[readline.find (start) + len (start):readline.rfind (end)]
					else:
						symbol = "_"
					start_pos_tmp = str(ar[1]).split(":")[1]
					start_pos=start_pos_tmp.split("-")[0]
					chrom=str(ar[1]).split(":")[0]
					cDNA_change = "c." + ar[0].split ("_")[-1].split ("/")[0] + str (ar[8]) + ar[0].split ("_")[-1].split ("/")[1]
					nm_Feature= str(ar[4])
					nm_Featuretype = str (ar[5])
					exon=readline[readline.find (";EXON=") + len (";HGVSc"):readline.rfind (";HGVSc")].split("/")[0]
					if "gnomAD_AF" in readline:
						start1 = ";gnomAD_AF="
						end1 = ";gnomAD_AFR_AF="
						gnomAD_AF_tmp = readline[readline.find (start1) + len (start):readline.rfind (end1)]
						if ";" in gnomAD_AF_tmp:
							gnomAD_AF = gnomAD_AF_tmp.split (";")[0]
						else:
							gnomAD_AF = gnomAD_AF_tmp
					else:
						gnomAD_AF = "_"
					if ";AF" in readline:
						start1 = ";AF="
						end1 = ";AFR_AF="
						onegenome_AF_tmp = readline[readline.find (start1) + 1:readline.rfind (end1)]
						if "_" in onegenome_AF_tmp:
							onegenome_AF="_"
						else:
							onegenome_AF=onegenome_AF_tmp
					else:
						onegenome_AF = "_"
					if ";AA_AF=" in readline:
						start1 = ";AA_AF="
						end1 = ";EA_AF="
						ESP_AA_AF = readline[readline.find (start1) + 1:readline.rfind (end1)]
					else:
						ESP_AA_AF = "_"
					if ";EA_AF=" in readline:
						start1 = ";EA_AF="
						end1 = ";gnomAD_AF="
						ESP_EA_AF_tmp = readline[readline.find (start1) + 1:readline.rfind (end1)]
						if ";" in ESP_EA_AF_tmp:
							ESP_EA_AF = ESP_EA_AF_tmp.split (";")[0]
						else:
							ESP_EA_AF = ESP_EA_AF_tmp
					else:
						ESP_EA_AF = "_"
					anno = get_dbsnp_cosmic (ar[12])
					fw.write ("".join (readline) + "\t" + str (seq_raw) + "\t" + str (seq_mut) + "\t" + str (symbol) +
								  "\t"+ "p." + str (ar[10]).split ("/")[0] + str (ar[9]) + str(ar[10]).split("/")[1] + "\t" + ar[0].split ("_")[-1].split ("/")[0] + "\t"
								  + str(loc)+ "\t" + str(len_mutation) + "\t" + str(start_pos)+ "\t"+ str(chrom)+ "\t" + str(cDNA_change)+ "\t" + str(nm_Feature)+"\t" + str(nm_Featuretype)+
								  "\t" + 'exon'+str(exon)+ "\t" + str(gnomAD_AF)+"\t" + str(onegenome_AF)+"\t" + str(ESP_AA_AF)+"\t" + str(ESP_EA_AF)+"\t"+anno+"\n")

	fw.close ()


def get_25mers_indel_rawseq(loc, aaChange_seq, seq_pep):
	if "-" not in loc:
		int_loc = int (loc)
	if "-" in loc:
		loc1 = loc.split ("-")[0]
		int_loc = int (loc1)
	if (len (seq_pep) - int_loc < 12):
		pos_pep = seq_pep[int_loc:len (seq_pep)]
	else:
		pos_pep = seq_pep[int_loc:int_loc + 12]
	if int_loc < 13:
		pre_pep = seq_pep[0:int_loc]
		newloc = int_loc
	else:
		pre_pep = seq_pep[int_loc - 13:int_loc]
		newloc = 13
	seq = pre_pep + pos_pep
	return (seq, newloc)

def get_25mers_inframe_insertion(loc, aaRaw_seq,aaChange_seq, seq_pep):
# This function gets 25mer sequence for frameshift variants
	len_mut=len(aaChange_seq)
	if "-" not in loc:
		int_loc = int (loc)
	if "-" in loc:
		loc1 = loc.split ("-")[0]
		int_loc = int (loc1)
	if "-" in aaRaw_seq:
		if int_loc < 13:
			pre_pep = seq_pep[0:int_loc]
			newloc = int_loc
		else:
			pre_pep = seq_pep[int_loc - 13:int_loc]
			newloc = 13
		if (len (seq_pep) - int_loc < 12):
			pos_pep = seq_pep[int_loc:len (seq_pep)]
		else:
			pos_pep = seq_pep[int_loc:(int_loc + 12)]
		seq = pre_pep + aaChange_seq + pos_pep
	if "-" not in aaRaw_seq:
		if int_loc < 13:
			pre_pep = seq_pep[0:int_loc - 1]
			newloc = int_loc
		else:
			pre_pep = seq_pep[int_loc - 13:int_loc - 1]
			newloc = 13
		if (len (seq_pep) - int_loc < 12):
			pos_pep = seq_pep[int_loc:len (seq_pep)]
		else:
			pos_pep = seq_pep[int_loc:(int_loc + 12)]
		seq = pre_pep + aaChange_seq + pos_pep
	return (seq, len_mut)


def get_25mers_inframe_deletion(loc, aaRaw_seq, aaChange_seq, seq_pep):
	aa_rawseq = Seq (aaRaw_seq)
	len_aa=len(aa_rawseq)

	if "-" not in loc:
		int_loc = int (loc)
	if "-" in loc:
		loc1=loc.split("-")[0]
		int_loc = int (loc1)
	if int_loc < 13:
		pre_pep = seq_pep[0:int_loc - 1]
		newloc = int_loc
	else:
		pre_pep = seq_pep[int_loc - 13:int_loc - 1]
		newloc = 13
	if (len (seq_pep) - int_loc < 12):
		pos_pep = seq_pep[(int_loc+len_aa-1):len (seq_pep)]
	else:
		pos_pep = seq_pep[(int_loc+len_aa-1):(int_loc + 12)]
	seq = pre_pep + aaChange_seq + pos_pep
	return (seq, newloc)

def get_frameshift(loc, aaRaw_seq,aaChange_seq, seq_pep, downstream_seq):# This function gets 25mer sequence for frameshift variants
	pos_pep = Seq (downstream_seq)
	aa_rawseq = Seq (aaRaw_seq)
	len_pos=len(pos_pep)
	if "-" not in loc:
		int_loc = int (loc)
	if "-" in loc:
		loc1=loc.split("-")[0]
		int_loc = int (loc1)

	if "-" in aaRaw_seq:
		if int_loc < 13:
			pre_pep = seq_pep[0:int_loc]
			newloc = int_loc
		else:
			pre_pep = seq_pep[int_loc - 13:int_loc]
			newloc = 13
		seq = pre_pep + pos_pep
	if "-" not in aaRaw_seq:
		if int_loc < 13:
			pre_pep = seq_pep[0:int_loc - 1]
			newloc = int_loc
		else:
			pre_pep = seq_pep[int_loc - 13:int_loc - 1]
			newloc = 13
		seq = pre_pep + pos_pep
	return (seq, len_pos)

def get_pep_MNP(filename, outfile):##extrace 25mer peptide from MNP annotations
	with open (outfile, 'w') as fw:
		with open (filename, 'r') as fp:
			for readline in fp:
				if '##' in readline:
					continue
				if '#Uploaded_variation' in readline:
					readline = readline.replace ('\n', '')
					ar = readline.strip ().split ('\t')
					ar.insert (len (readline), "25merpeptide")
					fw.write (''.join (readline) + "\t" + "WT_25mer_AA" + "\t" + "MT_25mer_AA" + "\t" + "GeneName" + "\t" + "Mutation"
							  + "\t" + "REF" + "\t" + "loc_mut" + "\t" + "len_mutation" + "\t" + "Start_Position"+ "\t" + "chrom"+ "\t" + "cDNA_Change"+ "\t" + "nm_Feature"
							  + "\t" + "nm_Featuretype" + "\t" + "EXON"+ "\t" + "gnomAD_AF" + "\t" +"1000genome_AF" + "\t" +"ESP_AA_AF" + "\t" +"ESP_EA_AF"
							  +"\t" + "dbSNP_ID" + "\t" + "COSMIC_ID" +"\n")
					continue
				if 'WildtypeProtein' not in readline:
					continue
				if 'intron_variant' in readline:
					continue
				if 'downstream_gene_variant' in readline:
					continue
				if 'regulatory_region_variant' in readline:
					continue
				if 'intergenic_variant' in readline:
					continue
				if 'splice_region_variant' in readline:
					continue
				if 'stop_gained' in readline:
					continue
				if 'stop_lost' in readline:
					continue
				if 'TF_binding_site_variant' in readline:
					continue
				if 'coding_sequence_variant' in readline:
					continue
				if 'upstream_gene_variant' in readline:
					continue
				if 'non_coding_transcript_exon_variant' in readline:
					continue
				if '5_prime_UTR_variant' in readline:
					continue
				if '3_prime_UTR_variant' in readline:
					continue
				else:
					mypep_seq = readline.split ("WildtypeProtein=", 1)[-1]
					mypep_seq = mypep_seq.replace ('\n', '')
					readline = readline.replace ('\n', '')
					ar = readline.strip ().split ('\t')

					aaChange_seq = str (ar[10]).split ("/")[1]
		#			aaRaw_seq = str (ar[10]).split ("/")[0]
					(seq_raw, loc) = get_25mers_MNP_rawseq (ar[9], mypep_seq)
					(seq_mut, len_mutation) = get_25mers_MNP_mut (ar[9], aaChange_seq, mypep_seq)
					if "SYMBOL" in readline:
						start = ";SYMBOL="
						end = ";SYMBOL_SOURCE="
						symbol = readline[readline.find (start) + len (start):readline.rfind (end)]
					else:
						symbol = "_"
					start_pos_tmp = str(ar[1]).split(":")[1]
					start_pos=start_pos_tmp.split("-")[0]
					chrom=str(ar[1]).split(":")[0]
					cDNA_change = "c." + ar[0].split ("_")[-1].split ("/")[0] + str (ar[8]) + ar[0].split ("_")[-1].split ("/")[1]
					nm_Feature= str(ar[4])
					nm_Featuretype = str (ar[5])
					exon=readline[readline.find (";EXON=") + len (";HGVSc"):readline.rfind (";HGVSc")].split("/")[0]
					if "gnomAD_AF" in readline:
						start1 = ";gnomAD_AF="
						end1 = ";gnomAD_AFR_AF="
						gnomAD_AF_tmp = readline[readline.find (start1) + len (start):readline.rfind (end1)]
						if ";" in gnomAD_AF_tmp:
							gnomAD_AF = gnomAD_AF_tmp.split (";")[0]
						else:
							gnomAD_AF = gnomAD_AF_tmp
					else:
						gnomAD_AF = "_"
					if ";AF" in readline:
						start1 = ";AF="
						end1 = ";AFR_AF="
						onegenome_AF_tmp = readline[readline.find (start1) + 1:readline.rfind (end1)]
						if "_" in onegenome_AF_tmp:
							onegenome_AF = "_"
						else:
							onegenome_AF = onegenome_AF_tmp
					else:
						onegenome_AF = "_"
					if ";AA_AF=" in readline:
						start1 = ";AA_AF="
						end1 = ";EA_AF="
						ESP_AA_AF = readline[readline.find (start1) + 1:readline.rfind (end1)]
					else:
						ESP_AA_AF = "_"
					if ";EA_AF=" in readline:
						start1 = ";EA_AF="
						end1 = ";gnomAD_AF="
						ESP_EA_AF_tmp = readline[readline.find (start1) + 1:readline.rfind (end1)]
						if ";" in ESP_EA_AF_tmp:
							ESP_EA_AF = ESP_EA_AF_tmp.split (";")[0]
						else:
							ESP_EA_AF = ESP_EA_AF_tmp
					else:
						ESP_EA_AF = "_"
					anno = get_dbsnp_cosmic (ar[12])
					fw.write ("".join (readline) + "\t" + str(seq_raw) + "\t" + str(seq_mut) + "\t" + str (symbol) +
								  "\t"+ "p." + str(ar[10]).split ("/")[0] + str(ar[9]) + str(ar[10]).split("/")[1] + "\t" + ar[0].split ("_")[-1].split ("/")[0] + "\t"
								  + str(loc)+ "\t" + str(len_mutation)+ "\t" + str(start_pos)+ "\t"+ str(chrom)+ "\t" + str(cDNA_change)+ "\t" + str(nm_Feature)+"\t" + str(nm_Featuretype)+
								  "\t" + 'exon'+str(exon)+ "\t" + str(gnomAD_AF)+"\t" + str(onegenome_AF)+"\t" + str(ESP_AA_AF)+"\t" + str(ESP_EA_AF)+"\t"+anno+"\n")

	fw.close ()


def get_25mers_MNP_rawseq(loc, seq_pep):
	if "-" not in loc:
		int_loc = int (loc)
	if "-" in loc:
		loc1 = loc.split ("-")[0]
		int_loc = int (loc1)
	if (len (seq_pep) - int_loc < 12):
		pos_pep = seq_pep[int_loc:len (seq_pep)]
	else:
		pos_pep = seq_pep[int_loc:int_loc + 12]
	if int_loc < 13:
		pre_pep = seq_pep[0:int_loc]
		newloc = int_loc
	else:
		pre_pep = seq_pep[int_loc - 13:int_loc]
		newloc = 13

	seq = pre_pep + pos_pep
	return (seq, newloc)


def get_25mers_MNP_mut(loc, aaChange_seq, seq_pep):
	len_aa=len(aaChange_seq)
	len_mut = len (aaChange_seq)
	if "-" not in loc:
		int_loc = int (loc)
	if "-" in loc:
		loc1 = loc.split ("-")[0]
		int_loc = int (loc1)
	if int_loc < 13:
		pre_pep = seq_pep[0:int_loc - 1]
		newloc = int_loc
	else:
		pre_pep = seq_pep[int_loc - 13:int_loc - 1]
		newloc = 13
	if (len (seq_pep) - int_loc < 12):
		pos_pep = seq_pep[(int_loc + len_aa - 1):len (seq_pep)]
	else:
		pos_pep = seq_pep[(int_loc + len_aa - 1):(int_loc + 12)]

	aa = Seq (aaChange_seq)
	seq = pre_pep + aa + pos_pep
	return (seq, len_mut)

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
	description='Obtain 25-mer AAchange up/downstream sequences')
	parser.add_argument ('-i', dest='input_dir', help='The inputdir path')
	parser.add_argument ('-a',dest = 'annovar_dir',help='Input annovar annotation dir', type=str)
	parser.add_argument ('-p', dest='pt_file', help='Tnput patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-t', '--type', dest='type', help="Type of somatic mutations", required=True,
						 choices=['all', 'SNP', 'MNP', 'INDEL'])
	parser.add_argument ('-e', dest='end', help='The end job id')
	parser.add_argument ('-o',dest = 'outdir', help='Output vep annotation file', type=str)

	args = parser.parse_args ()
	return args

def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	annovar_dir = args.annovar_dir
	outdir = args.outdir
	start = args.start
	end = args.end
	type = args.type
	pt_file = args.pt_file
	pts = open (pt_file)  # there is a header
	lns = pts.readlines ()
	if os.path.isdir(outdir):
		pass
	else :
		subprocess.call ('mkdir -p ' + outdir + '/', shell=True)
	for i in range ((int) (start), (int) (end) + 1):
		tmp = lns[i].strip ('\n').split (',')
		tumor_ID = tmp[1]  # tumor patient ID
		annovar_SNP = annovar_dir + "/" + tumor_ID + ".SNP.combine.txt"
		annovar_INDEL = annovar_dir + "/" + tumor_ID + ".INDEL.combine.txt"
		annovar_MNP = annovar_dir + "/" + tumor_ID + ".MNP.combine.txt"

		VEP_SNP = input_dir + "/" + tumor_ID + ".SNP.vepout.txt"
		VEP_INDEL = input_dir + "/" + tumor_ID + ".INDEL.vepout.txt"
		VEP_MNP = input_dir + "/" + tumor_ID + ".MNP.vepout.txt"

		VEP_SNP_c = input_dir + "/" + tumor_ID + ".SNP.vepout.canonical.txt"
		VEP_INDEL_c = input_dir + "/" + tumor_ID + ".INDEL.vepout.canonical.txt"
		VEP_MNP_c = input_dir + "/" + tumor_ID + ".MNP.vepout.canonical.txt"

		raw_SNP= outdir + "/" + tumor_ID + ".SNP.raw.txt"
		raw_INDEL = outdir + "/" + tumor_ID + ".INDEL.raw.txt"
		raw_MNP = outdir + "/" + tumor_ID + ".MNP.raw.txt"

		out_SNP= outdir + "/" + tumor_ID + ".SNP.25mer.txt"
		out_INDEL= outdir + "/" + tumor_ID + ".INDEL.25mer.txt"
		out_MNP= outdir + "/" + tumor_ID + ".MNP.25mer.txt"
		if type == 'all':
			if os.path.isfile (annovar_SNP) and os.path.isfile (VEP_SNP) and os.path.isfile (VEP_SNP_c):
				print ("SNP exist")
				get_pep_SNP (VEP_SNP, raw_SNP)
				get_map (raw_SNP, annovar_SNP, out_SNP, "SNP", outdir, tumor_ID, VEP_SNP_c)
			if os.path.isfile (annovar_MNP) and os.path.isfile (VEP_MNP) and os.path.isfile (VEP_MNP_c):
				print ("MNP exist")
				get_pep_MNP (VEP_MNP, raw_MNP)
				get_map_MNP (raw_MNP, annovar_MNP, out_MNP, "MNP", outdir, tumor_ID, VEP_MNP_c)
			if os.path.isfile (annovar_INDEL) and os.path.isfile (VEP_INDEL) and os.path.isfile (VEP_INDEL_c):
				print("INDEL exist")
				get_pep_INDEL (VEP_INDEL, raw_INDEL)
				get_map (raw_INDEL, annovar_INDEL, out_INDEL, "INDEL", outdir, tumor_ID, VEP_INDEL_c)
		if type == "SNP":
			get_pep_SNP (VEP_SNP, raw_SNP)
			get_map (raw_SNP, annovar_SNP, out_SNP, "SNP", outdir,tumor_ID,VEP_SNP_c)
		if type == "MNP":
			get_pep_MNP (VEP_MNP, raw_MNP)
			get_map_MNP (raw_MNP, annovar_MNP, out_MNP, "MNP", outdir,tumor_ID,VEP_MNP_c)
		if type == "INDEL":
			get_pep_INDEL (VEP_INDEL, raw_INDEL)
			get_map (raw_INDEL, annovar_INDEL, out_INDEL, "INDEL", outdir, tumor_ID,VEP_INDEL_c)


if __name__ == '__main__':
	main ()

'''
python /Users/cyu/Documents/work/mc3data/2020JUL14/vep/get_25mer.py \
-i /Users/cyu/Documents/work/mislabel_ICON/VEP/output/ \
-a /Users/cyu/Documents/work/mislabel_ICON/mafcombine/ \
-p /Users/cyu/Documents/work/mislabel_ICON/sample_mislabel_pairs1.txt \
-s 2 -e 2 -t SNP \
-o /Users/cyu/Documents/work/mislabel_ICON/VEP/25mer/


python /Users/cyu/Documents/work/mc3data/2020JUL14/vep/get_25mer.py \
-i /Users/cyu/Documents/work/mislabel_ICON/VEP/output/ \
-a /Users/cyu/Documents/work/mislabel_ICON/mafcombine/ \
-p /Users/cyu/Documents/work/mislabel_ICON/sample_mislabel_pairs1.txt \
-s 1 -e 2 -t INDEL \
-o /Users/cyu/Documents/work/mislabel_ICON/VEP/25mer/

python /Users/cyu/Documents/work/mc3data/2020JUL14/vep/get_25mer.py \
-i /Users/cyu/Documents/work/mislabel_ICON/VEP/output/ \
-a /Users/cyu/Documents/work/mislabel_ICON/mafcombine/ \
-p /Users/cyu/Documents/work/mislabel_ICON/sample_mislabel_pairs1.txt \
-s 2 -e 2 -t MNP \
-o /Users/cyu/Documents/work/mislabel_ICON/VEP/25mer/



/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/2_mutationcalling/get_25mer.py \
-i /Users/cyu/Documents/work/gbm/VEP/output/SNP/ \
-a /Users/cyu/Documents/work/mislabel_ICON/mafcombine/ \
-p /Users/cyu/Documents/work/mislabel_ICON/sample_mislabel_pairs1.txt \
-s 2 -e 2 -t MNP \
-o /Users/cyu/Documents/work/mislabel_ICON/VEP/25mer/



/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/2_mutationcalling/get_25mer.py \
-i /Users/cyu/Documents/work/gbm/VEP/output/ \
-a /Users/cyu/Documents/work/gbm/annovar/mafcombine/ \
-p /Users/cyu/Documents/work/gbm/fileinfo/sample_pairs.txt \
-s 1 -e 3 -t SNP \
-o /Users/cyu/Documents/work/gbm/VEP/25mer/


/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/2_mutationcalling/get_25mer.py \
-i /Users/cyu/Documents/test/VEP/output/ \
-a /Users/cyu/Documents/work/gbm/annovar/mafcombine/ \
-p /Users/cyu/Documents/work/gbm/fileinfo/sample_pairs.txt \
-s 2 -e 2 -t MNP \
-o /Users/cyu/Documents/test/VEP/25mer/




/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/2_mutationcalling/get_25mer.py \
-i /data/cyu/gbm/data/VEP/output/ \
-a /data/cyu/gbm/data/annovar/mafcombine/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 1 -e 5 -t SNP \
-o /data/cyu/gbm/data/VEP/25mer/

/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/2_mutationcalling/get_25mer.py \
-i /data/cyu/gbm/data/VEP/output/ \
-a /data/cyu/gbm/data/annovar/mafcombine/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 1 -e 5 -t INDEL \
-o /data/cyu/gbm/data/VEP/25mer/


/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/2_mutationcalling/get_25mer.py \
-i /data/cyu/gbm/data/VEP/output/ \
-a /data/cyu/gbm/data/annovar/mafcombine/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 2 -e 2 -t MNP \
-o /data/cyu/gbm/data/VEP/25mer/


/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/2_mutationcalling/get_25mer.py \
-i /data/cyu/gbm/data/VEP/output/ \
-a /data/cyu/gbm/data/annovar/mafcombine/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 2 -e 2 -t all \
-o /data/cyu/gbm/data/VEP/25mer/


'''
