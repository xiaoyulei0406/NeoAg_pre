import sys
import matplotlib
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
def get_df(in_file,col_name): ##read annovar file from Mutect,Mutect2, Muse, Strelka, Varscan2
	tumor_vaf_col = 'tumor_vaf_' + col_name
	tumor_ref_count_col = 'tumor_ref_count_' + col_name
	tumor_alt_count_col = 'tumor_alt_count_' + col_name
	normal_vaf_col = 'normal_vaf_' + col_name
	normal_ref_count_col = 'normal_ref_count_' + col_name
	normal_alt_count_col = 'normal_alt_count_' + col_name
	if os.path.isfile(in_file):
		df=pd.read_csv (in_file, sep='\t')
		df_mut=df[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
					   'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand',
					   'tumor_vaf','tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]
		df_mut.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand',tumor_vaf_col, tumor_ref_count_col, tumor_alt_count_col, normal_vaf_col,
						 normal_ref_count_col, normal_alt_count_col,
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
		df_mut[tumor_vaf_col] = df_mut[tumor_vaf_col].apply (lambda x: round (x, 3))
		df_mut[normal_vaf_col] = df_mut[normal_vaf_col].apply (lambda x: round (x, 3))
		df_mut['Chromosome']=df_mut['Chromosome'].astype(str)#change data types from int/object to object
		return(df_mut)
	else:
		df_mut = pd.DataFrame(columns=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
					 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
					 'cytoBand',tumor_vaf_col, tumor_ref_count_col, tumor_alt_count_col, normal_vaf_col,
						 normal_ref_count_col, normal_alt_count_col,
					 'Tumor_Sample_Barcode', 'snp142', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN',
					 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'avsnp147'])
		return (df_mut)

def get_combineMAF_SNP(ptID,input_dir,outDir):
	in_VS = input_dir + ptID + '.Varscan2.SNPs.txt'
	in_M1 = input_dir + ptID + '.Mutect.SNPs.txt'
	in_M2 = input_dir + ptID + '.Mutect2.SNPs.txt'
	in_MS = input_dir + ptID + '.Muse.SNPs.txt'
	in_SK = input_dir + ptID + '.Strelka.SNPs.txt'

	############################################################

	##column name
	sample_VS = "VS"
	sample_M2 = "M2"
	sample_M1 = "M1"
	sample_MS = "MS"
	sample_SK = "SK"

	df_mut_VS=get_df(in_VS,sample_VS)##Create new dataframe with selected columns by Varscan2
	df_mut_M2=get_df(in_M2,sample_M2)##Create new dataframe with selected columns by Mutect2
	df_mut_M1=get_df(in_M1,sample_M1)##Create new dataframe with selected columns by Mutect
	df_mut_MS=get_df(in_MS,sample_MS)##Create new dataframe with selected columns by MuSE
	df_mut_SK=get_df(in_SK,sample_SK)##Create new dataframe with selected columns by Strelka

	##Start merge all annovar results from 5 tools
	df_join = pd.merge (df_mut_VS, df_mut_M2,
						on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'],
						how='outer', indicator=False)
	df_join= pd.merge (df_join, df_mut_M1,
						on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'],
						how='outer', indicator=False)
	df_join= pd.merge (df_join, df_mut_MS,
						on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'],
						how='outer', indicator=False)

	df_join= pd.merge (df_join, df_mut_SK,
						on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'],
						how='outer', indicator=False)
	cols = []
	for indx, row in df_join.iterrows ():
		temp = []
		if (df_mut_VS[((df_mut_VS['Chromosome'] == row['Chromosome']) & (df_mut_VS['Start_Position'] == row['Start_Position']) &
					   (df_mut_VS['End_Position'] == row['End_Position']) & (df_mut_VS['Reference_Allele'] == row['Reference_Allele']) &
					   (df_mut_VS['Tumor_Seq_Allele2'] == row['Tumor_Seq_Allele2']))].shape[0] == 1):
			temp.append ("Varscan2")
		else:
			temp.append (np.nan)
		if (df_mut_M2[((df_mut_M2['Chromosome'] == row['Chromosome']) & (df_mut_M2['Start_Position'] == row['Start_Position']) &
					   (df_mut_M2['End_Position'] == row['End_Position']) & (df_mut_M2['Reference_Allele'] == row['Reference_Allele']) &
					   (df_mut_M2['Tumor_Seq_Allele2'] == row['Tumor_Seq_Allele2']))].shape[0] == 1):
			temp.append ("Mutect2")
		else:
			temp.append (np.nan)
		if (df_mut_M1[((df_mut_M1['Chromosome'] == row['Chromosome']) & (df_mut_M1['Start_Position'] == row['Start_Position']) &
					   (df_mut_M1['End_Position'] == row['End_Position']) & (df_mut_M1['Reference_Allele'] == row['Reference_Allele']) &
					   (df_mut_M1['Tumor_Seq_Allele2'] == row['Tumor_Seq_Allele2']))].shape[0] == 1):
			temp.append ("MuTect")
		else:
			temp.append (np.nan)
		if (df_mut_MS[((df_mut_MS['Chromosome'] == row['Chromosome']) & (df_mut_MS['Start_Position'] == row['Start_Position']) &
					   (df_mut_MS['End_Position'] == row['End_Position']) & (df_mut_MS['Reference_Allele'] == row['Reference_Allele']) &
					   (df_mut_MS['Tumor_Seq_Allele2'] == row['Tumor_Seq_Allele2']))].shape[0] == 1):
			temp.append ("MuSE")
		else:
			temp.append (np.nan)
		if (df_mut_SK[((df_mut_SK['Chromosome'] == row['Chromosome']) & (df_mut_SK['Start_Position'] == row['Start_Position']) &
					   (df_mut_SK['End_Position'] == row['End_Position']) & (df_mut_SK['Reference_Allele'] == row['Reference_Allele']) &
					   (df_mut_SK['Tumor_Seq_Allele2'] == row['Tumor_Seq_Allele2']))].shape[0] == 1):
			temp.append ("Strelka")
		else:
			temp.append (np.nan)

		cols.append (temp)
	df_join[sample_VS] = [x[0] for x in cols]
	df_join[sample_M2] = [x[1] for x in cols]
	df_join[sample_M1] = [x[2] for x in cols]
	df_join[sample_MS] = [x[3] for x in cols]
	df_join[sample_SK] = [x[4] for x in cols]

	### add all callers into a new column
	df_join['Callers']=df_join[[ "VS",'M2', 'M1', "MS", "SK"]].apply(lambda x: ';'.join(x.dropna()), axis=1)
	### add number of callers into a new column
	df_join['NCallers'] = df_join['Callers'].str.count (';') + 1
	df_join = df_join.fillna (0)
	# Get a max containing maximum value of each row
	maxtumor_f = df_join[["tumor_vaf_VS", "tumor_vaf_M2", "tumor_vaf_M1", "tumor_vaf_MS", "tumor_vaf_SK"]].max (axis=1)
	df_join['tumor_vaf'] = maxtumor_f


	df_join=df_join[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2','tumor_vaf','Callers','NCallers']]

	df_mut_VS.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
						 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
	df_join_VS = pd.merge (df_join, df_mut_VS,on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2','tumor_vaf'], how='inner', indicator=False)

	# select out new columns
	df_join_VS = df_join_VS[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							 'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
							 'cytoBand','Callers','NCallers','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf',
							 'normal_ref_count', 'normal_alt_count',
							 'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]

	df_mut_M2.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
						 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
	df_join_M2 = pd.merge (df_join, df_mut_M2,
						  on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							  'tumor_vaf'], how='inner', indicator=False)
	# select out new columns
	df_join_M2 = df_join_M2[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							 'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
							 'cytoBand','Callers','NCallers','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf',
							 'normal_ref_count', 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]

	df_mut_M1.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
						 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
	df_join_M1 = pd.merge (df_join, df_mut_M1,
						  on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							  'tumor_vaf'], how='inner', indicator=False)
	# select out new columns
	df_join_M1 = df_join_M1[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							 'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand',
							 'Callers','NCallers','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
							 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]

	df_mut_MS.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
						 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
	df_join_MS = pd.merge (df_join, df_mut_MS,
						  on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							  'tumor_vaf'], how='inner', indicator=False)
	# select out new columns
	df_join_MS = df_join_MS[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							 'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand',
							 'Callers','NCallers','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]

	df_mut_SK.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
						 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
	df_join_SK = pd.merge (df_join, df_mut_SK,
						  on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							  'tumor_vaf'], how='inner', indicator=False)
	df_join_SK = df_join_SK[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							 'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand',
							 'Callers','NCallers','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]
	#print ("join_SK.shape:")
	#print(df_join_SK.shape)

	# merge dataframe
	df_join = pd.DataFrame (np.concatenate ([df_join_VS.values, df_join_M2.values, df_join_M1.values, df_join_MS.values, df_join_SK.values]),columns=df_join_VS.columns)

	df_join_uniq = df_join.sort_values ('tumor_vaf').drop_duplicates ('Start_Position', keep='last')
	df_filter=df_join_uniq.loc[df_join_uniq['NCallers'] > 2] #select out mutations with 3 callers
	df_filter.to_csv (outDir + "/" + ptID +".SNP.combine.txt", index=False, header=True, sep='\t')
	print ("df_filter.shape:")
	print (df_filter.shape)
	return (df_filter)


def get_combineMAF_INDEL(ptID,input_dir,outDir):
	in_VS = input_dir + ptID + '.Varscan2.INDELs.txt'
	in_M2 = input_dir + ptID + '.Mutect2.INDELs.txt'
	in_SK = input_dir + ptID + '.Strelka.INDELs.txt'

	############################################################

	##column name
	sample_VS = "VS"
	sample_M2 = "M2"
	sample_SK = "SK"
	df_mut_VS = get_df (in_VS, sample_VS)  ##Create new dataframe with selected columns by Varscan2
	df_mut_M2 = get_df (in_M2, sample_M2)  ##Create new dataframe with selected columns by Mutect2
	df_mut_SK = get_df (in_SK, sample_SK)  ##Create new dataframe with selected columns by Strelka

	df_join = pd.merge (df_mut_VS, df_mut_M2,
						on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'],
						how='outer', indicator=False)
	df_join = df_join.fillna (0)
	df_join['Chromosome'] = df_join['Chromosome'].astype(str)

	df_join= pd.merge (df_join, df_mut_SK,
						on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'],
						how='outer', indicator=False)

	cols = [] ##create new dataframe for methods summary
	for indx, row in df_join.iterrows ():
		temp = []
		if (df_mut_VS[((df_mut_VS['Chromosome'] == row['Chromosome']) & (df_mut_VS['Start_Position'] == row['Start_Position']) &
					   (df_mut_VS['End_Position'] == row['End_Position']) & (df_mut_VS['Reference_Allele'] == row['Reference_Allele']) &
					   (df_mut_VS['Tumor_Seq_Allele2'] == row['Tumor_Seq_Allele2']))].shape[0] == 1):
			temp.append ("Varscan2")
		else:
			temp.append (np.nan)
		if (df_mut_M2[((df_mut_M2['Chromosome'] == row['Chromosome']) & (df_mut_M2['Start_Position'] == row['Start_Position']) &
					   (df_mut_M2['End_Position'] == row['End_Position']) & (df_mut_M2['Reference_Allele'] == row['Reference_Allele']) &
					   (df_mut_M2['Tumor_Seq_Allele2'] == row['Tumor_Seq_Allele2']))].shape[0] == 1):
			temp.append ("Mutect2")
		else:
			temp.append (np.nan)
		if (df_mut_SK[((df_mut_SK['Chromosome'] == row['Chromosome']) & (df_mut_SK['Start_Position'] == row['Start_Position']) &
					   (df_mut_SK['End_Position'] == row['End_Position']) & (df_mut_SK['Reference_Allele'] == row['Reference_Allele']) &
					   (df_mut_SK['Tumor_Seq_Allele2'] == row['Tumor_Seq_Allele2']))].shape[0] == 1):
			temp.append ("Strelka")
		else:
			temp.append (np.nan)

		cols.append (temp)
	df_join[sample_VS] = [x[0] for x in cols]
	df_join[sample_M2] = [x[1] for x in cols]
	df_join[sample_SK] = [x[2] for x in cols]

	### add all callers into a new column
	df_join['Callers']=df_join[['VS', 'M2',"SK",]].apply(lambda x: ';'.join(x.dropna()), axis=1)
	### add number of callers into a new column
	df_join['NCallers'] = df_join['Callers'].str.count (';') + 1
	# print(df_join)
	df_join = df_join.fillna (0)
	# Get a series containing maximum value of each row
	maxtumor_f = df_join[["tumor_vaf_VS", "tumor_vaf_M2",  "tumor_vaf_SK"]].max (axis=1)
	df_join['tumor_vaf'] = maxtumor_f


	df_join=df_join[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2','tumor_vaf','Callers','NCallers']]

	df_mut_VS.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
						 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
	df_join_VS = pd.merge (df_join, df_mut_VS,
						  on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							  'tumor_vaf'], how='inner', indicator=False)
	# select out new columns
	df_join_VS = df_join_VS[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							 'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand',
							 'Callers','NCallers','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]
	df_mut_M2.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
						 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
	df_join_M2 = pd.merge (df_join, df_mut_M2,
						  on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							  'tumor_vaf'], how='inner', indicator=False)
	# select out new columns
	df_join_M2 = df_join_M2[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							 'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand',
							 'Callers','NCallers','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]
	df_mut_SK.columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
						 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
						 'cytoBand','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count',
						 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']
	df_join_SK = pd.merge (df_join, df_mut_SK,
						  on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							  'tumor_vaf'], how='inner', indicator=False)
	# select out new columns
	df_join_SK = df_join_SK[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
							 'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand',
							 'Callers','NCallers','tumor_vaf', 'tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
					   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]

	# merge dataframe
	df_join = pd.DataFrame (np.concatenate ([df_join_VS.values, df_join_M2.values, df_join_SK.values]),columns=df_join_VS.columns)
	df_join_uniq = df_join.sort_values ('normal_vaf').drop_duplicates ('Start_Position', keep='last')
	df_join_uniq = df_join_uniq.sort_values ('tumor_vaf').drop_duplicates ('Start_Position', keep='last')
	df_filter = df_join_uniq.loc[df_join_uniq['NCallers'] > 2] #select out mutations with 3 callers
	df_filter.to_csv (outDir + "/" + ptID +".INDEL.combine.txt", index=False, header=True, sep='\t')
	print ("df_filter.shape:")
	print (df_filter.shape)
	return(df_filter)


def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Combine all the annotation results')
	parser.add_argument ('-i', dest='input_dir', help='The inputdir path')
	parser.add_argument ('-p', dest='pt_file', help='Input patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	args = parser.parse_args ()
	return args


def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	out_dir = args.out_dir
	pt_file = args.pt_file
	start = args.start
	end = args.end

	pts = open (pt_file)  # paired sample information
	lns = pts.readlines ()

	for i in range ((int) (start), (int) (end) + 1):
		print (lns[i])
		tmp = lns[i].strip('\n').split(',')
		ptID = tmp[1] #tumor patient ID
		subprocess.call ('mkdir -p ' + './' + out_dir + '/', shell=True)
		in_MNPs = input_dir + ptID + '.Mutect2.MNPs.txt'
		df_MNP_M2 = pd.read_csv (in_MNPs, sep='\t')
		df_MNP_M2 = df_MNP_M2[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand',
							   'Callers','NCallers','tumor_vaf','tumor_ref_count', 'tumor_alt_count', 'normal_vaf', 'normal_ref_count', 'normal_alt_count',
							   'Tumor_Sample_Barcode','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147']]
		df_SNP = get_combineMAF_SNP (ptID,input_dir,out_dir)
		df_INDEL = get_combineMAF_INDEL (ptID, input_dir, out_dir)
		df_join = pd.DataFrame (np.concatenate ([df_SNP.values, df_INDEL.values,df_MNP_M2.values]),
								columns=df_INDEL.columns)
		df_join.to_csv (out_dir + "/" + ptID + ".combine.all.txt", index=False, header=True, sep='\t')
		df_MNP_M2.to_csv (out_dir + "/" + ptID + ".MNP.combine.txt", index=False, header=True, sep='\t')

	pts.close ()


if __name__ == '__main__':
	main()

'''

python /Users/cyu/Documents/work/mc3data/2020JUL14/somaticCalling/combine_annovar_alltools.py \
-i /Users/cyu/Documents/work/gbm/annovar/single_mafs/ \
-p /Users/cyu/Documents/work/gbm/fileinfo/sample_pairs.txt \
-s 1 -e 5 -o /Users/cyu/Documents/work/gbm/annovar/mafcombine/

/home/ec2-user/anaconda3/bin/python /data/cyu/scripts/wes/combine_annovar_alltools.py \
-i /data/cyu/gbm/data/annovar/single_mafs/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 1 -e 5 -o /data/cyu/gbm/data/annovar/mafcombine1/


sample_pairs.txt:

Normal,Tumor
457001_13S01187880_S1,457001_13S00073824_S7
457003_13S33628480_S3,457003_13S01573060_S8
457004_13S34677088_S4,457004_13S00720972_S9
457005_13S01056776_S5,457005_S1300524484_S10
457006_13S01712288_S6,457006_13S01638476_S11


'''
