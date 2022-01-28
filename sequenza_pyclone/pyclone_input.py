import os, sys, getopt
from collections import defaultdict, namedtuple
import pandas as pd



def usage():
	usage = """
		This script takes as an input de _segments.txt file resulting from Sequenza and a vcf file containing variants.
		The scripts outputs a tsv file that serves as an input for PyClone with --prior parameters set to major_copy_number.
		These are the columns in the tsv output:
			- mutation_id - A unique ID to identify the mutation.
			- ref_counts - The number of reads covering the mutation which contain the reference (genome) allele.
			- var_counts - The number of reads covering the mutation which contain the variant allele.
			- normal_cn - The copy number of the cells in the normal population. For autosomal chromosomes this will be 2 and for sex chromosomes it could be either 1 or 2. For species besides human other values are possible.
			- minor_cn - The minor copy number of the cancer cells. 
			- major_cn - The major copy number of the cancer cells.
		The warning_sequenza_mutations.tsv file contains variants that are not present in the Sequenza segments.
		Required arguments:                Description                                             
		-i  --sequenza_input               _segments.txt file containing segements and copy numbers obtained by Sequenza
		-v  --vcf_input                    VCF file (MuTect2) containing variants that will serve as input for PyClone
		Optional arguments:
		-o  --output_prefix                 Prefix of the output tsv file
		-h, --help                          Print this help information and exit
		"""

	print (usage)


# Function - read in options
def read_options(argv):
	try:
		optlist, args = getopt.getopt (argv,
									   'i:v:r:o:h',
									   ['sequenza_input', 'vcf_input', 'raw_output','output', 'help'])
		if not optlist:
			print ('No options supplied')
			usage ()
	except getopt.GetoptError:
		usage ();
		sys.exit ("Input errors")
	# Create dictionary of long and short formats
	format_dict = {
		'-i': '--sequenza_input',
		'-v': '--vcf_input',
		'-r': '--raw_output',
		'-o': '--output',
		'-h': '--help'
	}

	# Create a dictionary of options and input from the options list
	opts = dict (optlist)

	# Use the long format dictionary to change the option to the short annotation, if long is given by the user.
	for short, long_ in format_dict.items ():
		if long_ in opts:
			opts[short] = opts.pop (long_)

	# Print usage help
	if '-h' in opts.keys ():
		usage ();
		sys.exit ()

	# Define values
	sequenza_input = opts['-i'] if '-i' in opts.keys () else None
	if sequenza_input == None:
		usage ();
		sys.exit ('Sequenza segement file missing')
	vcf_input = opts['-v'] if '-v' in opts.keys () else None
	if vcf_input == None:
		usage ();
		sys.exit ('VCF file missing')
	output = opts['-o'] if '-o' in opts.keys () else None
	if output == None:
		output = 'PyClone_input.tsv'
	raw_output = opts['-r'] if '-r' in opts.keys () else None
	if raw_output == None:
		output = 'PyClone_raw_input.tsv'
	# Create and fill input named-tuple
	Input = namedtuple ('input', ['sequenza_input', 'vcf_input', 'raw_output','output'])
	inputinfo = Input (sequenza_input, vcf_input, raw_output, output)

	return inputinfo

def getCN(cnt,lcn):
	"""
	Depending on the value of cnt and lcn, returns if there is a copy number aberration:(A)mplification, (D)eletion, Copy number (N)ormal, or Copy Number Neutral (L)oss of Heterozygosity
	Parameters :
		cnt (int) : Total copy number of the region
		lcn (int) : Low copy number of the region
	Returns :
		char : Amplification/deletion/normal/Copy Number Neutral
	"""
	ret = ""
	if cnt > 2:
		ret = "Amplification"
	elif cnt == 2 and lcn == 1:
		ret = "Normal"
	elif cnt == 2 and lcn == 0:
		ret = "Copy Number Neutral"
	elif cnt == 2 and lcn == 2:
		ret = "Copy Number Neutral"
	elif cnt < 2:
		ret = "Deletion"
	else:
		ret = "_"
	return ret
def main(args):
	# Read input
	input_ = read_options (args)

	sequenza = pd.read_csv (input_.sequenza_input, sep='\t')
	#print(sequenza['chromosome'])
	sequenza['chromosome']=sequenza['chromosome'].astype(str)
	sequenza['chromosome'] = sequenza['chromosome'].apply (lambda x: x[0:])


	inputfile = open (input_.vcf_input, 'r')
	outputfile = input_.output
	outputfile1 = input_.raw_output
	mutation_id_list = []
	ref_counts_list = []
	var_counts_list = []
	normal_cn_list = []
	minor_cn_list = []
	major_cn_list = []
	warning_list = []
	ref_g_list = []
	alt_g_list = []
	gene_list = []
	sample_list = []
	copy_list=[]

	for line in inputfile:

		if line.startswith('Chr'):
			continue
		else:
		# from vcf extract mutation_id, var_counts and ref_counts
			record=line.strip().split("\t")
#            print(record)
			chrom = record[0]
			pos = record[1]
			mutation_id = str (chrom) + ':' + str (pos)
			ref_counts = record[9]
			var_counts = record[10]
			t_vaf = record[7]
#            n_vaf = record[8]
#            n_ref_counts = record[11]
#            n_alt_counts = record[12]
			ref_g = record[3]
			alt_g = record[4]
			gene=record[5]
			sample=record[21]


		# extract normal_cn, minor_cn and major_cn from accucopy cnv.output.tsv

			total_cn_df = sequenza.loc[(sequenza['chromosome'] == str (chrom)) & (sequenza['start.pos'] <= int(pos))
									   & (sequenza['end.pos'] >= int(pos)), 'CNt']
#			print(total_cn_df)
			if sum (sequenza.chromosome.str.contains ("Y")) > 0:
				if chrom == 'Y' or chrom == 'X':
					normal_cn = 1
				else:
					normal_cn = 2
			else:
				normal_cn = 2

			major_cn_df = sequenza.loc[(sequenza['chromosome'] == str (chrom)) & (sequenza['start.pos'] <= int(pos))
								   & (sequenza['end.pos'] >= int(pos)), 'A']

			minor_cn_df = sequenza.loc[(sequenza['chromosome'] == str (chrom)) & (sequenza['start.pos'] <= int (pos))
								   & (sequenza['end.pos'] >= int (pos)), 'B']


#			print (mutation_id,total_cn_df.values[0], minor_cn_df.values[0])

			if not total_cn_df.empty:
				if not major_cn_df.empty:
					normal_cn_list.append (normal_cn)
					major_cn_list.append (major_cn_df.values[0])
					minor_cn_list.append (total_cn_df.values[0] - major_cn_df.values[0])
					mutation_id_list.append (mutation_id)
					ref_counts_list.append (ref_counts)
					var_counts_list.append (var_counts)
					ref_g_list.append (ref_g)
					alt_g_list.append(alt_g)
					gene_list.append(gene)
					sample_list.append(sample)
					ret = getCN (total_cn_df.values[0], minor_cn_df.values[0])
					copy_list.append(ret)

				else:
					minor_cn = 0
					normal_cn_list.append (normal_cn)
					major_cn_list.append (total_cn_df.values[0])
					minor_cn_list.append (minor_cn)
					mutation_id_list.append (mutation_id)
					ref_counts_list.append (ref_counts)
					var_counts_list.append (var_counts)
					ref_g_list.append (ref_g)
					alt_g_list.append(alt_g)
					gene_list.append(gene)
					sample_list.append(sample)
					ret = getCN (total_cn_df.values[0], minor_cn_df.values[0])
					copy_list.append (ret)
			else:  # if mutation lies not in a segment, the major_cn is set to the normal_cn, this mutation is flagged in the warning_sequenza_mutations.tsv file
				minor_cn = 0
				normal_cn_list.append (normal_cn)
				major_cn_list.append (1)
				minor_cn_list.append (1)
				mutation_id_list.append (mutation_id)
				ref_counts_list.append (ref_counts)
				var_counts_list.append (var_counts)
				ref_g_list.append (ref_g)
				alt_g_list.append (alt_g)
				gene_list.append (gene)
				sample_list.append (sample)
				ret = getCN (2, 1)
				copy_list.append (ret)
				warning_list.append (mutation_id)


	df = pd.DataFrame ({'mutation_id': mutation_id_list, 'ref_counts': ref_counts_list,
						'var_counts': var_counts_list, 'normal_cn': normal_cn_list, 'major_cn': major_cn_list,
						'minor_cn': minor_cn_list, 'Gene_Name': gene_list, 'Sample_ID':sample_list, 'Copy_number_Abberation':copy_list},
					   columns=['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'major_cn', 'minor_cn','Gene_Name','Sample_ID','Copy_number_Abberation'])

	df = df.dropna ()
	df['normal_cn'] = df['normal_cn'].astype (int)
	df['major_cn'] = df['major_cn'].astype (int)
	df['minor_cn'] = df['minor_cn'].astype (int)
	Genotype = []
	for i in range (len (df.major_cn)):
		if df.major_cn[i] > 0 and df.minor_cn[i] > 0:
			gt = 'A' * df.minor_cn[i] + 'B' * df.major_cn[i]
		elif df.major_cn[i] > 0 and df.minor_cn[i] == 0:
			gt = "B" * df.major_cn[i] + ''
		elif df.major_cn[i] == 0 and df.minor_cn[i] > 0:
			gt = "A" * df.minor_cn[i] + ''
		elif df.major_cn[i] == 0 and df.minor_cn[i] == 0:
			gt = '__'
		else:
			print ("Ops,this cause by chrosome Y total copynumber undetecting!")
		Genotype.append (gt)
	df['genotype'] = Genotype
##change column order
	df=df[['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'major_cn', 'minor_cn','genotype','Gene_Name','Sample_ID','Copy_number_Abberation']]


	df_warning = pd.DataFrame ({'mutation_id': warning_list})
	df_filter=df.loc[df['major_cn']>0]
	df.to_csv (outputfile1, sep='\t', index=False)
	df_filter.to_csv (outputfile, sep='\t', index=False)
	df_warning.to_csv ('warning_sequenza_mutations.tsv', sep='\t', index=False)


######################################################################
######                         RUNNER                           ######
######################################################################
if __name__ == '__main__':
	main (sys.argv[1:])
