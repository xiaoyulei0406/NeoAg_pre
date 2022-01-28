import argparse
import os
import pandas as pd
import glob
import openpyxl
import xlwt
def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
		description='Summarize all QC table from different tools')
	parser.add_argument ('-i', dest='indir', help='Input directionary of flagstat, qualimap and merits results', type=str)
	parser.add_argument ('-o', dest='output', help='Output for QC table', type=str)
	args = parser.parse_args ()
	return args


def calc_coverage(bed_file):
	if not os.path.exists (bed_file):
		print ("Not exist file, " + bed_file)
		return
	bed = pd.read_csv (bed_file, header=None, sep='\t', usecols=[1, 2, 3, 5, 6])
	# 100X coverage rates
	X100_targets = bed[bed[3] >= 100]
	X100_targets_size = X100_targets.shape[0]
	bases = bed.shape[0]
	target_100X = round (float (X100_targets_size) / bases * 100, 2)
	# 150X coverage rates
	X150_targets = bed[bed[3] >= 150]
	X150_targets_size = X150_targets.shape[0]
	target_150X = round (float (X150_targets_size) / bases * 100, 2)
	# 150X coverage rates
	X200_targets = bed[bed[3] >= 200]
	X200_targets_size = X200_targets.shape[0]
	target_200X = round (float (X200_targets_size) / bases * 100, 2)
	# 5000X coverage rates
	X500_targets = bed[bed[3] >= 500]
	X500_targets_size = X500_targets.shape[0]
	target_500X = round (float (X500_targets_size) / bases * 100, 2)
	return (target_100X, target_150X, target_200X, target_500X)


def calc_maprates(srtbamflagstat, rmdbamflagstat, rmdmerits):
	with open (rmdmerits, 'r') as f:
		for i in range (1, 8):
			tmp = f.readline ()
		for i in range (1, 2):
			line = f.readline ()
			line_tmp = line.strip ().split ("\t")
			dupRates = round (float (line_tmp[8]) * 100, 2)
	f.close
	with open (srtbamflagstat, 'r') as f1:
		for i in range (1, 2):
			line = f1.readline ()
			mapReads_tmp = line.split (" ")[0]
			mapReads=round(float(mapReads_tmp)/1000000,2)
	f1.close
	with open (rmdbamflagstat, 'r') as f2:
		for i in range (1, 5):
			tmp = f2.readline ()
		for i in range (1, 2):
			line = f2.readline ()
			mapRates_tmp = line.split (" ")[4]
			mapRates_tmp1 = str (mapRates_tmp).split ("(")[1]
			mapRates = str (mapRates_tmp1).split ("%")[0]
			continue
	f2.close
	return (mapReads, mapRates, dupRates)

def calc_qualimap(genome_results):
	with open (genome_results, 'r') as f3:
		for line in f3.readlines():
			if ("mean coverageData" in line):
				meancoverage_tmp = line.split (" = ")[1]
				meancoverage_tmp1 = meancoverage_tmp.split ("X")[0]
				meancoverage=round(float(meancoverage_tmp1),2)
	f3.close
	return meancoverage

def main():
	args = parse_arguments ()
	output = args.output
	inputDir = args.indir
	covArr=[]
	allSBstats=[]
	allRBstats=[]
	allMetrics=[]
	allGres=[]
	allcovs=[]
	for root, dirs, files in os.walk(inputDir):
		for filename in files:
			if filename.endswith (".coverage.csv"):
				coverage = os.path.join (root, filename)
				allcovs.append(coverage)
			if filename.endswith ("_genome_results.txt"):
				quali = os.path.join (root, filename)
				allGres.append(quali)
			if filename.endswith (".rmdup.bam.flagstat"):
				rBstat = os.path.join (root, filename)
				allRBstats.append(rBstat)
			if filename.endswith (".rmdup.metrics.txt"):
				merits = os.path.join (root, filename)
				allMetrics.append(merits)
			if filename.endswith(".sorted.bam.flagstat"):
				sBstat=os.path.join(root,filename)
				allSBstats.append(sBstat)



	allcovs1=allcovs.sort()
	allGres1=allGres.sort()
	allRBstats1=allRBstats.sort()
	allMetrics1=allMetrics.sort()
	allSBstats1=allSBstats.sort()

	for i in range(0,len(allGres)):
		print(allGres[i])
		ptID_tmp=allSBstats[i].split("/")[-1]
		print(ptID_tmp)
		ptID = str(ptID_tmp).split (".")[0]
		(target_100X, target_150X, target_200X, target_500X) = calc_coverage(allcovs[i])
		ave_target_base = calc_qualimap (allGres[i])
		(total_reads, mapping_rates, dup_rates) = calc_maprates (allSBstats[i], allRBstats[i], allMetrics[i])
		data_string1 = "{0},{1},{2},{3},{4},{5},{6},{7},{8}".format (ptID,total_reads, mapping_rates, dup_rates,
														 ave_target_base, target_100X, target_150X,
														 target_200X, target_500X)
		dataList = list (data_string1.split (","))
		covArr.append (dataList)
	y=pd.DataFrame(covArr,columns = ["sampleID","total_reads(million)","mapping_rates%","dup_rates%","average_cov",
									 "100X_targets%","150X_targets%","200X_targets%","500X_targets%"])
	#print(y)
	y.to_excel (output, index=False,float_format="%.2f",sheet_name='Sheet1', engine='xlsxwriter')


if __name__ == '__main__':
	main ()
