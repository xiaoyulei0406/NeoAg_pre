import glob
import os
import argparse
import subprocess

def get_peploc (inputdir, outputdir):
	for filepath in glob.glob(os.path.join(inputdir,'*.SNP.uniq.25mer.txt')):
		#print(filepath)
		filename_tmp = filepath.split("/")[-1]
		filename = filename_tmp.split(".")[0] + ".pep.loc"
		#print(filename)
		with open ('{}/{}'.format (outputdir, filename), 'w') as fw:
			fw.writelines("identity"+"\t"+"peptide"+"\t"+"loc"+"\n")
			pep = open(filepath)
			for l in pep.readlines():
				if l.startswith("ID"):
					continue
				else:
					gene=l.strip().split("\t")
					#print (gene[8])
					can_AA=gene[8]
					if str(can_AA).startswith("_"):
						AA=str(gene[7]).split(":")[3]
						#print(AA)
						if ("," not in AA):
							AA_M=AA
							fw.writelines(str(gene[1])+":"+ AA_M +":"+str(gene[30])+":"+str(gene[0])+"\n")
							fw.writelines(str(gene[20]) + + '\t'+ str(gene[30]) + "\n")
						if ("," in AA):
							AA_tmp=AA.split(",")[0]
							AA_M=AA_tmp
							fw.writelines(">"+str(gene[1])+":"+ AA_M +":"+str(gene[30])+":"+str(gene[0])+"\n")
							fw.writelines(str(gene[20]) + "\n")
					else:
						AA_M = str(can_AA)
						AA = str (gene[7]).split (":")[3]
						if ("," not in AA):
							if (AA==AA_M):
								fw.writelines (str (gene[1]) + ":" + AA_M + "\t"+ str (gene[20]) + "\t" +str(gene[30]) +"\n")
							else:
								fw.writelines (str (gene[1]) + ":" + AA + "\t" + str (gene[20]) + "\t" + str (gene[30]) + "\n")
						else:
							if AA_M in str (gene[7]):
								fw.writelines (str (gene[1]) + ":" + AA_M + "\t" + str (gene[20]) + "\t" + str (gene[30]) + "\n")
							else:
								AA_tmp = AA.split (",")[0]
								fw.writelines ( str (gene[1]) + ":" + AA_tmp + "\t" + str (gene[20]) + "\t" + str (gene[30]) + "\n")
		fw.close()

def get_classIfasta (inputdir, outputdir):
	for filepath in glob.glob(os.path.join(inputdir,'*.SNP.uniq.25mer.txt')):
		#print(filepath)
		filename_tmp = filepath.split("/")[-1]
		filename = filename_tmp.split(".")[0] + ".classI.fasta"
		#print(filename)
		with open ('{}/{}'.format (outputdir, filename), 'w') as fw:
			pep = open(filepath)
			for l in pep.readlines():
				if l.startswith("ID"):
					continue
				else:
					gene=l.strip().split("\t")
					#print (gene[8])
					can_AA=gene[8]
					if str(can_AA).startswith("_"):

						AA=str(gene[7]).split(":")[3]
						#print(AA)
						if ("," not in AA):
							AA_M=AA
							fw.writelines(">"+str(gene[1])+":"+ AA_M +"\n")
							fw.writelines(str(gene[20]) + "\n")
						if ("," in AA):
							AA_tmp=AA.split(",")[0]
							AA_M=AA_tmp
							fw.writelines(">"+str(gene[1])+":"+ AA_M +"\n")
							fw.writelines(str(gene[20]) + "\n")
					else:
						AA_M = str(can_AA)
						AA = str (gene[7]).split (":")[3]
						if ("," not in AA):
							if (AA==AA_M):
								fw.writelines(">"+str(gene[1])+":"+ AA_M +"\n")
								fw.writelines (str (gene[20]) + "\n")
							else:
								fw.writelines (">" + str (gene[1]) + ":" + AA +  "\n")
								fw.writelines (str (gene[20]) + "\n")
						else:
							if AA_M in str (gene[7]):
								fw.writelines (">" + str (gene[1]) + ":" + AA_M + "\n")
								fw.writelines (str (gene[20]) + "\n")
							else:
								AA_tmp = AA.split (",")[0]
								fw.writelines (">" + str (gene[1]) + ":" + AA_tmp + "\n")
								fw.writelines(str(gene[20]) + "\n")
		fw.close()


def get_classIIfasta(inputdir,outputdir):
	for filepath in glob.glob (os.path.join (inputdir, '*.SNP.uniq.25mer.txt')):
		print (filepath)
		filename_tmp = filepath.split ("/")[-1]
		filename = filename_tmp.split (".")[0] + ".classII.fasta"
		#print (filename)
		with open ('{}/{}'.format (outputdir, filename), 'w') as fw:
			pep = open (filepath)
			for l in pep.readlines ():
				if l.startswith ("ID"):
					continue
				else:
					gene = l.strip ().split ("\t")
					if (len(str(gene[19]))<15):
						continue
					else:

						#print (gene[8])
						can_AA = gene[8]
						if str (can_AA).startswith ("_"):

							AA = str (gene[7]).split (":")[3]
							#print (AA)
							if ("," not in AA):
								AA_M = AA
								fw.writelines (">" + str (gene[1]) + ":" + AA_M + "\n")
								fw.writelines (str (gene[19]) + "\n")
							if ("," in AA):
								AA_tmp = AA.split (",")[0]
								AA_M = AA_tmp
								fw.writelines (">" + str (gene[1]) + ":" + AA_M + "\n")
								fw.writelines (str (gene[19]) + "\n")
						else:
							AA_M = str (can_AA)
							AA = str (gene[7]).split (":")[3]
							if ("," not in AA):
								if (AA == AA_M):
									fw.writelines (">" + str (gene[1]) + ":" + AA_M + "\n")
									fw.writelines (str (gene[19]) + "\n")
								else:
									fw.writelines (
									">" + str (gene[1]) + ":" + AA +  "\n")
									fw.writelines (str (gene[19]) + "\n")
							else:
								if AA_M in str (gene[7]):
									fw.writelines (">" + str (gene[1]) + ":" + AA_M + "\n")
									fw.writelines (str (gene[19]) + "\n")
								else:
									AA_tmp = AA.split (",")[0]
									fw.writelines (">" + str (gene[1]) + ":" + AA_tmp + "\n")
									fw.writelines (str (gene[19]) + "\n")
		fw.close ()


def main():
	input_dir = os.getcwd()
	print("Work directionary:")
	print(input_dir)
	input_dir_VEP = input_dir + '/data/VEP/25mer/'
	out_dir_pep=input_dir + '/data/neoantigen/pep_fasta/'
	if os.path.isdir (input_dir_VEP):
		pass
	else:
		subprocess.call ('mkdir -p ' + input_dir_VEP, shell=True)
	if os.path.isdir (out_dir_pep):
		pass
	else:
		subprocess.call ('mkdir -p ' + out_dir_pep, shell=True)
	get_classIIfasta(input_dir_VEP,out_dir_pep)
	get_classIfasta (input_dir_VEP, out_dir_pep)
	get_peploc(input_dir_VEP, out_dir_pep)


if __name__ == '__main__':
	main ()

'''
python /data/cyu/scripts/wes/4_neoantigen/vep2fasta.py

python /data/cyu/scripts/wes/4_neoantigen/vep2fasta.py
'''
