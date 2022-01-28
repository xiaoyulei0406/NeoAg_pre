import os, sys, subprocess
import argparse
import glob


PyClone='/home/ec2-user/anaconda3/envs/pyclone/bin/PyClone'

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='Run PyClone')
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

    sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
    log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed

    fout = sh_dir + 'Lohhla_' + start + '_' + end + '.sh'
    sh = open (fout, 'w')
    pts = open (pt_file)  # note there is a header
    log = log_dir + 'Lohhla_' + start + '_' + end + '.log'
    lns = pts.readlines ()
    subprocess.call ('sudo mkdir -p ' + out_dir + '/data/lohhla/input_bam/', shell=True)
    subprocess.call ('sudo chmod 777 ' + out_dir + '/data/lohhla/input_bam/', shell=True)
    subprocess.call ('sudo chown -R ec2-user:ec2-user ' + out_dir + '/data/lohhla/', shell=True)
    cmd = ''

    for i in range ((int) (start), (int) (end) + 1):
        print (lns[i])
        tmp= lns[i].strip("\n").split (',')
        normal_sample_id = tmp[0]  # patient ID
        tumor_sample_id = tmp[1]
        print (normal_sample_id)
        subprocess.call ('sudo mkdir -p ' + out_dir + '/data/lohhla/input_bam/'+ tumor_sample_id +'/', shell=True)
        subprocess.call ('sudo mkdir -p ' + out_dir + '/data/lohhla/output/' + tumor_sample_id + '/', shell=True)
        subprocess.call ('sudo chmod 777 ' + out_dir + '/data/lohhla/output/', shell=True)
        input_normal_bam = input_dir + normal_sample_id + '.recal.bam'
        out_normal_bam= out_dir + '/data/lohhla/input_bam/'+ tumor_sample_id + '/' + normal_sample_id +'_recal_chr.bam'
        out_normal_bai = out_dir + '/data/lohhla/input_bam/'+ tumor_sample_id +'/'  + normal_sample_id + '_recal_chr.bam.bai'
        input_tumor_bam = input_dir + tumor_sample_id + '.recal.bam'
        out_tumor_bam = out_dir + '/data/lohhla/input_bam/'+ tumor_sample_id + '/' + tumor_sample_id +'_recal_chr.bam'
        out_tumor_bai = out_dir + '/data/lohhla/input_bam/'+ tumor_sample_id +'/'  + tumor_sample_id + '_recal_chr.bam.bai'

        cmd += 'echo Start adding chr to normal bam files \n\n'
        #cmd += 'sh /data/cyu/scripts/wes/hlatyping/samtools4chr.sh ' + input_normal_bam + ' ' + out_normal_bam + ' ' + out_normal_bai
        cmd += '/home/ec2-user/anaconda3/bin/samtools view -H ' + input_normal_bam + ' | '
        cmd +='sed -e \'s/SN:\\([0-9XY]\\)/SN:chr\\1/\' -e \'s/SN:MT/SN:chrM/\' | '
        cmd +='/home/ec2-user/anaconda3/bin/samtools reheader - ' + input_normal_bam + ' > '+ out_normal_bam
        cmd += ' \n\n'
        cmd +='/home/ec2-user/anaconda3/bin/samtools index '+ out_normal_bam + ' > '+ out_normal_bai
        cmd += '\n\n'


        cmd += '\n\n'
        cmd += 'echo Start adding chr to tumor bam files \n\n'
        #cmd += 'sh /data/cyu/scripts/wes/hlatyping/samtools4chr.sh ' + input_tumor_bam + ' '+ out_tumor_bam + ' ' + out_tumor_bai
        #cmd += '\n\n'
        cmd += '/home/ec2-user/anaconda3/bin/samtools view -H ' + input_tumor_bam + ' | '
        cmd +='sed -e \'s/SN:\\([0-9XY]\\)/SN:chr\\1/\' -e \'s/SN:MT/SN:chrM/\' | '
        cmd +='/home/ec2-user/anaconda3/bin/samtools reheader - ' + input_tumor_bam + ' > '+ out_tumor_bam
        cmd += ' \n\n'
        cmd +='/home/ec2-user/anaconda3/bin/samtools index '+ out_tumor_bam + ' > '+ out_tumor_bai
        cmd += '\n\n'

        cmd += 'echo Running lohhla ... \n\n'
        cmd +='docker run -v ' + out_dir + ':/data/ pici/lohhla:1.0 Rscript /root/lohhla/LOHHLAscript.R '
        cmd +='--patientId ' + tumor_sample_id + ' --outputDir /data/data/lohhla/output/' + tumor_sample_id +'/ '
        cmd +='--normalBAMfile /data/data/lohhla/input_bam/' + tumor_sample_id + '/'+ normal_sample_id +'_recal_chr.bam '
        cmd +='--BAMDir ' + '/data/data/lohhla/input_bam/' + tumor_sample_id + '/ '
        cmd +='--hlaPath /data/data/HLAtyping/summary/' + tumor_sample_id+ '.hla.consensus.txt '
        cmd +='--HLAfastaLoc /data/refs/hla_abc.fasta '
        cmd +='--HLAexonLoc /data/refs/hla.dat '
        cmd +='--CopyNumLoc /data/data/sequenza/' + tumor_sample_id + '/solutions.txt '
        cmd +='--mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE '
        cmd +='--gatkDir /picard --novoDir /novocraft/ '
        cmd +='\n\n'

    sh.write (cmd)
    sh.close ()
    pts.close ()

    cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
    subprocess.call (cmd, shell=True)

if __name__ == '__main__':
    main ()

'''
python /data/cyu/scripts/wes/3_hlatyping/run_lohhla.py \
-i /data/cyu/gbm/data/bam/ \
-p /data/cyu/gbm/sample_pairs.txt \
-s 1 -e 1 -o /data/cyu/gbm/

#run with ec2-user
'''
