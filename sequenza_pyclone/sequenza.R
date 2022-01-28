library(sequenza)

args = commandArgs(trailingOnly = TRUE)

if(length(args) ==0){
  cat("This script run sequenza for one patient.\n\n")
  cat("Usage: Rscript sequenza.R patient <input_dir> \n")
  quit()
}

pt = args[1]
input_dir= args[2]


cat('###Processing')
cat(pt)
cat('#######\n\n')
data.file <-paste0( input_dir ,'/', pt , '.binning.seqz.gz')
cat('Running sequenza.extract....\n')
seqz <- sequenza.extract(data.file, chromosome.list=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
cat('sequenza.extract done. \n\n')


cat('Running sequenza.fit...\n\n')
CP <- sequenza.fit(seqz)
cat('sequenza.fit done.\n\n')

#out_dir <- paste0(input_dir,'/plot/', pt)
#cat('Making an directory for the patient\n\n')
#system(paste0('mkdir -p ', out_dir))
#cmd = paste0('mv ',input_dir, '/', pt, '*.seq*.gz* ', out_dir)

#system(cmd)

cat('Directory has been made. \n')


cat('Running sequenza.results...\n\n')
sequenza.results(sequenza.extract = seqz, cp.table = CP, sample.id = pt, out.dir = input_dir)
cat('sequenza.results done.\n\n')
cat('The result is at ')
cat(out_dir)
cat('\n\n')

cat('#########Done with patient ')
cat(pt)
cat(' #########\n\n')
