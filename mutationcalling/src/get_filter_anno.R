#!/usr/bin/env Rscript
suppressMessages(library("optparse"))
optionList = list(make_option(c("-i", "--input_dir"), type="character", default=NULL, help="inputdir: file or working folder.\n"),
                  make_option(c("-p", "--patientID"), type="character", default=NULL, help="patient ID information"),
                  make_option(c("-c", "--caller"), type="character", default=NULL, help="caller: Varscan2,Mutect2,Mutect,Strelka,Muse"),
                  make_option(c("-o","--outdir"),type="character",default=NULL,help="output folder")
                  );

opt_parser = OptionParser(option_list=optionList)
opt = parse_args(opt_parser)

if (is.null(opt$input_dir)||is.null(opt$outdir)||is.null(opt$caller)){
  print_help(opt_parser)
  stop("Parameter errors!", call.=FALSE)
}

suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(gtable))
suppressMessages(library(maftools))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))

fileDir=opt$input_dir
outdir=opt$outdir
caller=opt$caller

merged<-data.frame(matrix(ncol=62, nrow=0))
mfMerged<-data.frame(matrix(ncol=63, nrow=0))
colnames(merged)<-c('Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene',
                    'AAChange.refGene','cytoBand','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE',
                    'ExAC_OTH','ExAC_SAS','avsnp147','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
                    'Polyphen2_HVAR_score','Polyphen2_HVAR_pred','LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred',
                    'MutationAssessor_score','MutationAssessor_pred','FATHMM_score','FATHMM_pred','PROVEAN_score','PROVEAN_pred',
                    'VEST3_score','CADD_raw','CADD_phred','DANN_score','fathmm.MKL_coding_score','fathmm.MKL_coding_pred','MetaSVM_score',
                    'MetaSVM_pred','MetaLR_score','MetaLR_pred','integrated_fitCons_score','integrated_confidence_value',
                    'GERP.._RS','phyloP7way_vertebrate','phyloP20way_mammalian','phastCons7way_vertebrate',
                    'phastCons20way_mammalian','SiPhy_29way_logOdds','tumor_vaf','tumor_ref_count','tumor_alt_count','normal_vaf',
                    'normal_ref_count','normal_alt_count','Tumor_Sample_Barcode')
colnames(mfMerged)<-c('Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2',
                      'Hugo_Symbol','Variant_Classification','tx','exon','txChange','aaChange','Variant_Type','Func63025520.refGene',
                      'Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand','snp142','ExAC_ALL',
                      'ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','avsnp147','SIFT_score',
                      'SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_pred',
                      'LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred','MutationAssessor_score',
                      'MutationAssessor_pred','FATHMM_score','FATHMM_pred','PROVEAN_score','PROVEAN_pred','VEST3_score',
                      'CADD_raw','CADD_phred','DANN_score','fathmm-MKL_coding_score','fathmm-MKL_coding_pred','MetaSVM_score',
                      'MetaSVM_pred','MetaLR_score','MetaLR_pred','integrated_fitCons_score','integrated_confidence_value',
                      'GERP++_RS','phyloP7way_vertebrate','phyloP20way_mammalian','phastCons7way_vertebrate',
                      'phastCons20way_mammalian','SiPhy_29way_logOdds')

get_mergeSNPIndels<-function(annovar_INDEL_Files,annovar_SNP_Files,vcf_INDEL_Files,vcf_SNP_Files,caller){
  ##check if dir exists
  sub_dir<-"filtered"
  output_dir1<- file.path(fileDir, sub_dir)
  if (!dir.exists(output_dir1)){
    dir.create(output_dir1)
  } else {
    print("Dir already exists!")
  }

  ##merge SNP and INDEL annovar files first
  for(annovar_INDEL_File in annovar_INDEL_Files){
    annovar_INDEL_NameTMP<- unlist(strsplit(annovar_INDEL_File, split = '/'))[length(unlist(strsplit(annovar_INDEL_File, split = '/')))]
    annovar_INDEL_Name=gsub(".indels.refGene.hg19_multianno.txt", "", annovar_INDEL_NameTMP)
#    print(annovar_INDEL_Name)
    for(annovar_SNP_File in annovar_SNP_Files){
      annovar_SNP_NameTMP<- unlist(strsplit(annovar_SNP_File, split = '/'))[length(unlist(strsplit(annovar_SNP_File, split = '/')))]
      annovar_SNP_Name=gsub(".snps.refGene.hg19_multianno.txt", "", annovar_SNP_NameTMP)
      if (annovar_SNP_Name==annovar_INDEL_Name){
#        print(annovar_SNP_Name)
        annovar_SNPF=read.table(annovar_SNP_File,skip = 1, header = FALSE,sep ="\t")
        ano_SNPF<-annovar_SNPF[,c(1:55)]
        annovar_INDELF=read.table(annovar_INDEL_File,skip = 1, header = FALSE,sep = "\t")
        ano_INDELF<-annovar_INDELF[,c(1:55)]
        merged_annovar<-rbind(ano_SNPF,ano_INDELF)
        names(merged_annovar) <- c('Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene',
                                   'AAChange.refGene','cytoBand','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE',
                                   'ExAC_OTH','ExAC_SAS','avsnp147','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
                                   'Polyphen2_HVAR_score','Polyphen2_HVAR_pred','LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred',
                                   'MutationAssessor_score','MutationAssessor_pred','FATHMM_score','FATHMM_pred','PROVEAN_score','PROVEAN_pred',
                                   'VEST3_score','CADD_raw','CADD_phred','DANN_score','fathmm.MKL_coding_score','fathmm.MKL_coding_pred','MetaSVM_score',
                                   'MetaSVM_pred','MetaLR_score','MetaLR_pred','integrated_fitCons_score','integrated_confidence_value',
                                   'GERP.._RS','phyloP7way_vertebrate','phyloP20way_mammalian','phastCons7way_vertebrate',
                                   'phastCons20way_mammalian','SiPhy_29way_logOdds')
        write.table(merged_annovar,paste(fileDir,annovar_SNP_Name,".refGene.hg19_multianno.merged.txt",sep = ""),sep="\t",quote=F,row.names = F)
        
      }
    }
  }
  ##merge SNP and INDEL depth files first
  for (vcf_INDEL_File in vcf_INDEL_Files){
    vcf_INDEL_NameTMP<- unlist(strsplit(vcf_INDEL_File, split = '/'))[length(unlist(strsplit(vcf_INDEL_File, split = '/')))]
    vcf_INDEL_Name=gsub(".indels.adddepth.vcf", "", vcf_INDEL_NameTMP)
#    print(vcf_INDEL_Name)
    for(vcf_SNP_File in vcf_SNP_Files){
      vcf_SNP_NameTMP<- unlist(strsplit(vcf_SNP_File, split = '/'))[length(unlist(strsplit(vcf_SNP_File, split = '/')))]
      vcf_SNP_Name=gsub(".snps.adddepth.vcf", "", vcf_SNP_NameTMP)
      if (vcf_SNP_Name==vcf_INDEL_Name){
#        print(vcf_SNP_Name)
        vcf_SNPF=read.table(vcf_SNP_File,header =F, sep ="\t")
        vcf_INDELF=read.table(vcf_INDEL_File,header =F,sep = "\t")
        merged_vcf<-rbind(vcf_SNPF,vcf_INDELF)
        write.table(merged_vcf,paste(fileDir,vcf_SNP_Name,".adddepth.merged.vcf",sep = ""),sep="\t",quote=F,row.names = F)
        
      }
    }
  }  
  
  vcfFiles<-list.files(fileDir, pattern = paste0("*.",caller,".adddepth.merged.vcf",sep=""), full.names = TRUE) 
  annovarFiles<-list.files(fileDir, pattern = paste0("*.",caller,".refGene.hg19_multianno.merged.txt",sep=""), full.names = TRUE) 
  
  for(annovarFile in annovarFiles){
    annovarNameTMP<- unlist(strsplit(annovarFile, split = '/'))[length(unlist(strsplit(annovarFile, split = '/')))]
    annovarName=gsub(".refGene.hg19_multianno.merged.txt", "", annovarNameTMP)
    for (vcfFile in vcfFiles){
      vcfNameTMP<- unlist(strsplit(vcfFile, split = '/'))[length(unlist(strsplit(vcfFile, split = '/')))]
      vcfName = gsub(".adddepth.merged.vcf", "", vcfNameTMP)
      caller1=paste(".",caller,sep="")
      vcfName1=gsub(caller1, "", vcfName)
#      print(vcfName)
#      print(annovarName)
      if (vcfName==annovarName){
        annovarF=read.table(annovarFile,header =T, sep ="\t")
        vcfF=read.table(vcfFile,,header =T,sep = "\t")
        print(annovarName)
        df<-cbind(annovarF,vcfF)
        
        df_new<-df[,c(1:55,61:66)]
        names(df_new)[56]<-"tumor_vaf"
        names(df_new)[57]<-"tumor_ref_count"
        names(df_new)[58]<-"tumor_alt_count"
        names(df_new)[59]<-"normal_vaf"
        names(df_new)[60]<-"normal_ref_count"
        names(df_new)[61]<-"normal_alt_count"
        df_1<-dplyr::filter(df_new,tumor_vaf>=0.05)##Select out tumor VAF >=0.05
        df_2<-dplyr::filter(df_1,normal_alt_count<2)##Filter out normal alt count <2
        df_3<-dplyr::filter(df_2,tumor_alt_count>3)##select out normal alt count >3
        df_3$Tumor_Sample_Barcode=vcfName
        df_4<-dplyr::filter(df_3,Func.refGene=='exonic')##select out only exonic
        df_5<-dplyr::filter(df_4,ExonicFunc.refGene!='synonymous SNV') ##delete synonymous SNV
        df_6<-dplyr::filter(df_5,ExonicFunc.refGene!='unknown')##delete unknow
        df_7<-dplyr::filter(df_6,tumor_alt_count + tumor_ref_count >10) ##select out mutations with coverage >10
        df_8<-dplyr::filter(df_7,normal_alt_count + normal_ref_count >10) 
        #add tumor sample barcode
        df_8$Tumor_Sample_Barcode=vcfName1
        write.table(df_8,paste(fileDir,"/filtered/",vcfName,".merged.refGene.hg19_multianno.filtered.txt",sep = ""),sep="\t",quote=F,row.names = F)
        merged<-rbind(merged, df_8)
        
      }
    }
  }
  write.table(merged,paste(fileDir,"/filtered/",caller,".merged.refGene.hg19_multianno.filtered.all.txt",sep = ""),sep="\t",quote=F,row.names = F)
  fileFiltered_Dir=paste(fileDir,"/filtered/",sep="")
  annovar_filteredFiles<-list.files(fileFiltered_Dir, pattern = paste0("*.",caller,".merged.refGene.hg19_multianno.filtered.txt",sep=""), full.names = TRUE) 
  
  for(annovarFile in annovar_filteredFiles){
        annovarNameTMP<- unlist(strsplit(annovarFile, split = '/'))[length(unlist(strsplit(annovarFile, split = '/')))]
        print(annovarName)
        annovarName=gsub(paste0(".",caller,".merged.refGene.hg19_multianno.filtered.txt",sep=""), "", annovarNameTMP)
        annovar_maf = lapply(annovarFile, annovarToMaf) #convert to MAFs using annovarToMaf
        annovar_maf = data.table::rbindlist(l = annovar_maf, fill = TRUE) #Merge into single MAF
        write.table(merged,paste(fileDir,"/filtered/",annovarName,".",caller,".merged.refGene.hg19_multianno.filtered.maf.txt",sep = ""),sep="\t",quote=F,row.names = F)
        
        ##get indels mafs
        DEL<-annovar_maf[which(with(annovar_maf,Variant_Type=='DEL'))]
        INS<-annovar_maf[which(with(annovar_maf,Variant_Type=='INS'))]
        INDEL<-rbind(INS,DEL)
        ##get SNV mafs
        SNP<-annovar_maf[which(with(annovar_maf,Variant_Type=='SNP'))]
        
        write.table(SNP,paste(outdir,annovarName,".",caller,".SNPs.txt",sep = ""),sep="\t",quote=F,row.names = F)
        write.table(INDEL,paste(outdir,annovarName,".",caller,".INDELs.txt",sep = ""),sep="\t",quote=F,row.names = F)
  }
}
get_mergeSNPs<-function(annovar_SNP_Files,vcf_SNP_Files,caller){
  sub_dir<-"filtered"
  output_dir1<- file.path(fileDir, sub_dir)
  if (!dir.exists(output_dir1)){
    dir.create(output_dir1)
  } else {
    print("Dir already exists!")
  }
  ##combine SNP annovar and adddepth files first
  
  for(annovarFile in annovar_SNP_Files){
#    print(annovarFile)
    annovarNameTMP<- unlist(strsplit(annovarFile, split = '/'))[length(unlist(strsplit(annovarFile, split = '/')))]
    annovarName=gsub(".snps.refGene.hg19_multianno.txt", "", annovarNameTMP)
    for (vcfFile in vcf_SNP_Files){
      vcfNameTMP<- unlist(strsplit(vcfFile, split = '/'))[length(unlist(strsplit(vcfFile, split = '/')))]
      vcfName = gsub(".snps.adddepth.vcf", "", vcfNameTMP)
      caller1=paste(".",caller,sep="")
      vcfName1=gsub(caller1, "", vcfName)

      if (vcfName==annovarName){
#        print(annovarFile)
        annovar_SNPF=read.table(annovarFile,skip = 1, header = FALSE,sep ="\t")
        annovarF<-annovar_SNPF[,c(1:55)]
        

        names(annovarF) <- c('Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene',
                                   'AAChange.refGene','cytoBand','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE',
                                   'ExAC_OTH','ExAC_SAS','avsnp147','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
                                   'Polyphen2_HVAR_score','Polyphen2_HVAR_pred','LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred',
                                   'MutationAssessor_score','MutationAssessor_pred','FATHMM_score','FATHMM_pred','PROVEAN_score','PROVEAN_pred',
                                   'VEST3_score','CADD_raw','CADD_phred','DANN_score','fathmm.MKL_coding_score','fathmm.MKL_coding_pred','MetaSVM_score',
                                   'MetaSVM_pred','MetaLR_score','MetaLR_pred','integrated_fitCons_score','integrated_confidence_value',
                                   'GERP.._RS','phyloP7way_vertebrate','phyloP20way_mammalian','phastCons7way_vertebrate',
                                   'phastCons20way_mammalian','SiPhy_29way_logOdds')
        vcfF=read.table(vcfFile,header =F,sep = "\t")

        df<-cbind(annovarF,vcfF)
        df_new<-df[,c(1:55,61:66)]
        names(df_new)[56]<-"tumor_vaf"
        names(df_new)[57]<-"tumor_ref_count"
        names(df_new)[58]<-"tumor_alt_count"
        names(df_new)[59]<-"normal_vaf"
        names(df_new)[60]<-"normal_ref_count"
        names(df_new)[61]<-"normal_alt_count"
        ##Select out tumor VAF >=0.05
        df_1<-dplyr::filter(df_new,tumor_vaf>=0.05)
        ##Filter out normal alt count <2
        df_2<-dplyr::filter(df_1,normal_alt_count<2)
        ##select out normal alt count >3
        df_3<-dplyr::filter(df_2,tumor_alt_count>3)
        ##select out only exonic
        df_4<-dplyr::filter(df_3,Func.refGene=='exonic')
        ##select out only exonic
        df_5<-dplyr::filter(df_4,ExonicFunc.refGene!='synonymous SNV')
        df_6<-dplyr::filter(df_5,ExonicFunc.refGene!='unknown')
        df_7<-dplyr::filter(df_6,tumor_alt_count + tumor_ref_count >10) ##select out mutations with coverage >10
        df_8<-dplyr::filter(df_7,normal_alt_count + normal_ref_count >10) 
        if (dim(df_8)[1] == 0){
          #write.table(df_8,paste(fileDir,"/filtered/",vcfName,".snps.refGene.hg19_multianno.filtered.txt",sep = ""),sep="\t",quote=F,row.names = F)
          #merged<-rbind(merged, df_8)
        }else{
          df_8$Tumor_Sample_Barcode=vcfName1
          print(annovarFile)
          write.table(df_8,paste(fileDir,"/filtered/",vcfName,".snps.refGene.hg19_multianno.filtered.txt",sep = ""),sep="\t",quote=F,row.names = F)
          merged<-rbind(merged, df_8)
        }
        #add tumor sample barcode
        
        
      }
    }
  }
  write.table(merged,paste(fileDir,"/filtered/",caller,".snps.refGene.hg19_multianno.filtered.all.txt",sep = ""),sep="\t",quote=F,row.names = F)
  fileFiltered_Dir<-paste(fileDir,"/filtered/",sep="")
  annovar_filteredFiles<-list.files(fileFiltered_Dir, pattern = paste0("*.",caller,".snps.refGene.hg19_multianno.filtered.txt",sep=""), full.names = TRUE) 
  for(annovarFile in annovar_filteredFiles){
        annovarNameTMP<- unlist(strsplit(annovarFile, split = '/'))[length(unlist(strsplit(annovarFile, split = '/')))]
        annovarName=gsub(".snps.refGene.hg19_multianno.filtered.txt", "", annovarNameTMP)
        annovar_maf = lapply(annovarFile, annovarToMaf) #convert to MAFs using annovarToMaf
        annovar_maf = data.table::rbindlist(l = annovar_maf, fill = TRUE) #Merge into single MAF
        ##get SNV mafs
        SNP<-annovar_maf[which(with(annovar_maf,Variant_Type=='SNP'))]
        write.table(SNP,paste(outdir,annovarName,".SNPs.txt",sep = ""),sep="\t",quote=F,row.names = F)
    }
}

get_merge<- function(annovar_Files,vcf_Files,caller){
  sub_dir<-"filtered"
  output_dir1<- file.path(fileDir, sub_dir)
  if (!dir.exists(output_dir1)){
    dir.create(output_dir1)
  } else {
    print("Dir already exists!")
  }
  ##combine SNP annovar and adddepth files first
  
  vcfFiles<-list.files(fileDir, pattern = paste0("*.",caller,".adddepth.vcf",sep=""), full.names = TRUE) 
  annovarFiles<-list.files(fileDir, pattern = paste0("*.",caller,".refGene.hg19_multianno.txt",sep=""), full.names = TRUE) 
  for(annovarFile in annovarFiles){
#    print(annovarFile)
    annovarNameTMP<- unlist(strsplit(annovarFile, split = '/'))[length(unlist(strsplit(annovarFile, split = '/')))]
    annovarName=gsub(".refGene.hg19_multianno.txt", "", annovarNameTMP)
    for (vcfFile in vcfFiles){
      vcfNameTMP<- unlist(strsplit(vcfFile, split = '/'))[length(unlist(strsplit(vcfFile, split = '/')))]
      vcfName = gsub(".adddepth.vcf", "", vcfNameTMP)
      caller1=paste(".",caller,sep="")
      vcfName1=gsub(caller1, "", vcfName)
      if (vcfName==annovarName){
#        print(vcfName)
        annovar_F=read.table(annovarFile,skip = 1, header = FALSE,sep ="\t")
        annovarF<-annovar_F[,c(1:55)]
        names(annovarF) <- c('Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene',
                             'AAChange.refGene','cytoBand','snp142','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE',
                             'ExAC_OTH','ExAC_SAS','avsnp147','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
                             'Polyphen2_HVAR_score','Polyphen2_HVAR_pred','LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred',
                             'MutationAssessor_score','MutationAssessor_pred','FATHMM_score','FATHMM_pred','PROVEAN_score','PROVEAN_pred',
                             'VEST3_score','CADD_raw','CADD_phred','DANN_score','fathmm.MKL_coding_score','fathmm.MKL_coding_pred','MetaSVM_score',
                             'MetaSVM_pred','MetaLR_score','MetaLR_pred','integrated_fitCons_score','integrated_confidence_value',
                             'GERP.._RS','phyloP7way_vertebrate','phyloP20way_mammalian','phastCons7way_vertebrate',
                             'phastCons20way_mammalian','SiPhy_29way_logOdds')
        vcfF=read.table(vcfFile,sep = "\t")
        df<-cbind(annovarF,vcfF)
        df_new<-df[,c(1:55,61:66)]
        names(df_new)[56]<-"tumor_vaf"
        names(df_new)[57]<-"tumor_ref_count"
        names(df_new)[58]<-"tumor_alt_count"
        names(df_new)[59]<-"normal_vaf"
        names(df_new)[60]<-"normal_ref_count"
        names(df_new)[61]<-"normal_alt_count"
        ##Select out tumor VAF >=0.05
        df_1<-dplyr::filter(df_new,tumor_vaf>=0.05)
        ##Filter out normal alt count <2
        df_2<-dplyr::filter(df_1,normal_alt_count<2)
        ##select out normal alt count >3
        df_3<-dplyr::filter(df_2,tumor_alt_count>3)
        df_3$Tumor_Sample_Barcode=vcfName
        ##select out only exonic
        df_4<-dplyr::filter(df_3,Func.refGene=='exonic')
        ##select out only exonic
        df_5<-dplyr::filter(df_4,ExonicFunc.refGene!='synonymous SNV')
        df_6<-dplyr::filter(df_5,ExonicFunc.refGene!='unknown')
        df_7<-dplyr::filter(df_6,tumor_alt_count + tumor_ref_count >10) ##select out mutations with coverage >10
        df_8<-dplyr::filter(df_7,normal_alt_count + normal_ref_count >10) 
        #add tumor sample barcode
        df_8$Tumor_Sample_Barcode=vcfName1
#        print(vcfName1)
        write.table(df_8,paste(fileDir,"/filtered/",vcfName,".refGene.hg19_multianno.filtered.txt",sep = ""),sep="\t",quote=F,row.names = F)
        merged<-rbind(merged, df_8)
        
      }
    }
  }
  write.table(merged,paste(fileDir,"/filtered/",caller,".merged.refGene.hg19_multianno.filtered.all.txt",sep = ""),sep="\t",quote=F,row.names = F)
  fileFiltered_Dir<-paste(fileDir,"/filtered/",sep="")
  annovar_filteredFiles<-list.files(fileFiltered_Dir, pattern = paste0("*.",caller,".refGene.hg19_multianno.filtered.txt",sep=""), full.names = TRUE) 
  for(annovarFile in annovar_filteredFiles){
    print(annovarFile)
    annovarNameTMP<- unlist(strsplit(annovarFile, split = '/'))[length(unlist(strsplit(annovarFile, split = '/')))]
    annovarName=gsub(".refGene.hg19_multianno.filtered.txt", "", annovarNameTMP)
    annovar_maf = lapply(annovarFile, annovarToMaf) #convert to MAFs using annovarToMaf
    annovar_maf = data.table::rbindlist(l = annovar_maf, fill = TRUE) #Merge into single MAF
    write.table(annovar_maf,paste(fileDir,"/filtered/",annovarName,".",caller,".merged.refGene.hg19_multianno.filtered.maf.txt",sep = ""),sep="\t",quote=F,row.names = F)

    ##get indels mafs
    DEL<-annovar_maf[which(with(annovar_maf,Variant_Type=='DEL'))]
    INS<-annovar_maf[which(with(annovar_maf,Variant_Type=='INS'))]
    INDEL<-rbind(INS,DEL)
    ##get SNV mafs
    SNP<-annovar_maf[which(with(annovar_maf,Variant_Type=='SNP'))]
    ##get MNP mafs
    MNP<-annovar_maf[which(with(annovar_maf,Variant_Type=='MNP'))]
    
    MNP$Callers="Mutect2"
    MNP$NCallers="1"
    write.table(SNP,paste(outdir,annovarName,".SNPs.txt",sep = ""),sep="\t",quote=F,row.names = F)
    write.table(INDEL,paste(outdir,annovarName,".INDELs.txt",sep = ""),sep="\t",quote=F,row.names = F)
    write.table(MNP,paste(outdir,annovarName,".MNPs.txt",sep = ""),sep="\t",quote=F,row.names = F)
  }
}

if (caller=="Varscan2" || caller=="Strelka"){
  annovar_SNP_Files <-list.files(fileDir, pattern = paste0("*",caller,".snps.refGene.hg19_multianno.txt",sep=""), full.names = TRUE)
  annovar_INDEL_Files <-list.files(fileDir, pattern = paste0("*",caller,".indels.refGene.hg19_multianno.txt",sep=""), full.names = TRUE)
#  print(annovar_SNP_Files)
  vcf_SNP_Files<-list.files(fileDir, pattern = paste0("*",caller,".snps.adddepth.vcf",sep=""), full.names = TRUE)
  vcf_INDEL_Files<-list.files(fileDir, pattern = paste0("*",caller,".indels.adddepth.vcf",sep=""), full.names = TRUE)
#  print(vcf_SNP_Files)
  get_mergeSNPIndels(annovar_INDEL_Files,annovar_SNP_Files,vcf_INDEL_Files,vcf_SNP_Files,caller)

}
if (caller=="Mutect" || caller=="Muse"){
  annovar_SNP_Files <-list.files(fileDir, pattern = paste0("*",caller,".snps.refGene.hg19_multianno.txt",sep=""), full.names = TRUE)

  vcf_SNP_Files<-list.files(fileDir, pattern = paste0("*",caller,".snps.adddepth.vcf",sep=""), full.names = TRUE)
  
  get_mergeSNPs(annovar_SNP_Files,vcf_SNP_Files,caller)
}
if (caller=="Mutect2"){
  annovarFiles <-list.files(fileDir, pattern = paste0("*",caller,".refGene.hg19_multianno.txt",sep=""), full.names = TRUE)
  
  vcfFiles<-list.files(fileDir, pattern = paste0("*",caller,".adddepth.vcf",sep=""), full.names = TRUE)
  
  get_merge(annovar_Files,vcf_Files,caller)
}

