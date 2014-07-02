library(lineprof)
library(seqtools)
library(affinity)
library(toxicity)
library(targetsurveyorr)
library(WriteXLS)
setwd('~/Documents//Santaris//Assignment')


col_reformat <- function(cname){
  NC <- length(cname)
  c(1:5,8:11,NC,NC-2,12:17,7,6,NC-1,grep('_d',cname),grep('site_present_in',cname),18:22)
}

################################################
#IRX3 codes for a protein that regulates other genes, and is present both in and outside the 
#brain, in organs and cells such as fat cells.
#IRX3 <- 'ENSG00000177508'
#IRX3_survey <- make_targetsurvey_for(IRX3)
#save(list=c('IRX3','IRX3_survey'),file='~/Documents/Santaris/Assignment/IRX3.RData')
load('~/Documents/Santaris/Assignment/IRX3.RData')
IRX3_select <- select_oligo_library(annotated_sites = IRX3_survey$annotated_sites,site_present_in = 'mmusculus',n=list(exon=50),dg_opt = -20,w_dg = 1,w_rf = 1)

IRX3_index <- col_reformat(colnames(IRX3_survey$annotated_sites))
file_out1 <- as.data.frame(IRX3_select[,c(IRX3_index,length(IRX3_index)+1:4)])
file_out2 <- as.data.frame(IRX3_survey$annotated_sites)[,IRX3_index]
file_name = '~/Documents/Santaris//Assignment/IRX3_oligos.xls'
WriteXLS(c('file_out1','file_out2'),
         ExcelFileName = file_name,row.names = F,
         SheetNames=c('50_selected_oligos','targetsurvey'))
rm(list=c('file_out1','file_out2','file_name'))
rm(list=ls()[grep('IRX',ls())])

################################################
#Smad7 <- 'ENSG00000101665'
#Smad7_survey <- make_targetsurvey_for(Smad7)
#save(list=c('Smad7','Smad7_survey'),file='Smad7.RData')
load('~/Documents/Santaris/Assignment/Smad7.RData')
Smad7_select <- select_oligo_library(annotated_sites = Smad7_survey$annotated_sites,
                                     n=list(exon=25, intron=25),site_present_in = 'mmusculus')

Smad7_index <- col_reformat(colnames(Smad7_survey$annotated_sites))
file_out1 <- as.data.frame(Smad7_select[,c(Smad7_index,length(Smad7_index)+1:4)])
file_out2 <- as.data.frame(Smad7_survey$annotated_sites)[,Smad7_index]
file_name = '~/Documents/Santaris//Assignment/SMAD7_oligos.xls'
WriteXLS(c('file_out1','file_out2'),
         ExcelFileName = file_name,row.names = F,
         SheetNames=c('50_selected_oligos','targetsurvey'))
rm(list=c('file_out1','file_out2','file_name'))
rm(list=ls()[grep('Smad',ls())])

################################################
#NFKB1 <- 'ENSG00000109320' 
#NFKB1_survey <- make_targetsurvey_for(NFKB1)
#save(list=c('NFKB1','NFKB1_survey'),file='NFKB1.RData')
#load('~/Documents/Santaris/Assignment/NFKB1.RData')
NFKB1_select <- select_oligo_library(annotated_sites = NFKB1_survey$annotated_sites,
                                     n=list(exon=25, intron=25),site_present_in = 'mmusculus')
NFKB1_index <- col_reformat(colnames(NFKB1_survey$annotated_sites))
file_out1 <- as.data.frame(NFKB1_select[,c(NFKB1_index,length(NFKB1_index)+1:4)])
N <- nrow(NFKB1_survey$annotated_sites)
file_out2_1 <- as.data.frame(NFKB1_survey$annotated_sites)[1:(N/4),NFKB1_index]
file_out2_2 <- as.data.frame(NFKB1_survey$annotated_sites)[(N/4+1):(N/2),NFKB1_index]
file_out2_3 <- as.data.frame(NFKB1_survey$annotated_sites)[(N/2+1):(N*3/4),NFKB1_index]
file_out2_4 <- as.data.frame(NFKB1_survey$annotated_sites)[(N*3/4+1):N,NFKB1_index]
file_name = '~/Documents/Santaris//Assignment/NFKB_oligos.xls'
WriteXLS(c('file_out1','file_out2_1','file_out2_2','file_out2_3','file_out2_4'),
         ExcelFileName = file_name,row.names = F,
         SheetNames=c('50_selected_oligos','targetsurvey_1','targetsurvey_2',
                      'targetsurvey_3','targetsurvey_4'))
rm(list=c('file_out1','file_out2_1','file_out2_2','file_out2_3','file_out2_4','file_name'))
rm(list=ls()[grep('NFK',ls())])

################################################
#P65 <- 'ENSG00000173039'
#P65_survey <- make_targetsurvey_for(P65)
#save(list=c('P65','P65_survey'),file='P65.RData')
load('~/Documents/Santaris/Assignment/P65.RData')
P65_select <- select_oligo_library(annotated_sites = P65_survey$annotated_sites,
                                     n=list(exon=49, intron=1),site_present_in = 'mmusculus')

P65_index <- col_reformat(colnames(P65_survey$annotated_sites))
file_out1 <- as.data.frame(P65_select[,c(P65_index,length(P65_index)+1:4)])
file_out2 <- as.data.frame(P65_survey$annotated_sites)[,P65_index]
file_name <- '~/Documents/Santaris//Assignment/P65_oligos.xls'
WriteXLS(c('file_out1','file_out2'),
         ExcelFileName = file_name,row.names = F,
         SheetNames=c('50_selected_oligos','targetsurvey'))
rm(list=c('file_out1','file_out2','file_name'))
rm(list=ls()[grep('P65',ls())])

################################################
#ALDH1A1 <- 'ENSG00000165092'
#ALDH_survey <- make_targetsurvey_for(ALDH1A1)
#save(list=c('ALDH1A1','ALDH_survey'),file='ALDH1A1.RData')
load('~/Documents/Santaris/Assignment/ALDH1A1.RData')
ALDH_select <- select_oligo_library(annotated_sites = ALDH_survey$annotated_sites,
                                    n=list(exon=50),site_present_in = 'mmusculus')

ALDH_index <- col_reformat(colnames(ALDH_survey$annotated_sites))
file_out1 <- as.data.frame(ALDH_select[,c(ALDH_index,length(ALDH_index)+1:4)])
N <- nrow(ALDH_survey$annotated_sites)
file_out2_1 <- as.data.frame(ALDH_survey$annotated_sites)[1:(N/6),ALDH_index]
file_out2_2 <- as.data.frame(ALDH_survey$annotated_sites)[(N/6+1):(N*2/6),ALDH_index]
file_out2_3 <- as.data.frame(ALDH_survey$annotated_sites)[(N*2/6+1):(N*3/6),ALDH_index]
file_out2_4 <- as.data.frame(ALDH_survey$annotated_sites)[(N*3/6+1):(N*4/6),ALDH_index]
file_out2_5 <- as.data.frame(ALDH_survey$annotated_sites)[(N*4/6+1):(N*5/6),ALDH_index]
file_out2_6 <- as.data.frame(ALDH_survey$annotated_sites)[(N*5/6+1):N,ALDH_index]
file_name <- '~/Documents/Santaris/Assignment/ALDH_oligos.xls'
WriteXLS(c('file_out1','file_out2_1','file_out2_2','file_out2_3','file_out2_4',
           'file_out2_5','file_out2_6'),
         ExcelFileName = file_name,row.names = F,
         SheetNames=c('50_selected_oligos','targetsurvey_1','targetsurvey_2',
                      'targetsurvey_3','targetsurvey_4','targetsurvey_5','targetsurvey_6'))
rm(list=c('file_out1','file_out2_1','file_out2_2','file_out2_3','file_out2_4','file_out2_5',
          'file_out2_6','file_name'))
rm(list=ls()[grep('ALDH',ls())])

################################################
#Hepcidin <- 'ENSG00000105697'
#HAMP_survey <- make_targetsurvey_for(Hepcidin)
#save(list=c('Hepcidin','HAMP_survey'),file='~/Documents/Santaris/Assignment/HAMP.RData')
load('~/Documents/Santaris/Assignment/HAMP.RData')
HAMP_select <- select_oligo_library(annotated_sites = HAMP_survey$annotated_sites,
                                    n=list(exon=25, intron=25),site_present_in = NULL)

HAMP_index <- col_reformat(colnames(HAMP_survey$annotated_sites))
file_out1 <- as.data.frame(HAMP_select[,c(HAMP_index,length(HAMP_index)+1:4)])
file_out2 <- as.data.frame(HAMP_survey$annotated_sites)[,HAMP_index]
file_name <- '~/Documents/Santaris//Assignment/HAMP_oligos.xls'
WriteXLS(c('file_out1','file_out2'),
         ExcelFileName = file_name,row.names = F,
         SheetNames=c('50_selected_oligos','targetsurvey'))
rm(list=c('file_out1','file_out2','file_name'))
rm(list=ls()[grep('HAMP',ls())])