library(tidyverse)
library(data.table)

dat <- read.csv('dataset/ADNIMERGE.csv')
dat2 <- dat %>% select(PTID, VISCODE, AGE, PTGENDER, PTEDUCAT, APOE4, 
                       ADAS13, RAVLT_immediate, RAVLT_learning, MMSE, CDRSB,
                       DX, Years_bl, Month, ORIGPROT) %>% arrange(PTID, Years_bl)
dat4 <- dat2

## filter out baseline observations: 1018 patients with MCI ##
bl.dat <- dat4[dat4$Years_bl==0, ]
bl.MCI.dat <- bl.dat[bl.dat$DX=='MCI', ] ## 1018 patients with MCI
sum(is.na(bl.MCI.dat$APOE4))  ## Remove 86 patients with missing APOE4 ##
sum(is.na(bl.MCI.dat$AGE))  ## Remove 1 patients with missing Age ##
bl.MCI.dat <- bl.MCI.dat[which(!is.na(bl.MCI.dat$APOE4)), ]
bl.MCI.dat <- bl.MCI.dat[which(!is.na(bl.MCI.dat$AGE)), ]
uID <- unique(bl.MCI.dat$PTID)

## filter out MCI patients with ID in baseline MCI dataset, 
MCI.dat <- dat4[dat4$PTID %in% uID, ]
max(MCI.dat$Month)/12

phase <- c('ADNI1', 'ADNIGO', 'ADNI2', 'ADNI3')
tab <- list()
for (i in 1:length(phase)){
  tmp.dat <- MCI.dat[which(MCI.dat$ORIGPROT==phase[i]), ]
  tab[[i]] <- table(tmp.dat$Month)
}

month <- (0:13)*12
month.ADNI1 <- sort(c(month, 6, 18))
month.ADNIGO2 <- sort(c(month, 6))
month.ADNI3 <- month

MCI.dat2 <- MCI.dat[0, ]
for (i in 1:length(uID)){
  tmp.dat <- MCI.dat[which(MCI.dat$PTID==uID[i]), ]
  tmp.prot <- tmp.dat$ORIGPROT[1]
  if (tmp.prot=='ADNI1'){
    MCI.dat2 <- rbind(MCI.dat2, tmp.dat[which(tmp.dat$Month %in% month.ADNI1), ])
  }
  if (tmp.prot=='ADNIGO' | tmp.prot=='ADNI2'){
    MCI.dat2 <- rbind(MCI.dat2, tmp.dat[which(tmp.dat$Month %in% month.ADNIGO2), ])
  }
  if (tmp.prot=='ADNI3'){
    MCI.dat2 <- rbind(MCI.dat2, tmp.dat[which(tmp.dat$Month %in% month.ADNI3), ])
  }
}

MCI.dat3 <- MCI.dat2[which(MCI.dat2$DX!=''), ]  ## 268 observations without diagnosis
table(MCI.dat3$Month)  ## 18 visits after Month 120
MCI.dat3 <- MCI.dat3[which(MCI.dat3$Month<=120), ]
bl.MCI.dat <- MCI.dat3[which(MCI.dat3$Month==0), ]
uID <- unique(MCI.dat3$PTID)  

## Use diagnosis to filter out survival dataset: 
## The event time is the time of first diagnosis of Dementia, otherwise is last observation time
surv.dat <- data.frame(ID = uID,
                       surv_time = rep(NA, length(uID)),
                       status = rep(NA, length(uID)),
                       Age = bl.MCI.dat$AGE,
                       Gender = bl.MCI.dat$PTGENDER,
                       Education = bl.MCI.dat$PTEDUCAT,
                       APOE4 = bl.MCI.dat$APOE4)
long.dat <- MCI.dat3[0, ]
for (i in 1:length(uID)){
  tmp.long <- MCI.dat3[MCI.dat3$PTID==uID[i], ]
  if (sum(tmp.long$DX=='Dementia')==0){ ## No dementia in later visits (Censoring indicator = 0)
    tmp.index <- nrow(tmp.long)
    tmp.surv.time <- tmp.long[tmp.index, ]$Years_bl
    tmp.status <- 0
  } else {  ## Dementia in later visits (Censoring indicator = 1)
    tmp.index <- min(which(tmp.long$DX=='Dementia'))
    tmp.surv.time <- tmp.long[tmp.index, ]$Years_bl
    tmp.status <- 1
  }
  long.dat <- rbind(long.dat, tmp.long[1:tmp.index, ])
  surv.dat$surv_time[i] <- tmp.surv.time
  surv.dat$status[i] <- tmp.status
}

surv.dat <- surv.dat %>% 
  mutate(Gender = ifelse(Gender=='Female', 1, 0)) %>% 
  select(ID, surv_time, status, Age, Gender, Education, APOE4)

## remove the dementia observation in longitudinal dataset ##
long.dat2 <- long.dat[long.dat$DX!='Dementia', ]
long.dat3 <- long.dat2 %>% mutate(ID = PTID,
                                  Time = Years_bl, 
                                  y1 = ADAS13, y2 = RAVLT_immediate, y3 = RAVLT_learning, y4 = MMSE, y5 = CDRSB) %>% 
  select(ID, Time, y1, y2, y3, y4, y5) 

## Keep subjects with valid baseline MRI data ##
adni1.mri <- read.csv('dataset/MRI_ADNI_1_GO_2.csv')
adni3.mri <- read.csv('dataset/MRI_ADNI3.csv')

mri.dat <- rbind(adni1.mri, adni3.mri) %>% mutate(visit_code='')

## represent each visit to the visit code, this can be seen in the image collections
for (i in 1:nrow(mri.dat)){
  tmp.visit <- mri.dat$Visit[i]
  if (tmp.visit==1) mri.dat$visit_code[i] <- 'sc' ## ADNI1 screening
  if (tmp.visit==2) mri.dat$visit_code[i] <- 'bl' ## ADNI1 baseline
  if (tmp.visit==3) mri.dat$visit_code[i] <- 'm06'
  if (tmp.visit==4) mri.dat$visit_code[i] <- 'm12'
  if (tmp.visit==5) mri.dat$visit_code[i] <- 'm18'
  if (tmp.visit==6) mri.dat$visit_code[i] <- 'm24'
  if (tmp.visit==7) mri.dat$visit_code[i] <- 'm30' ## no visit==7
  if (tmp.visit==8) mri.dat$visit_code[i] <- 'm36'
  if (tmp.visit==9) mri.dat$visit_code[i] <- 'm42' ## no visit==9
  if (tmp.visit==10) mri.dat$visit_code[i] <- 'm48'
  if (tmp.visit==11) mri.dat$visit_code[i] <- 'uns1'
  if (tmp.visit==13) mri.dat$visit_code[i] <- 'nv'
  
  if (tmp.visit==14) mri.dat$visit_code[i] <- 'scmri' ## ADNIGO screening
  if (tmp.visit==15) mri.dat$visit_code[i] <- 'm03'
  if (tmp.visit==17) mri.dat$visit_code[i] <- 'm60'
  
  if (tmp.visit==22) mri.dat$visit_code[i] <- 'v02' ## ADNI2 Screening MRI-New Pt (V02)
  if (tmp.visit==24) mri.dat$visit_code[i] <- 'v04' ## ADNI2 Month 3 MRI-New Pt (V04)
  if (tmp.visit==25) mri.dat$visit_code[i] <- 'v05' ## ADNI2 Month 6-New Pt (V05)
  if (tmp.visit==26) mri.dat$visit_code[i] <- 'v06' ## ADNI2 Initial Visit-Cont Pt (V06)
  if (tmp.visit==28) mri.dat$visit_code[i] <- 'v11' ## ADNI2 Year 1 Visit (V11)
  if (tmp.visit==30) mri.dat$visit_code[i] <- 'v21' ## ADNI2 Year 2 Visit (V21)
  if (tmp.visit==32) mri.dat$visit_code[i] <- 'v31' ## ADNI2 Year 3 Visit (V31)
  if (tmp.visit==34) mri.dat$visit_code[i] <- 'v41' ## ADNI2 Year 4 Visit (V41)
  if (tmp.visit==36) mri.dat$visit_code[i] <- 'v51' ## ADNI2 Year 5 Visit (V51)
  
  if (tmp.visit==101) mri.dat$visit_code[i] <- 'INIT' ## ADNI3 Initial Visit-Cont Pt (INIT)
  if (tmp.visit==102) mri.dat$visit_code[i] <- 'Y1' ## 	ADNI3 Year 1 Visit (Y1)
  if (tmp.visit==103) mri.dat$visit_code[i] <- 'Y2' ## 	ADNI3 Year 2 Visit (Y2)
}

## Capture baseline MRI information ##
bl <- mri.dat[which(mri.dat$visit_code %in% c('sc', 'bl', 'scmri', 'v02')), ]
bl.mci.mri <- bl[which(bl$Subject %in% uID), ] %>% select(visit_code, Subject, Acq.Date) %>% 
  distinct() %>% mutate(Date=as.Date(as.character(Acq.Date), "%m/%d/%Y")) %>% 
  select(-Acq.Date) %>% mutate(information=paste0(Subject, '_', Date)) %>% 
  arrange(Subject, desc(Date))
bl.mci.mri.nodup <- bl.mci.mri[which(!duplicated(bl.mci.mri$Subject)), ]
information <- bl.mci.mri.nodup$information

result_dir <- "E:/ADNI_result/all"
outdir <- "E:/MMFPCA_MRI/ADNI_MRI"

mri.files <- list.files(result_dir, full.names = T) ## Don't run
#save(list = 'mri.files', file = 'RData/mri.files.RData')

tmp <- str_extract(mri.files, '[0-9][0-9][0-9]_S_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]')
tmp2 <- mri.files[which(tmp %in% information)]
file.copy(tmp2, outdir) ## Don't run

mci.files <- list.files(outdir, full.names = T) ## Don't run
# save(list = 'mci.files', file = 'RData/mci.files.RData')

tmp <- str_extract(mci.files, '[0-9][0-9][0-9]_S_[0-9][0-9][0-9][0-9]')
dup <- which(duplicated(tmp))
file.remove(mci.files[dup])  ## Don't run

mci.files <- list.files(outdir, full.names = F) ## Don't run
# save(list = 'mci.files', file = 'RData/mci.filename.RData')
mci.id <- str_extract(mci.files, '[0-9][0-9][0-9]_S_[0-9][0-9][0-9][0-9]')

long.dat4 = long.dat3[which(long.dat3$ID %in% mci.id), ]
long.dat4$ID = rep(1:length(mci.id), table(long.dat4$ID))
surv.dat2 = surv.dat[which(surv.dat$ID %in% mci.id), ]
surv.dat2$ID = 1:length(mci.id)

## 
l <- NULL
for (j in 1:5){
  l <- c(l, sum(is.na(long.dat4[, j+2])))
}

l

## Missing data: 29 ADAS-Cog-13, 15 RAVLT-immediate, 15 RAVLT-learning, 8 MMSE, 35 CDR-SB
## 742 patients, 3323 visits ##
max.visit.time <- NULL
for (i in 1:length(mci.id)){
  tmp.dat <- long.dat4[which(long.dat4$ID==i), ]
  max.visit.time <- c(max.visit.time, max(tmp.dat$Time))
}

mean(max.visit.time)  ## 2.70 Years
sd(max.visit.time) ## 2.47 Years


write.csv(long.dat4, 'dataset/ADNI_long.csv', row.names = F)
write.csv(surv.dat2, 'dataset/ADNI_surv.csv', row.names = F)
