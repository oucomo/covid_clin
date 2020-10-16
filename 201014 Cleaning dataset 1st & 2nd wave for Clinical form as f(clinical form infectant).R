################################################################################
# this script aims at cleaning dataset from 1st & 2nd wave (Da Nang) to        #
# analyze infectors-infectee pairs to see whether there                        #
# is a relationship between the clinical form in the infector and the clinical #
# form in the infectee                                                         #
# Written by Thomas Kesteman on 13 Oct 2020                                    #
################################################################################
library(xlsx)
###### loading datasets ####
llr_wave1 <- read.xlsx("C:/Users/thoma/Documents/3.7 COVID-19/Modelling - Marc/raw data/Pos COVID-19 270 update.xlsx",1,encoding = "UTF-8") #sheetName="Database"
llr_wave2 <- readRDS("C:/Users/thoma/Documents/3.7 COVID-19/Modelling - Marc/clean data/Da Nang/danang_update.rds")

###### libraries & functions ####
library(descr)
options(descr.plot = F)
`%nin%` <- function (x, table) match(x, table, nomatch = 0L) == 0L
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggExtra)
library(lmtest)
library(forcats)
library(dplyr)
library(car)
library(gtools)
library(gee)

### function to get OR & 95%CI from glm models
ORtable_glm = function(modl,ndec=2,small=T) {
  res_or = as.data.frame(format(round(exp(modl$coef),ndec),nsmall=ndec))
  res_or[,2] =  format(round(exp(confint(modl))[,1],ndec),nsmall=ndec)
  res_or[,3] =  format(round(exp(confint(modl))[,2],ndec),nsmall=ndec)
  res_or[,4] =  coef(summary(modl))[,4]
  res_or[,5] =  cut(res_or[,4],breaks=c(0,0.001,0.01,0.05,1),labels=c("***","**","*",""),include.lowest=T)     # <=0.001 *** ; (0.001-0.01] ** ; (0.01-0.05]  *  ; >0.05  " "
  colnames(res_or) = c("OR","Inf 95%CI","Sup 95%CI","p_value","Signif")
  res_or$OR_CI95 = paste0(res_or[,1]," [",res_or[,2],"-",res_or[,3],"]")
  res_or[,4] = as.character(ifelse(res_or[,4]<0.001,"<0.001",round(res_or[,4],3)))
  res_or = res_or[-1,]
  res_or$categ=rownames(res_or)
  if(small) { res_or = res_or[,c("categ","OR_CI95","p_value","Signif")]}
  return(res_or)
}

# same for GEE #
geeORtable = function(fit,ref_txt = "1.00", sepCI=" - ",interact=F) { 
  fit_res = as.data.frame(paste0(format(round(exp(fit$coef[-1]),2), nsmall=2)," [",
                                 sub(" ","",format(round(exp(coef(fit)-summary(fit)$coef[,4]*qnorm(0.975)),2), nsmall=2)[-1]),sepCI,
                                 sub(" ","",format(round(exp(coef(fit)+summary(fit)$coef[,4]*qnorm(0.975)),2), nsmall=2)[-1]),"]"))
  colnames(fit_res) = "ORCI"
  fit_res = cbind(fit_res,pval =  as.character(ifelse(2*pnorm(-abs(summary(fit)$coef[-1,5]))<0.001,"<0.001",round(2*pnorm(-abs(summary(fit)$coef[-1,5])),3))))
  fit_res = cbind(fit_res,Sign = cut(2*pnorm(-abs(summary(fit)$coef[-1,5])),breaks=c(0,0.001,0.01,0.05,1),labels=c("***","**","*",""),include.lowest=T))
  fit_res = cbind(Categ =labels(fit$coef[-1]), fit_res)
  if (!interact) {
    varfit = gsub(" ","",unlist(strsplit(substr(as.character(fit$call)[2],regexpr("~",as.character(fit$call)[2])[1]+2,nchar(as.character(fit$call)[2])),"[+]")))
    labels_varfit = lapply(varfit,FUN=function(x) {eval(parse(text=paste0("levels(",as.character(fit$call)[4],"$",x,")")))})
    varfit_rep = rep(varfit,unlist(lapply(labels_varfit,FUN=function(x){length(x)})))
    Categ = paste0(varfit_rep,unlist(labels_varfit))
    order_var= 	rep(c(1:length(varfit)),unlist(lapply(labels_varfit,FUN=function(x){length(x)})))
    tabnames_varfit = cbind(as.data.frame(Categ),as.data.frame(order_var))
    fit_res = merge(tabnames_varfit, fit_res,all=T)
    colnames(fit_res) = c("Categ.","order_var","Adjusted OR [95%IC]","pval","Signif.")
    fit_res[,5] = ifelse(is.na(fit_res[,5]),"",as.character(fit_res[,5]))
    fit_res[,4] = ifelse(is.na(fit_res[,4]),"",as.character(fit_res[,4]))
    fit_res[,3] = ifelse(is.na(fit_res[,3]),ref_txt,as.character(fit_res[,3]))
    fit_res  = fit_res[order(fit_res$order_var, fit_res$Categ.),]
    fit_res  = fit_res[,-2]
  }
  return(fit_res) 
}

# graphic parameters
orangecb = rgb(230,159,0,maxColorValue=255)
bluecb = rgb(0,114,178,maxColorValue=255)
yellcb = rgb(240,228,65,maxColorValue=255)
greencb = rgb(0,158,115,maxColorValue=255)
redcb = rgb(213,94,0,maxColorValue=255)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### 1st wave #####
llr = llr_wave1
llr <- llr[!is.na(llr$Patients.number),]
# summary(llr)

######## creation of outcome variables #########

# create variable Symptomatic (the individual has developed a clinical form, vs asymptomatic)
# freq(llr$Onset)
llr$Symptomatic = ifelse(llr$Onset %in% "CÃ³" | llr$Onset %in% "Có",T,
                         ifelse(llr$Onset %in% "KhÃ´ng" | llr$Onset %in% "Không",F,NA))
# freq(llr$Symptomatic)
# with(llr,table(Onset,Symptomatic,useNA = "ifany"))

# create variable Severe (the individual has developed a severe clinical form, vs not severe or asymptomatic)
# freq(llr$Severe.1)
llr$Severe = llr$Severe.1 %in% 1
# freq(llr$Severe.2)
# any(llr$Severe.2 %in% 2 & !llr$Severe)
# any(llr$Severe.3 %in% 3 & !llr$Severe)
# with(llr,table(Severe,Severe.1,useNA = "ifany"))
# freq(llr$Severe)
# freq(llr$Severe[llr$Symptomatic])

llr$Clinical_form = ifelse(llr$Symptomatic,"Symptomatic","Asymptomatic")
llr$Clinical_form = ifelse(llr$Severe,"Severe",llr$Clinical_form)
# freq(llr$Clinical_form)

######## clean contact #########

llr$id = paste0("NB",as.numeric(substr(llr$Patients.number,3,6)))

# freq(llr$Code.of.contact)
llr$contact = gsub("BN","NB",llr$Code.of.contact)
llr$contact = gsub(", ",",",llr$contact)
llr$contact = gsub("\n",",",llr$contact)
# freq(llr$contact)

# adding Tan's information in the dataset
llr[llr$Patients.number=="NB126","contact"] <- "NB125"
llr[llr$Patients.number=="NB158","contact"] <- "NB124"

llr$contact_id = llr$contact

llr$n_contacts = str_count(pattern = ",",as.character(llr$contact))+1
llr$n_contacts = ifelse(is.na(llr$contact),0,llr$n_contacts)
# freq(llr$n_contacts)
# freq(llr$n_contacts>0)
# freq(llr$n_contacts[llr$n_contacts>0])

#### clean date contact ##
# summary(llr$Last.contact)
# head(format(as.Date(llr$Last.contact),"%Y"))
llr$last_contact = as.Date(ifelse(format(as.Date(llr$Last.contact),"%Y") %nin% "2020",
                          paste0("2020-",format(as.Date(llr$Last.contact),"%m-%d")),
                          as.character(llr$Last.contact)),"%Y-%m-%d")
# summary(llr$last_contact)
# any(!is.na(llr$contact) & is.na(llr$last_contact))
# summary(llr$Date.of.Quarantine.statted)
llr$last_contact_noNA = as.Date(ifelse(!is.na(llr$contact) & is.na(llr$last_contact),
                                       as.character(llr$Date.of.Quarantine.statted-1),
                                       as.character(llr$last_contact)),"%Y-%m-%d")
# summary(llr$last_contact_noNA)
# any(!is.na(llr$contact) & is.na(llr$last_contact_noNA))
# summary(llr$date.of.admision)
llr$last_contact_noNA = as.Date(ifelse(!is.na(llr$contact) & is.na(llr$last_contact_noNA),
                                       as.character(llr$date.of.admision-1),
                                       as.character(llr$last_contact_noNA)),"%Y-%m-%d")
# summary(llr$last_contact_noNA)
# any(!is.na(llr$contact) & is.na(llr$last_contact_noNA))
# summary(llr$Date.of.upload)
llr$last_contact_noNA = as.Date(ifelse(!is.na(llr$contact) & is.na(llr$last_contact_noNA),
                                       as.character(llr$Date.of.upload-1),
                                       as.character(llr$last_contact_noNA)),"%Y-%m-%d")
# summary(llr$last_contact_noNA)
# any(!is.na(llr$contact) & is.na(llr$last_contact_noNA))
# date onset
# summary(llr$Date.onset)


######## creation of other covariables #########
llr$Gender_FM = ifelse(llr$Gender %in% "Nam","M","F")

# freq(llr$Age.group)
llr$Age_group = substr(llr$Age.group,1,1)
llr$Age_group = ifelse(llr$Age_group=='1','2',llr$Age_group)
# freq(llr$Age_group)
# llr$Age_group = relevel(factor(llr$Age_group),4)
# summary(llr$Age_group)
# with(llr,table(Age.group,Age_group))
llr$Age_group2 = substr(llr$Age.group,1,1)
llr$Age_group2 = ifelse(llr$Age_group=='1','3',llr$Age_group2)
llr$Age_group2 = ifelse(llr$Age_group=='2','3',llr$Age_group2)
# freq(llr$Age_group2)
# llr$Age_group2 = relevel(factor(llr$Age_group2),3)
# summary(llr$Age_group)
# with(llr,table(Age.group,Age_group2))

### cleaning comorbidities ##
# freq(llr$underlineing.medical.condition)
llr$Diabetes = grepl("\u0110T\u0110",llr$underlineing.medical.condition) | #ĐTĐ Type II = diabetes type 2
  grepl("\u0110ái tháo",llr$underlineing.medical.condition)  #Đái tháo = DM  
# freq(llr$Diabetes)
llr$Cancer = grepl("Ung th\u01B0",llr$underlineing.medical.condition)  #Ung thư  = cancer
# freq(llr$Cancer)
llr$Obesity = grepl("Béo phì",llr$underlineing.medical.condition)  # Béo phì độ = obesity
# freq(llr$Obesity)
llr$Cardiac_injury = grepl("Suy tim",llr$underlineing.medical.condition) # Suy tim = Heart failure
# freq(llr$Cardiac_injury)
llr$Hypertension = grepl("THA",llr$underlineing.medical.condition) # THA = Arterial HyperTension
# freq(llr$Hypertension)
# with(llr,table(underlineing.medical.condition,Hypertension))

llr$Comorbidity = apply(llr[,c("Diabetes","Cancer","Obesity","Cardiac_injury","Hypertension")],1,any)
# with(llr,table(underlineing.medical.condition,Comorbidity))
# freq(llr$Comorbidity)

## health care workers and hospital clusters ##
#HCW
llr$HCW = grepl("Bác s\u0129",llr$Occupation) #Bác sĩ = Doctor
llr$HCW = grepl("\u0110i\u1EC1u d\u01B0\u1EE1ng",llr$Occupation) | #Điều dưỡng  = nursing 
          grepl("\u0111i\u1EC1u d\u01B0\u1EE1ng",llr$Occupation) | llr$HCW #điều dưỡng = nursing
llr$HCW = grepl("Ng\u01B0\u1EDDi ch\u0103m",llr$Occupation) | llr$HCW #Người chăm = caregiver
# with(llr,table(Occupation,HCW,useNA = "ifany"))

#hospitalized before onset of infection
# llr$id[llr$date.of.admision < llr$date.of.first.sample & (llr$date.of.admision < llr$Date.onset)]
# identifies people placed in quarantine in hospitals after contact with case, but not hospital-acquired SARS-CoV-2 infections

#Hospital cluster
llr$Hospital_cluster = grepl("Ng\u01B0\u1EDDi b\u1EC7nh",llr$Occupation) #Người bệnh = patient
llr$Hospital_cluster = grepl("Ng\u01B0\u1EDDi \u0111i kh\u00E1m",llr$Occupation) | llr$Hospital_cluster #Người đi khám  = consulting patient
llr$Hospital_cluster = grepl("Nhân viên Y t\u1EBF",llr$kind.of.contact) | llr$Hospital_cluster #Nhân viên Y tế = medical staff
llr$Hospital_cluster = grepl("Bác s\u0129",llr$relation.with.target.contact) | llr$Hospital_cluster #Bác sĩ = Doctor
llr$Hospital_cluster = grepl("Ch\u0103m ng\u01B0\u1EDDi b\u1EC7nh",llr$relation.with.target.contact) | llr$Hospital_cluster #Chăm người bệnh = caregiver
# with(llr,table(relation.with.target.contact,Hospital_cluster,useNA = "ifany"))


llr$wave = 1
llr_wave1 = llr

########### 2nd wave ######
llr = llr_wave2
# summary(llr)

######## creation of outcome variables #########

# create variable Symptomatic (the individual has developed a clinical form, vs asymptomatic)
llr$Symptomatic = !is.na(llr$date_onset) | !is.na(llr$symptom) | llr$symptom_admission %in% "severe" | llr$symptom_admission %in% "mild" | llr$dead
# with(llr,table(Symptomatic,symptom_admission,useNA = "ifany"))
# freq(llr$Symptomatic)
# freq(llr$Symptomatic[llr$n_contacts>0])

# # create variable Severe (the individual has developed a severe clinical form, vs not severe or asymptomatic)
llr$Severe = llr$symptom_admission %in% "severe" | llr$dead
# freq(llr$Severe)
# with(llr,table(Severe,symptom_admission,useNA = "ifany"))
# freq(llr$Severe)
# freq(llr$Severe[llr$n_contacts>0])
# freq(llr$Severe[llr$Symptomatic])

llr$Clinical_form = ifelse(llr$Symptomatic,"Symptomatic","Asymptomatic")
llr$Clinical_form = ifelse(llr$Severe,"Severe",llr$Clinical_form)
# freq(llr$Clinical_form)

######## clean contact #########

llr$contact_id <- ifelse(llr$infector %in% "",NA,llr$infector)

#to be cleaned:
# "BN833)"
# "BN862,"
# "BN897,BN953,BN954,BN955"
# "BN680,BN736" 
# "BN998,BN1006" etc
# "BN958. BN730"
# "BN686. BN780" etc
# "BN 825"

if (is.list(llr$contact_id)) {
  llr$contact_id = gsub(' \\| ', ",",llr$contact_id)
  llr$contact_id = gsub(". BN", ",BN",llr$contact_id)
  llr$contact_id = gsub("c\\(\"BN", "BN",llr$contact_id)
  llr$contact_id = gsub("\"\\)", "",llr$contact_id)
  llr$contact_id = gsub("\", \"", ",",llr$contact_id)
  llr$contact_id = str_replace_all(string=llr$contact_id, pattern="BN 825", repl="BN825")
  llr$contact_id = gsub("\\)", "",llr$contact_id)
  # any( grepl("\\.",llr$contact_id))
} else {
  llr$contact_id = str_replace_all(string=llr$contact_id, pattern=",BN", repl=" | BN")
  llr$contact_id = str_replace_all(string=llr$contact_id, pattern="BN 825", repl="BN825")
  llr$contact_id = gsub(". BN", " | BN",llr$contact_id)
  llr$contact_id = gsub(")", "",llr$contact_id)
  llr$contact_id = gsub("  ", " ",llr$contact_id)
  llr$contact_id = str_replace_all(string=llr$contact_id, pattern=",", repl="")
  llr$contact_id = gsub(' \\| ', ",",llr$contact_id)
}
# tail(unique(llr$contact_id),20)

llr$n_contacts = str_count(pattern = ",",as.character(llr$contact_id))+1
llr$n_contacts = ifelse(is.na(llr$contact_id),0,llr$n_contacts)
# freq(llr$n_contacts)
# freq(llr$n_contacts>0)
# freq(llr$n_contacts[llr$n_contacts>0]) 

#### clean date contact ##
# summary(llr$date_last_contact)
# unique(format(as.Date(llr$date_last_contact),"%Y"))
llr$last_contact = as.Date(as.character(llr$date_last_contact),"%Y-%m-%d")
# summary(llr$last_contact)
# any(!is.na(llr$contact_id) & is.na(llr$last_contact))
llr$last_contact_noNA = as.Date(ifelse(!is.na(llr$contact_id) & is.na(llr$last_contact),
                                       as.character(llr$date_pos_test_1-1),
                                       as.character(llr$last_contact)),"%Y-%m-%d")
# summary(llr$last_contact_noNA)
# any(!is.na(llr$contact_id) & is.na(llr$last_contact_noNA))

# ######## creation of other covariables #########
# freq(llr$gender)
llr$Gender_FM = ifelse(llr$gender %in% "male","M","F")
llr$Gender = llr$gender
llr$Age = llr$age

llr$Age_group = cut(llr$age,breaks = c(0,16,26,41,60,71,200),right = F)
llr$Age_group = as.numeric(llr$Age_group)+1
# llr$Age_group = relevel(factor(llr$Age_group),4)
# summary(llr$Age_group)
llr$Age_group2 = as.character(llr$Age_group)
llr$Age_group2 = ifelse(llr$Age_group=='2','3',llr$Age_group2)
# llr$Age_group2 = relevel(factor(llr$Age_group2),3)
# summary(llr$Age_group)
# with(llr,table(Age_group,Age_group2))


# ### cleaning comorbidities ##
llr$Diabetes = grepl("diabetes",as.character(llr$comorbidity_name)) 
# freq(llr$Diabetes)
llr$Cancer = grepl("cancer",llr$comorbidity_name) | 
  grepl("metastasis",llr$comorbidity_name) | 
  grepl("myeloma",llr$comorbidity_name) 
# freq(llr$Cancer)
llr$Cardiovascular = grepl("myocardial ischemia",llr$comorbidity_name) |
  grepl("heart failure",llr$comorbidity_name) |
  grepl("heart disease",llr$comorbidity_name) |
  grepl("coronary artery disease",llr$comorbidity_name) |
  grepl("myocardial infarction",llr$comorbidity_name) |
  grepl("cardiovascular disease",llr$comorbidity_name) |
  grepl("cerebrovascular accident",llr$comorbidity_name) |
  grepl("stroke",llr$comorbidity_name) 
# freq(llr$Cardiovascular)
llr$Kidney_failure = grepl("kidney failure",llr$comorbidity_name) | 
  grepl("chronic kidney disease",llr$comorbidity_name) 
# freq(llr$Kidney_failure)
llr$Hypertension = grepl("hypertension",llr$comorbidity_name) 
# freq(llr$Hypertension)

llr$Comorbidity = apply(llr[,c("Diabetes","Cancer","Cardiovascular","Hypertension","Kidney_failure")],1,any)
# with(llr,table(as.character(comorbidity_name),Comorbidity))
# freq(llr$Comorbidity)

## health care workers and hospital clusters ##
#HCW
llr$HCW = grepl("doctor",llr$profession)
llr$HCW = grepl("healthcare staff",llr$profession) | llr$HCW 
llr$HCW = grepl("nurse",llr$profession) | llr$HCW 
llr$HCW = grepl("paramedic",llr$profession) | llr$HCW 
llr$HCW = grepl("pharmacist",llr$profession) | llr$HCW 
with(llr,table(profession,HCW,useNA = "ifany"))

#Hospital cluster suspected_source
llr$Hospital_cluster = grepl("hospital",llr$sus_infect_location) 
llr$Hospital_cluster = grepl("medical center",llr$sus_infect_location) | llr$Hospital_cluster 
llr$Hospital_cluster = grepl("hospital",llr$suspected_source) | llr$Hospital_cluster 
llr$Hospital_cluster = grepl("medical center",llr$suspected_source) | llr$Hospital_cluster 
with(llr,table(suspected_source,Hospital_cluster,useNA = "ifany"))

### end wave 1

llr$wave = 2
llr_wave2 = llr

########### merge wave 1 & wave 2 ######

llr_all = smartbind(as.data.frame(llr_wave2),
                    as.data.frame(llr_wave1))
dim(llr_all)
colnames(llr_all)
# table(llr_all$wave)
llr = llr_all

## identical variables wave 1/wave 2 ####

llr$contact = llr$contact_id
llr$Date_of_onset = llr$Date.onset
llr$Date_of_onset = ifelse(is.na(llr$Date.onset),llr$date_onset, llr$Date_of_onset)
# unique(llr$Date_of_onset)

# problem of pairs of infector-infectee

llr$pair_infectorinfectee_YN = F
llr$pair_infectorinfectee = NA
for (i in which(llr$n_contacts %in% 1)) {
  c = unlist(strsplit(llr[i,"contact_id"], "\\,"))
  ci = which(llr$id %in% c)
  cc = unlist(strsplit(llr[ci,"contact_id"], "\\,"))
  if (length(cc)==1 & llr$id[i] %in% cc) { # if this is a pair
    doi = llr$Date_of_onset[i]
    doc = llr$Date_of_onset[ci]
    if (!is.na(doi) & !is.na(doc) ) { # and if we can identify the individual infected first
      if (doi < doc) { # then this one is the infector
        llr$contact_id[i] <- NA # and the other one the infectee
      }
    } else { # otherwise we label the pair as such
      llr$pair_infectorinfectee_YN[i] = T
      llr$pair_infectorinfectee[i] = paste0(c,"_",cc)
      llr$pair_infectorinfectee[ci] = paste0(c,"_",cc)
    }
  }
}
# freq(llr$pair_infectorinfectee_YN)
# freq(!is.na(llr$pair_infectorinfectee))
# sum(llr$pair_infectorinfectee_YN)/2
# llr[llr$id %in% c("BN749","BN750"),c("wave","id","contact","Clinical_form","date_onset","last_contact_noNA","date_sampling","pair_infectorinfectee_YN")]
# llr[llr$pair_infectorinfectee_YN,c("wave","id","contact_id","Clinical_form","date_onset","last_contact_noNA","date_sampling"
#                                    #,"Date.onset","date.of.admision"
#                                    )]
# I don't know how to deal with that. I can't think of another rule to sort out the 
# "real" infector 

# recount N contacts
llr$n_contacts_before_clean = llr$n_contacts
llr$n_contacts = str_count(pattern = ",",as.character(llr$contact_id))+1
llr$n_contacts = ifelse(is.na(llr$contact_id),0,llr$n_contacts)
# with(llr,table(n_contacts,n_contacts_before_clean))

########## contact with symptomatic ###################

llr$contact_with_symptomatic = NA
for (i in which(llr$n_contacts>0)) {
  c = unlist(strsplit(as.character(llr[i,"contact_id"]), "\\,"))
  ci = which(llr$id %in% c)
  llr$contact_with_symptomatic[i] <- any(llr$Symptomatic[ci])
}
# freq(llr$contact_with_symptomatic)
# with(llr,table(wave,contact_with_symptomatic))

llr$contact_with_sympt_attimeofcontact = NA
for (i in which(llr$n_contacts>0)) {
  c = unlist(strsplit(as.character(llr[i,"contact_id"]), "\\,"))
  ci = which(llr$id %in% c)
  if (any(llr$Symptomatic[ci])) {
    cis = which((llr$id %in% c) & llr$Symptomatic & !is.na(llr$Date_of_onset))
    datec = llr$last_contact_noNA[i]
    dateo = llr$Date_of_onset[cis]
    llr$contact_with_sympt_attimeofcontact[i] <- any(dateo<=datec)
  } else {
    llr$contact_with_sympt_attimeofcontact[i] <- F
  }
}
# freq(llr$contact_with_sympt_attimeofcontact)
# with(llr,table(wave,contact_with_sympt_attimeofcontact))

########## contact with severe case ##

llr$contact_with_severe = NA
for (i in which(llr$n_contacts>0)) {
  c = unlist(strsplit(as.character(llr[i,"contact_id"]), "\\,"))
  ci = which(llr$id %in% c)
  llr$contact_with_severe[i] <- any(llr$Severe[ci])
}
# freq(llr$contact_with_severe)

llr$contact_with_severe_sympt_attimeofcontact = NA
for (i in which(llr$n_contacts>0)) {
  c = unlist(strsplit(as.character(llr[i,"contact_id"]), "\\,"))
  ci = which(llr$id %in% c)
  if (any(llr$Severe[ci])) {
    cis = which((llr$id %in% c) & llr$Severe & !is.na(llr$Date_of_onset))
    datec = llr$last_contact_noNA[i]
    dateo = llr$Date_of_onset[cis]
    llr$contact_with_severe_sympt_attimeofcontact[i] <- any(dateo<=datec)
  } else {
    llr$contact_with_severe_sympt_attimeofcontact[i] <- F
  }
}
# freq(llr$contact_with_severe_sympt_attimeofcontact)

#### contact with child ##########
llr$contact_with_child = NA
for (i in which(llr$n_contacts>0)) {
  c = unlist(strsplit(llr[i,"contact_id"], "\\,"))
  ci = which(llr$id %in% c)
  llr$contact_with_child[i] <- all(llr$Age[ci]<=18)
}
# freq(llr$contact_with_child)

########## age of contact ###############

llr$age_contact = NA
for (i in which(llr$n_contacts %in% 1)) {
  c = unlist(strsplit(llr[i,"contact_id"], "\\,"))
  ci = which(llr$id %in% c)
  llr$age_contact[i] <- llr$Age[ci]
}
# freq(!is.na(llr$age_contact))
# length(unique(llr$contact[llr$n_contacts %in% 1]))
# freq(llr$age_contact)


llr$age_gp = cut(llr$Age,breaks = seq(0,100,10),right = F)
llr$age_gp = ifelse(llr$Age >=100,"[90,100)",as.character(llr$age_gp))
# freq(llr$age_gp)
llr$age_contact_gp = cut(llr$age_contact,breaks = seq(0,100,10),right = F)
# freq(llr$age_contact_gp)
for (i in which(llr$n_contacts>1)) {
  c = unlist(strsplit(llr[i,"contact_id"], "\\,"))
  ci = which(llr$id %in% c)
  ag = llr$age_gp[ci]
  if (length(unique(ag))==1) {
    llr$age_contact_gp[i] <- ag[1]
  }
}
# freq(!is.na(llr$age_contact_gp))
# freq(llr$age_contact_gp)


llr$Ageover40 = llr$Age>40
llr$age_contact_over40 = NA
for (i in which(llr$n_contacts>0)) {
  c = unlist(strsplit(as.character(llr[i,"contact_id"]), "\\,"))
  ci = which(llr$id %in% c)
  ag = llr$Ageover40[ci]
  if (length(unique(ag))==1) {
    llr$age_contact_over40[i] <- ag[1]
  }
}
# freq(!is.na(llr$age_contact_over40))
# freq(llr$age_contact_over40)

#### export environment ####
save.image(file="C:/Users/thoma/Documents/3.7 COVID-19/Modelling - Marc/Rmarkdown/201014 both datasets.RData")

