rm(list=ls())

library(janitor)
library(xlsx)
library(readxl)
library(reshape2)
library(tidyverse)


setwd("Experiment_1")


#Read in neuroactive chemicals from Kosnik et al. 2020
Kosnik <- read_csv("input/Rager_Urine_Metabolomics_Priority_Substances.csv")
Kosnik <- Kosnik %>% select(DTXSID, PREFERRED_NAME, CASRN) %>% mutate(InVitro_Neuroactive_Kosnik=TRUE)

#Read in EPA list of neurotoxic substances in vivo
invivo_neurotox <- read_csv("input/DNTREF_2022_10_31.csv")
invivo_neurotox <- invivo_neurotox %>%
  select(DTXSID, `PREFERRED NAME`, CASRN) %>%
  mutate(DTXSID=str_extract(DTXSID,"DTXSID\\d+")) %>% 
  mutate(InVivo_Neurotox_List=TRUE) %>% 
  rename(PREFERRED_NAME="PREFERRED NAME")

#Read in EPA list of lit mined neuorotoxic substances
pubmed_neurotox <- read_csv("input/LITMINEDNEURO_2022_10_31.csv")
pubmed_neurotox <- pubmed_neurotox %>%
  select(DTXSID, `PREFERRED NAME`, CASRN) %>%
  mutate(DTXSID=str_extract(DTXSID,"DTXSID\\d+")) %å>% 
  mutate(PubMedLitReview_Neurotox_List=TRUE) %>% 
  rename(PREFERRED_NAME="PREFERRED NAME")

#Read in EPA list of chemicals characterizing the blood exposome
blood_expo <- read_csv("input/BLOODEXPOSOME_2022_12_09.csv")
blood_expo <- blood_expo %>%
  select(DTXSID, `PREFERRED NAME`, CASRN) %>%
  mutate(DTXSID=str_extract(DTXSID,"DTXSID\\d+")) %>% 
  mutate(Blood_Exposome=TRUE) %>% 
  rename(PREFERRED_NAME="PREFERRED NAME")

#Read in EPA list of chemicals of emerging concern
cec <- read_csv("input/CECSCREEN_2022_10_31.csv")
cec <- cec %>%
  select(DTXSID, `PREFERRED NAME`, CASRN) %>%
  mutate(DTXSID=str_extract(DTXSID,"DTXSID\\d+")) %>% 
  mutate(Chemicals_EmergConcern_Metabolites=TRUE) %>% 
  rename(PREFERRED_NAME="PREFERRED NAME")

#Read in chemicals from CPDat selected children's products PUCs
cpdat <- read_csv("input/CPDAT_child_prod_chems.csv")
cpdat <- cpdat %>%
  mutate(CPDat_Child_Prod=TRUE) %>% 
  rename(PREFERRED_NAME="PREFERRED NAME")

#Read in list of chemicals for biomonitoring in Pellizzari et al. 2019
echo <- read_csv("input/dashboard_Pellizzari_2019_chems_DTXSIDs_120922.csv")
echo <- echo %>%
  mutate(DTXSID=str_extract(DTXSID,"DTXSID\\d+")) %>% 
  select(DTXSID, PREFERRED_NAME, CASRN) %>%
  mutate(ECHO_Prioritization=TRUE)

#Create a mapping of chemical names/DTXSID/CASRN
id_map <-rbind(invivo_neurotox %>% select(!InVivo_Neurotox_List),
               pubmed_neurotox %>% select(!PubMedLitReview_Neurotox_List),
               blood_expo %>% select(!Blood_Exposome),
               cec %>% select(!Chemicals_EmergConcern_Metabolites),
               echo %>% select(!ECHO_Prioritization))

id_map <- id_map %>% distinct()


#Outer join  all these lists together retaining as many chemicals as possible, except for CPDat which is a left join in attempt to drop uninformative
#chemicals e.g. water. Also keep track of which chemicals are in which lists.
merge1 <- merge(Kosnik %>% select(DTXSID,InVitro_Neuroactive_Kosnik),invivo_neurotox %>% select(DTXSID,InVivo_Neurotox_List), by="DTXSID", all=TRUE)
merge2 <- merge(merge1, pubmed_neurotox %>% select(DTXSID, PubMedLitReview_Neurotox_List), by="DTXSID", all = TRUE)
merge3 <- merge(merge2, blood_expo %>% select(DTXSID, Blood_Exposome), by="DTXSID", all=TRUE)
merge4 <- merge(merge3, cec %>% select(DTXSID, Chemicals_EmergConcern_Metabolites), by="DTXSID", all = TRUE)
merge5 <- merge(merge4, echo %>% select(DTXSID, ECHO_Prioritization), by="DTXSID", all=TRUE)
merge6 <- merge(merge5, cpdat %>% select(DTXSID, CPDat_Child_Prod), by="DTXSID", all.x = TRUE)


merge7 <- merge(merge6, id_map, by="DTXSID", all.x = TRUE)
merge7 <- merge7 %>%
  select(DTXSID, CASRN, PREFERRED_NAME, InVitro_Neuroactive_Kosnik, InVivo_Neurotox_List, PubMedLitReview_Neurotox_List, Blood_Exposome, Chemicals_EmergConcern_Metabolites, CPDat_Child_Prod, ECHO_Prioritization) %>%
  replace(is.na(.), "")
  


#Read in total ExpoCast Predictions which span multiple files 
files <- list.files("input/dashboard_exposure_preds")

expo_preds <- data.frame()
for(f in files){
  temp <- read_csv(paste0("input/dashboard_exposure_preds/",f))
  expo_preds <- rbind(expo_preds,temp)
  
}

#If exposure estimates did not exist, gap fill with 1.2e-6 per expert judgment (standard SEEM3 output if lacking data)
expo_preds <- expo_preds %>%
  mutate(ExpoCastPredictions_Total=`EXPOCAST_MEDIAN_EXPOSURE_PREDICTION_MG/KG-BW/DAY`) %>% 
  select(DTXSID, PREFERRED_NAME, ExpoCastPredictions_Total) %>% 
  mutate(ExpoCastPredictions_Total=as.numeric(ExpoCastPredictions_Total)) %>% 
  mutate(ExpoCastPredictions_Total=if_else(is.na(ExpoCastPredictions_Total), 1.2e-6, ExpoCastPredictions_Total))

#Read in predictions just for 3-5 year olds for whatever chemicals there exists predictions
NHANES_preds_3to5 <- read_xlsx("input/NHANES_inferences_3_to_5_ZS.xlsx", sheet="Table S4.", skip=1)
NHANES_preds_3to5 <- NHANES_preds_3to5 %>% filter(Population=="3 - 5 years") %>% mutate(Median_Exposure_Preds_3to5=Median_Exposure) %>%  select(DTXSID, Median_Exposure_Preds_3to5)


#Merge in ExpoCast predictions into  dataframe
merge8 <- merge(merge7, expo_preds %>% select(DTXSID,ExpoCastPredictions_Total), by="DTXSID", all.x = TRUE)
merge9 <- merge(merge8, NHANES_preds_3to5, by="DTXSID", all.x=TRUE)


#Read in EPA list of chemicals that are DMSO insoluble
dmso_ins <- read_csv("input/CHEMINV_dmsoinsolubles_2022_11_08.csv")
dmso_ins <- dmso_ins %>%
  select(DTXSID, `PREFERRED NAME`, CASRN) %>%
  mutate(DTXSID=str_extract(DTXSID,"DTXSID\\d+")) %>% 
  mutate(DMSO_Insolubles=TRUE) %>% 
  rename(PREFERRED_NAME="PREFERRED NAME")


#Read in EPA list of chemicals that are too volatile to test
vols <- read_csv("input/CHEMINV_volatiles_2022_11_08.csv")
vols <- vols %>%
  select(DTXSID, `PREFERRED NAME`, CASRN) %>%
  mutate(DTXSID=str_extract(DTXSID,"DTXSID\\d+")) %>% 
  mutate(Too_Volatile_toTest=TRUE) %>% 
  rename(PREFERRED_NAME="PREFERRED NAME")


#Merge DMSO and Volatility status 
merge10 <- merge(merge9, dmso_ins %>% select(DTXSID, DMSO_Insolubles), by="DTXSID", all.x=TRUE)
merge11 <- merge(merge10, vols %>% select(DTXSID, Too_Volatile_toTest), by="DTXSID", all.x = TRUE)


#Organize final dataframe 
final_df <- merge11 %>% 
  replace(is.na(.), "") %>% 
  select(DTXSID,CASRN,PREFERRED_NAME, InVitro_Neuroactive_Kosnik, InVivo_Neurotox_List, PubMedLitReview_Neurotox_List,
         Blood_Exposome, Chemicals_EmergConcern_Metabolites, CPDat_Child_Prod,ECHO_Prioritization, ExpoCastPredictions_Total,
         Median_Exposure_Preds_3to5, DMSO_Insolubles, Too_Volatile_toTest)

write_csv(final_df, "output/neuroactive_chems_no_filt_summary_092323.csv")

############################################
# Manual curation of previous iterations
############################################

#List of chemicals to retain that after previous manual curation
keep_1827 <- read_xlsx("input/1827_chems_for_prioritization_092323.xlsx")


#Apply filters to further refine chemical list
chem_filt <- final_df %>%
  filter(InVitro_Neuroactive_Kosnik==TRUE | InVivo_Neurotox_List==TRUE | PubMedLitReview_Neurotox_List==TRUE) %>% #must be in a least one neuro list
  filter(Blood_Exposome==TRUE | Chemicals_EmergConcern_Metabolites==TRUE | CPDat_Child_Prod==TRUE | ECHO_Prioritization==TRUE) %>% #must be in at least one exposure list
  filter(DMSO_Insolubles!=TRUE) %>% #must NOT be DMSO insoluble
  filter(Too_Volatile_toTest!=TRUE) %>% #must NOT be too volatile to test
  filter(DTXSID %in% keep_1827$DTXSID)

#Apply prioritization schema
prior <- chem_filt %>% 
  mutate(across(.cols=-c(DTXSID, CASRN, PREFERRED_NAME, ExpoCastPredictions_Total, Median_Exposure_Preds_3to5), function(x) replace(x,x=="", FALSE))) %>%  #for all 7 binary lists, if the value is missing make the value FALSE
  mutate(across(.cols=-c(DTXSID, CASRN, PREFERRED_NAME, ExpoCastPredictions_Total, Median_Exposure_Preds_3to5), function(x) as.integer(as.logical(x)))) %>% #Convert TRUE values to 1 and FALSE values to 0 for binary lists
  mutate(list_score=InVitro_Neuroactive_Kosnik+InVivo_Neurotox_List+PubMedLitReview_Neurotox_List+Blood_Exposome+Chemicals_EmergConcern_Metabolites+CPDat_Child_Prod+ECHO_Prioritization) %>% #add up the number of selected lists that a chemical is on to form a "list score"
  mutate(list_score_scale= (list_score-min(list_score))/(max(list_score)-min(list_score))) %>% #min/max scale the list score
  mutate(ExpoCastPredictions_Total=as.numeric(ExpoCastPredictions_Total)) %>% #make sure the ExpoCast predictions are numeric
  arrange(desc(ExpoCastPredictions_Total)) %>% #arrange the ExpoCast predictions from highest to lowest
  mutate(expocast_rank_score=nrow(chem_filt):1) %>% #assign an "expocast rank score" based on the order of ExpoCast predictions
  mutate(expocast_rank_score_scale= (expocast_rank_score-min(expocast_rank_score))/(max(expocast_rank_score)-min(expocast_rank_score))) %>% #min/max scale the expocast rank score
  mutate(overall_score=list_score_scale+expocast_rank_score_scale) %>% #add the scaled list score and the scaled expocast rank score to produce an "overall score"
  arrange(desc(overall_score)) %>% #arrange by overall score
  mutate(overall_rank=1:nrow(chem_filt)) #assign a rank


#Export ranking results
write.xlsx(prior, "output/ranked_chemicals_092323.xlsx", row.names=FALSE)

