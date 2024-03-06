rm(list=ls())

#load libraries
library(tidyverse)
library(httk)
library(readxl)

#set working directory
setwd("Experiment_2")

#read in prioritized chemicals
prior <- read_xlsx("input/ranked_chemicals_092323.xlsx")


#get lists of chemicals from httk and Dawson et al
dawson <- load_dawson2021()
httk_chems <- get_cheminfo(info="DTXSID", model = "pbtk")

#get list of overlapping chemicals to evaluate
pbtk_chems <- httk_chems
overlap <- intersect(prior$DTXSID, pbtk_chems)
overlap_df <- data.frame(DTXSID=overlap)


#modify chem's Clint and corresponding pval for 10 micromolar from wetmore et al 2015 Table S2 for chem.physical_and_invitro.data 
chem_in <- chem.physical_and_invitro.data
chem_in_filt <- chem_in %>%
  mutate(idx=1:nrow(chem_in)) %>%
  filter(DTXSID %in% overlap_df$DTXSID) %>%
  filter(Human.Clint==0) %>%
  filter(Human.Clint.pValue>0.05) %>%
  filter(Human.Clint.Reference=="Wetmore 2015")


wetmore_data <- data.frame(DTXSID=c("DTXSID3024368","DTXSID1021956","DTXSID5024059","DTXSID0020529","DTXSID3021986","DTXSID6027052","DTXSID4020537","DTXSID5021097","DTXSID9043938","DTXSID4032459","DTXSID8022828","DTXSID0020523","DTXSID1022053","DTXSID6020561","DTXSID3020207","DTXSID5023877","DTXSID7020348"),
                           Human.Clint=pmax(0,c(-0.8196, 1.0312, 17.124, 5.704, 12.014, 5.932, -1.8806, 5.496, -0.10688, 20.66, 1.9992, 1.237, 2.936, 4.098, -0.343, 5.07, 15.856)),
                           Human.Clint.Pval=c(0.188, 0.6577, 0.0125, 0.0881, 0.0001, 0.0168, 0.0326, 0.0451, 0.6209, 0.0001, 0.2528, 0.1503, 0.1183, 0.1861, 0.2274, 0.2835, 0.0001))

dtxsids <- chem.physical_and_invitro.data %>% filter(DTXSID %in% wetmore_data$DTXSID) %>% dplyr::select(DTXSID, CAS)

wetmore_data <- merge(wetmore_data, dtxsids, by="DTXSID", all.x = TRUE)

chem.physical_and_invitro.data <- add_chemtable(new.table = wetmore_data,
              current.table = chem.physical_and_invitro.data,
              data.list = list(DTXSID="DTXSID", CAS="CAS", Clint="Human.Clint", Clint.pValue="Human.Clint.Pval"),
              species = "Human",
              reference = "Wetmore 2015",
              overwrite = TRUE)

ref_tab <- chem.physical_and_invitro.data
ref_tab <- ref_tab %>% filter(DTXSID %in% overlap)

#export updated chemical and in vivo physical properties used in analysis
write_csv(chem_in, "output/updated_chem_invivo_phys_props_021924.csv")

  

#vectors to store results
dtxsid_list <- c() #store chemical DTXSID
chem_name_list <- c() #store chemical names
day_list <- c() #store the number of days for which the model simulates
mw_list <- c() #store chemical molecular weight in g/mol
expo_list <- c() #store predicted daily exposure
frac_ex_list <- c() #store the calculated fraction of the parent compound excreted
mg_p_kg_p_day_list <- c() #store the amount of the parent compound excreted per day
fgutabs_list <- c() #store the fraction of the chemical absorbed in the gut
clint_list <- c() #store the hepatic clearance rate
aven_list <- c() #store the amount in venous blood
rblood2plasma_list <- c() # store the ratio of the concentration of chemical in the blood to the concentration of the chemical in the plasma
fup_list <- c() #store the chemical fraction unbound in the presence of plasma proteins


#establish counter for storing results
i <- 1

#make dose matrix indicating 1 mg/kgBW/day for 5 day exposure window
dm <- cbind(seq(0,5),1)
colnames(dm) <- c("time","dose")

#iterate over each chemical by DTXSID in the list of overlapping chemicals
for(c in overlap){
  print(paste0(i,"  ",c))
  
  #store DTXSID and chemical name of substance being evaluated in respective lists
  chem <- c
  name <- prior %>% filter(DTXSID==chem) %>% pull(PREFERRED_NAME)
  dtxsid_list[i] <- chem
  chem_name_list[i] <- name
  params <- parameterize_pbtk(dtxsid=chem, clint.pvalue.threshold = 0.2)
  

  # acute 5 day exposure
  day_list[i] <- "5"
  
  #run the solve_pbtk function for the chemical for 5 days at 1 dose per day of 1 mg/kg BW. 
  mod <- solve_model(model="pbtk", dtxsid = chem, dosing=list(dosing.matrix=dm), parameterize.arg.list=list(clint.pvalue.threshold=0.2))
  mod_df <- as.data.frame(mod)
  
  #isolate results for the 5th day of exposure and plot
  day_5 <- mod_df %>% filter(time>=4 & time<=5)
  ggplot(data=day_5, aes(x=time, y=Atubules))+geom_point(size=.5)
  
  #calculate the amount excreted on day 5 in micromoles by subtracting tge amount in the tubules at the beginning of the 5th day from the amount in the tubules on the end of the 5th day
  amount_exec <- day_5[nrow(day_5),"Atubules"]- day_5[1,"Atubules"]
  
  #store the molecular weight
  mw <- params$MW
  mw_list[i] <- mw

  #convert micromoles/day to  mgs/day excreted
  mgs_excreted <- amount_exec*(10^-6)*mw*(10^3)
  
  #divide by dose on last day (1 mg/kgBW/day) with assumed BW of 70kg
  frac_excreted <- mgs_excreted/70
  frac_ex_list[i] <- frac_excreted
  
  #get the estimated daily exposure 
  expo <- as.numeric(prior %>% filter(DTXSID==chem) %>% pull(ExpoCastPredictions_Total))
  expo_list[i] <- expo

  #calculate the total amount excreted and store in list
  final_mpkpd <- frac_excreted*expo
  mg_p_kg_p_day_list[i] <- final_mpkpd

  
  #extract pbtk parameter values for decision tree
  fgutabs_list[i] <- params$Fabsgut
  clint_list[i] <- params$Clint
  aven_list[i] <- (tail(day_5$Cven, n=1))*(params$Vvenc)
  rblood2plasma_list[i] <- params$Rblood2plasma
  fup_list[i] <- params$Funbound.plasma


  #adjust counter
  i <- i+1
  
}
 Sys.time()
 
#compile results into dataframe
results <- data.frame(DTXSID=dtxsid_list,
                      name=chem_name_list,
                      Fgutabs=fgutabs_list,
                      Clint=clint_list,
                      Aven=aven_list,
                      Rblood2plasma=rblood2plasma_list,
                      Fup=fup_list,
                      days_evaluated=day_list,
                      MW=mw_list,
                      frac_excreted=frac_ex_list,
                      daily_expo=expo_list,
                      mg_p_kg_p_day_excreted=mg_p_kg_p_day_list)

#export results
write_csv(results, "output/httk_urinary_excretion_results_021924.csv")




