rm(list=ls())

#load libraries
library(tidyverse)
library(rpart)
library(rpart.plot)
library(RColorBrewer)
library(xlsx)
library(readxl)


#set working directory
setwd("Experiment_2")

#read in httk results
res <- read_csv("output/httk_urinary_excretion_results_021924.csv")

#rename cols for downstrean interpretation
colnames(res)[3:12] <- lapply(colnames(res)[3:12], function(x) paste0("httk_",x))

#make df for decision tree
tree_data <- res %>% column_to_rownames("DTXSID") %>%
  dplyr::select(httk_frac_excreted, httk_mg_p_kg_p_day_excreted, httk_Fgutabs, httk_Clint, httk_Aven, httk_Rblood2plasma, httk_Fup) %>%
  mutate(across(.cols=everything(),.fns=as.numeric))

set.seed(17)

#make decision tree for fraction excreted dependening on all toxicokinetic parameters for the chemical
tree <- rpart(httk_frac_excreted~httk_Fgutabs+httk_Clint+httk_Aven+httk_Rblood2plasma+httk_Fup, data=tree_data)
rpart.plot(tree, type=1, extra = 1, digits = 4)
rpart.rules(tree)

#make a list of chemicals to loop through
chems <- unlist(res$DTXSID)

#lists to store loop results
dtxsid_list <- c() #store chemical DTXSID for downstream merging
cat_list <- c() #store which branch of decision tree chemical fall in

#set counter
i <- 1

#iterate over all chemicals and store which of the 6 final bins of the decision tree they end up in, reflective of the fraction of the parent compound excreted
for(chem in chems){
  print(paste0(i,"  ",chem))
  dtxsid_list[i] <- chem
  res_chem <- res %>% filter(DTXSID==chem)
  if(res_chem$httk_Clint>=0.82){
    if(res_chem$httk_Clint>=4.548){
      cat_list[i] <- 1
    }
    else{
      cat_list[i] <- 2
    }
  }
  else{
    if(res_chem$httk_Fup<0.05718){
      if(res_chem$httk_Aven<0.4609){
        cat_list[i] <- 3
      }
      else{
        cat_list[i] <- 4
      }
    }
    else{
      if(res_chem$httk_Fgutabs<0.8709){
        cat_list[i] <- 5
      }
      else{
        cat_list[i] <- 6
      }
    }
  }
  
  i <- i+1
}
  
  
#make a dataframe of the results 
results <- data.frame(DTXSID=dtxsid_list,
                      httk_result_category=cat_list)

#merge in additional chemical information
results_merge <- merge(res, results, by="DTXSID")

#export results
write_csv(results_merge, "output/chemical_decision_tree_results_022124.csv")


