### The Environmental Neuroactive Chemicals List of Prioritized Substances for Human Biomonitoring and Neurotoxicity Testing: A Database and High-Throughput Toxicokinetics Approach 

#### Experiment 2 Goal
First, this analysis set out to investigate the fraction of a parent chemical that reaches the urine after a 5 day acute exposure through daily oral ingestion.  Second, we derived a decision tree based on important toxicokinetic parameters to inform which parameters are important in determining how much of the parent compound reaches the urine. Additionally, we categorized the chemicals for which high-throughput toxicokinetic analyses were run into groups based on the fraction of the parent compound excreted.


#### Experiment 2 Summary

The resulting list of 1,827 chemicals identified in **Experiment_1** were read into R. Then, the package httk was leveraged to predict the amount of the parent compound in the kidney tubules, a surrogate for the amount excreted in urine, over a 5-day period assuming an oral daily dose of 1mg/kgBW/day for chemicals with appropriate data.  In addition to the chemicals with measured toxicokinetic parameters, computationally derived parameters from Dawson et al. 2021 were utilized to expand the chemical coverage. Also, the hepatic clearance rate for 17 chemicals was manually adjusted to match Wetmore et al. 2015. 

Overall, 1,207 of the 1,827 chemicals had measured or computationally predicted toxicokinetic parameters that allowed us to predict a fraction excreted. A pbtk model was run for each of these chemicals and the difference in the amount of the parent compound in the tubules from the beginning to the end of day 5 of exposure, in micromoles, was extracted. This value was converted in the amount excreted in mgs/day using the chemical’s molecular weight, then to the fraction of the dose excreted assuming a bodyweight (BW) of 70kg.

After a fraction excreted value was calculated for each of the 1,207 chemicals, the r package rpart was used to generate a decision tree of the toxicokinetic parameters used by httk to inform the fraction excreted. Specifically, the parameters considered were the hepatic clearance rate (Clint), the fraction of the chemical absorbed by the gut (Fgutabs), the amount in venous blood (Aven), the ratio of the concentration of chemical in the blood to the concentration of the chemical in the plasma (Rblood2plasma), and the chemical fraction unbound in the presence of plasma proteins (Fup). After the decision tree was generated, a loop was written to walk through the tree and assign each of the 1,207 chemicals to one of the 6 terminal leaves, which reflect varying levels of fraction excreted.


#### Experiment 2 Results

The high-throughput toxicokinetic modeling is run in the script **httk_analysis.R**. The file **updated_chem_invivo_phys_props_021924.csv** contains the full chemical and in vivo  parameters used in the analysis after including the Dawson predictions and modifying the hepatic clearance rate for the aforementioned 17 select chemicals . **httk_urinary_excretion_results_021924.csv** contains the resulting fraction excreted and all parameters and intermediate values using in relevant calculations. 
The decision tree is produced in the script **decision_tree.R**. This script utilizes the results file **httk_urinary_excretion_results_021924.csv** from the httk script. The file **chemical_decision_tree_results_022124.csv** contains the “bin”, or terminal leaf of the decision tree, the chemical falls under.


