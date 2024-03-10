### The Environmental Neuroactive Chemicals List of Prioritized Substances for Human Biomonitoring and Neurotoxicity Testing: A Database and High-Throughput Toxicokinetics Approach 

#### Experiment 3 Goal
This analysis set out to understand possible near-field sources of exposure to the 250 identified ENRICH chemicals with a specific focus on the household environment.


#### Experiment 3 Summary
Product Use Category (PUC) data from the US EPA’s Chemicals and Products Database (CPDat) was [acquired](https://comptox.epa.gov/chemexpo). This hierarchical system categorizes consumer products and tracks chemicals used in the manufacturing of these products (Isaacs et al. 2020). 

The PUC data was filtered to show unique pairs of PUCs and the 250 ENRICH chemicals yielding a dataset providing insight into possible exposure sources for substances in this chemical set. A heatmap was then generated showing the ENRICH chemicals that were associated with at least three unique PUCs and bar graph was made that showed the number of total unique chemicals associated with each PUC.

#### Experiment 3 Results
All analysis is conducted in the script Enrich_exposures.R. This script reads in the 250 ENRICH chemicals as well as the CPDat Bulk Composition data containing chemicals and associated PUCs. Resulting figures are **ENRICH_PUC_HM_030724.png** 
<img src="https://github.com/Ragerlab/Script_for_The-Environmental-Neuroactive-Chemicals-List-of-Prioritized-Substances/blob/main/Experiment_3/figures/ENRICH_PUC_HM_030724.png" alt="drawing" width="800" height="650"/>

and **ENRICH_PUC_bar_chart_030724.png**.
<img src="https://github.com/Ragerlab/Script_for_The-Environmental-Neuroactive-Chemicals-List-of-Prioritized-Substances/blob/main/Experiment_3/figures/ENRICH_PUC_bar_chart_030724.png" alt="drawing" width="800" height="400"/>
