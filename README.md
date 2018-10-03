Comprehensive-analysis-of-mutational-processses-in-Primary-and-recurrent-Medulloblastoma
Introduction:
Key out comes - 1) Perform Mutational Analysis in R 2) Identify the mutational signatures in the Primary and Recurrent Medulloblastoma 3) Identify the similarity between the mutational profiles and COSMIC signatures 4) Analyze the result to detect any differences in the mutational profile of primary and resurrent medulloblastoma. 
#####
R Packages Used 
1) deConstructSigs
2) MutationalPatterns
3)ggplot2
4)BSgenome
5)Maftools
######
Base subsitutions represents records the different mutational processes that have the cell have been exposed to during its life time. DNA in the cell is under constant stress by sources of DNA damage , such as UV-light, nascent oxygen species etc. The cell also harbors repair mechnisms to balance these insults and repair the damage caused by the insults. Somatic Mutations are acquired by two broad mechanisms 1) when the extent of damage by the exogenous agent cannot be completely be repaired by canonical repair mechanisms or 2) there are defect in the canonical DNA repair mechanisms like MMR. Acquired smatic mutations can have functional consequences such as cancer. The sudy of these mutational processes could inform a great deal about evolution of a disease. In current scenario, we are keen on learning about difference in these mutaional processes between primary and recurrent medulloblastoma. 
#####
Data Source: The WGS of Primary , Recurrent medulloblastoma and germline control (in some cases) was downloaded from  ICGC database. 

Data Processing:




#####
Alexandrove and collegaues develpoed an alogrithm using non-negative matric factorization (NMF) and model selection to extract the signatures of mutational processes presnt in 5 million mutations in over 7000 cancer genomes to identify about 30 signatures presnt accross 30 tumor types. In simulations it was found that atleat 200 WGS were required to determine the signatures of 20 mutational processes. With exome sequencing covering only % of the human genome , resulting in fewer mutations identified , it would require thousands of samples. Therefore, application of there method was not feasible in this case of about 50 samples, here we used 2 establisded R-packages deconstructSigs and MuttaionalPatterns

###
Step1 :
