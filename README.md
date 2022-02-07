# dia-benchmarking

## File information:

### Descriptive characterisation of the 17 benchmark datasets:
For the descriptive analysis of the original DIA datasets 
diaWorkflowResults_allDilutions.rds was used, which included protein intensity information to all four spike-in conditions (human only, 1:25, 1:12, 1:6)

### Reference protein lists:
'Intersection' intersectProteinNames.rds
'Combined' combinedProteinNames.rds
'DiaWorkflow' was derived directly from sublist of respective dialects workflow in diaWorkflowResults.rds (which included only information to the spike-in conditions 1:25 and 1:12 and only includes proteins, for which intensities were measured in these two conditions) by extracting the row names. 

### R scripts in the order usage: 
bootstrapping.R --> Cleaning and unification of the data matrices of the investigated 17 DIA workflows, and generation of list of 2100 vectors of column indices, which were used to generate bootstrap datasets from the data matrices in diaWorkflowResults.rds. 

benchmark_analysis.R --> Running distinct sparsity reduction - normalization - statistical test
combinations on the 2100 bootstrap datasets accessed via the list of column indices. Different prediction performance measures (partial AUC, RMSE) based on the three reference protein lists and data characteristics are calculated.
	For each of the 1428 analysis combinations generation of:
	- .csv file ("benchmark_results....csv") with 2100 rows representing the results of the 2100 	bootstrap datasets \
	- .RData file ("resultlist_....Rdata") with the information in the .csv plus a table with the output of the statistical tests and the detected fold-change for each protein 

visualization.R --> Descriptive characterisation of the original 17 DIA benchmark datasets as well as analysis of the data generated via benchmark_analysis.R.

