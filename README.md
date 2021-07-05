# dia-benchmarking

Processing steps:

1) Generation of bootstrap datasets (1_bootstrapping.R)

2) Generation of result tables containing protein name, p-value, log2(fold-change) for each benchmark setting combination applied to each bootstrap dataset (+ some data characteristics) (2-pValFCTable_generation.R)

3a) Partial AUC (pAUC) calculation based on p-values in result tables (3a_pAUC_calculation.R)

3b) Root Mean Square Error (RMSE) calculation based on log2(fold-change) in result tables (3b_RMSE_calculation.R)

3c) Calculation of additional data characteristics of bootstrap datasets (3c_additionalDataCharacteristics.R)

4) Visualizations (4_visualizations.R)
