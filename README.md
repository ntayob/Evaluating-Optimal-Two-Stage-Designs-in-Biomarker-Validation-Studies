# Evaluating-Optimal-Two-Stage-Designs-in-Biomarker-Validation-Studies
We have provided a framework to evaluate two-stage biomarker study designs and have provided R code for researchers to facilitate their more widespread use in validation studies.

These R scripts allow for recreation and manipulation of results presented in "Evaluating Optimal Two-Stage Designs in Biomarker Validation Studies". 

The code presented is based on the reference set of simulations. 

- The file "Parallel_Simulation_Submit.R" runs each simulation in parallel on a computing cluster. 

- The file "Simulation_Calculations_Summarize.R" compiles the results to calculate the values present in the manuscript.

The changeable parameters are primarily in the middle section of the "Parallel_Simulation_Submit.R" script. See below for line numbers corresponding with various parameters:

Line 210: True values of ROC(t)
Line 211: Value of t
Lines 214-216: Parameters of correlation matrix
Line 217: Definition of correlation matrix (note dimensions)
Lines 222-233: Definition of mean vector
Lines 239-245: Setup of biomarker panel
Lines 247-248: Setup of bootstrap replications
Lines 251-256: Object construction (note dimensions in lines 255 and 256)
Line 591: Null hypothesis for power calculations



