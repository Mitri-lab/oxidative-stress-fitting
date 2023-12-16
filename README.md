# oxidative-stress-fitting
Code for Di Martino et al 2023 "Oxidative stress changes interactions between two bacterial
species from competitive to facilitative" (https://www.biorxiv.org/content/10.1101/2023.05.24.542164v1.abstract). 

This repository contains the R Markdown scripts, notebooks and C++ code used to fit mathematical models to experimental data and simulate the growth of At and Ct in different linoleic acid concentrations.

"Fitting_DiMartinoetal2023.Rmd" contains the script and an HTML version of the notebook that follows the fitting steps in 2 models. In the first model, the toxicity of the environment is implicit. Data corresponding to the first model is loaded through the "CFU_data_model1.RData file". In the second model, we include data on measured Reactive Oxygen Species dynamics. Data corresponding to the second model is loaded through the "CFUROS_data_model2.RData" file. 

Once the parameters have been obtained, we simulated serial transfers (see end of R Markdown file for an example). In order to run transfers in a broader parameter space to explore coexistence (Fig 4 in the paper), we implement a parameter sweeps in C++. The code for the transfers for model 1 can be found in "transfersmodel1.cpp" and for model 2 in "transfersmodel2.cpp". 

To compile the C++ scripts, run the following command in a terminal: "g++ transfersmodel2.cpp -o nameoftheex -lgsl -lgslcblas" then execute "./nameoftheex". An example of the output of this C++ script is given in "Model2_cocult_5.txt"
