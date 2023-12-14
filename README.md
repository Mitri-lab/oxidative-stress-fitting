# oxidative-stress-fitting
Fitting code for Di Martino et al 2023 - fitting bacterial growth in different oxidative stress concentrations (https://www.biorxiv.org/content/10.1101/2023.05.24.542164v1.abstract). 

This repository contains the R Markdown script that follows the steps used to fit the growth of At and Ct in different linoleic acid concentrations: Fitting_DiMartinoetal2023.Rmd as well as the HTLM version of the notebook. 
In the first model, the toxicity is implicit. Data corresponding to the first model is loaded through the "CFU_data_model1.RData file". 
In the second model, we include information on Reactive Oxygen Species dynamics. Data corresponding to the second model is loaded through the "CFUROS_data_model2.RData" file. 

Once the parameters have been obtained, we can simulated serial transfers (see end of R Markdown file for an example). In order to run transfers in a broader parameter area to explore coexistence (Fig 4), we implement the parameter sweep in C++. The code for the transfers for model 1 can be found in "transfersmodel1.cpp" and for model 2 in "transfersmodel2.cpp". 

To compile the C++ scripts, run the following command in a terminal: "g++ transfersmodel2.cpp -o nameoftheex -lgsl -lgslcblas" then execute "./nameoftheex".
