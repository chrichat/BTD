# BTD
Blind fMRI Source Unmixing via Higher-Order Tensor Decompositions

All the scripts presented in this repository are used in the paper "Blind fMRI Source Unmixing via Higher-Order Tensor Decompositions" submitted in Journal of Neuroscience methods.

The data used in the journal paper (as explained in the dataprofile.xml) was either generated or obtained by other researchers and we have used them (after applying some modifications). In order to acknowledge those researchers and highlight their job, we do not upload the data, but they can be obtained either via their repositories or by contacting them. 
Simulation of the perception study (Section 3.1.1)- used laso in [1] - Contact the author N. Helwig    helwig@umn.edu 
Simulation with artifacts (Section 3.1.2) - used also in [2] - Download via the link http://mlsp.umbc.edu/resources.html
Multislice simulation (Section 3.1.3) - used also in [3] - Contact the author A. Stegeman stegeman.alwin@gmail.com
Real data (Section 3.2) - used laso in [4] - Download via the link    https://openneuro.org/datasets/ds000157/versions/00001

As prequisite for the Tensorlab 3 [5] must be downloaded and included in the path of Matlab. Furthermore all the folders of the repository must be included in the path. In order to compare also with Sparce ICA used in [6] the code can be downloaded via http://mlsp.umbc.edu/SparseICA_EBM.html

In order to reproduce the figures 8,9,18,19,20 of the paper the file create_figures.m must be executed (the correlations saved in file journalbtd_correlations.m are used, those are the tables also presented in Section 3 of the Supplementary Material).

The code of Heuristic one is included in the function zstat.m . It's use L=zstat(U,p) returns the estimated L value. (Input, U is the cell arry which results of function lp_nls algorithm using an overestimation of the initial L and a P-value usually 0,05 or 0,01).

In order to obtain the correlation values used in Section 3.1.1 execute file section_311_perception.m
For the correlation values used in Section 3.1.2 execute file section_312_simtb.m
While for the correlation values used in Section 3.1.3 the function section_313_multislicesim.m can be used. 
[corr_tpica, corr_cpd, corr_btd, corr_pfac2, corr_btd2]=section_313_multislicesim(SNR, dataset,iter)  as an input are defined the SNR. the dataset ('A','C', 'G', 'I') and the number of iterations used.


For the analysis of the real data after downloading the data from openfmri and changing the niftiifile from type 2 to type 1 (the code is inlcuded also but commented) the file run_btdanalysis2_REAL  (in the Real data folder) must be executed. In order to create the augmented datasets of Sections 3.2.0.1 the file make_augemented_dataset can be used.

[1] N. E. Helwig and S. Hong, "A critique of Tensor Probabilistic Independent Component Analysis: implications and recommendations for multi-subject fMRI data analysis,"  J. Neuroscience Methods, vol. 213, no. 2, pp. 263-273, Mar. 2013.

[2]E. B. Erhardt, S. Rachakonda, E. J. Bedrick, E. A. Allen, T. Adali, and V. D. Calhoun, "Comparison of multi-subject ICA methods for analysis of fMRI data," Human Brain Mapping, vol. 32, no. 12, pp. 2075-2095, 2011.12.

[3] A. Stegeman, "Comparing Independent Component Analysis and the PARAFAC model for artificial multisubject fMRI data," Unpublished Technical Report, Univ. of Groeningen, Feb. 2007.

[4] P. A. Smeets, F. M. Kroese, C. Evers, and D. T. de Ridder, "Allured or alarmed: counteractive control responses to food temptations in the brain," Behavioural Brain Research, vol. 248, no. 1, pp. 41-45, Jul. 2013.

[5] N. Vervliet, O. Debals, L. Sorber, M. Van Barel, and L. De Lathauwer, "Tensorlab 3.0," Mar. 2016. https://www.tensorlab.net/

[6] Z. Boukouvalas, Y. Levin-Schwartz, Vince D. Calhoun, and T. Adali, "Sparsity and Independence: Balancing of two Objectives in Optimization for Source Separation with Application to fMRI Analysis," Elsevier, Journal of the Franklin Institute (JFI), Engineering and Applied Mathematics, 2017. 
