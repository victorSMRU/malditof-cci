# Identification of Southeast Asian _Anopheles_ mosquito species with matrix-assisted laser desorption/ionization time-of-flight mass spectrometry using a modified cross-correlation approach

This repository provides data and code for the analysis of a study carried out to develop an open-access data analysis pipeline for the identification of _Anopheles_ mosquito species identification with matrix-assisted laser desorption/ionization time-of-flight mass spectrometry using a modified cross-correlation approach.

This analysis was carried out by Victor Chaumeau* who developed the analysis framework and wrote the R code for mass spectral fingerprint matching to mosquito species using a cross-correlation algorithm based on the original work of Arnold and Reilly.

_Arnold RJ, Reilly JP. 1998. Fingerprint matching of E. coli strains with matrix-assisted laser desorption/ionization time-of-flight mass spectrometry of whole cells using a modified correlation approach. Rapid Commun Mass Spectrom 12:630–636. doi:10.1002/(SICI)1097-0231(19980529)12:10<630::AID-RCM206>3.0.CO;2-0_

\*For correspondence: victor@tropmedres.ac

The preprint with the results: [to be completed]

## Data

The data are provided in the folder data.

## Code

The code used for this analysis is provided in the folder code and in the R Markdown documents _data_preparation.Rmd_, _msl_construction.Rmd_, _analysis.Rmd_, _msl_upgrade.Rmd_ and _sensitivity_analysis.Rmd_.

The generated files are not included in the GitHub repository to keep it reasonably sized. The R Markdown documents should be run in this order to generate the files and reproduce the analysis. 

## Main result

The results in the paper are generated in the R Markdown _analysis.Rmd_.

### Principle of the cross-correlation algorithm

Spectra were visualized and processed using R software version 4.2 and the MALDIquant package (Gibb and Strimmer, 2012; R Core Team, 2013). Processing of raw mass spectra included intensity transformation and smoothing, baseline removal, normalization of intensity values and spectra alignment (Figure 1). A cross-correlation index (CCI) was then calculated for each pairwise spectra comparison using a custom algorithm adapted from the original work of Arnold and Reilly (Arnold and Reilly, 1998). The algorithm outputs the maximum value of the cross-correlation function for two input mass spectra over specified mass intervals, and the CCI is given by the product of the local cross-correlation value of each mass interval. If no local maximum is detected on a given mass interval, the algorithm returns 0. Therefore, the range of possible CCI values on the log scale (log10CCI) is a real number varying between a negative infinite limit (too dissimilar spectra) and 0 (identical spectra). In this assessment, the CCI algorithm was parameterized with mass intervals of 500 Da spanning the mass range 3000-12000 Da.

![__Figure 1. Principle of the cross-correlation algorithm.__ (A) A typical raw MALDI-TOF mass spectrum of _An. minimus_, (B) the corresponding processed spectrum after intensity normalization and baseline removal, (C) comparison of the processed mass spectra of two _An. minimus_ specimens over the 5000-5500 kDa mass interval showing high similarity between spectra, (D) comparison of the processed mass spectra of an _An. minimus_ specimen and of an _An. maculatus_ specimen over the 5000-5500 kDa mass interval showing limited similarity between spectra, (E) the cross-correlation function of the two _An. minimus_ spectra over the 5000-5500 kDa mass interval gives a local maximum of 0.982 and (F) the cross-correlation function of the _An. minimus_ and _An. maculatus_ spectra over the 5000-5500 kDa mass interval gives a local maximum of 0.540. If no local maximum of the cross-correlation function is detected, the algorithm is parameterized to return 0. The resulting cross-correlation index on the log scale (log10CCI) over the 3000-12000 kDa mass range is -4.9 for the two _An. minimus_ spectra and -Infinite for the _An. minimus_ and _An. maculatus_ spectra.](analysis_files/figure-html/figure1_CCI_principle-1.png)

<br/>

### Construction of the reference mass spectra database

2535 mass spectra of the 254 Anopheles specimens selected for inclusion in the reference panel and identified with PCR (25 taxa including 20 _sensu stricto_ species and 5 sibling species pairs or complexes) were acquired, yielding 3,211,845 pairwise comparisons of distinct spectra pairs. A CCI value reflecting similarity between spectra was computed for every spectra pairs. The repeatability, reproducibility and specificity of mass spectra was assessed using subsets of the comparisons (Figure 2): the median log10CCI was -7.9 (inter-quartile range [IQR]: -9.2 to -6.8) for comparisons of technical replicates of the same specimen, -10.7 (IQR: -12.6 to -9.4) for comparisons of different specimens of the same species and -Infinite (IQR: - Infinite to -Infinite) for comparisons of different species. 

![__Figure 2. Repeatability, reproducibility and specificity of the mass spectra.__ (A) median log10CCI of pairwise comparisons between technical replicates of the same specimen collated by mass spectrum and (B) corresponding density function, (C) median log10CCI of pairwise comparisons between spectra of different specimens of the same species collated by mass spectrum and (D) corresponding density function, (E) median log10CCI of pairwise comparisons between spectra of different species collated by mass spectrum and (F) corresponding density function.](analysis_files/figure-html/figure3_spectra_reproducibility-1.png)

<br/>

### Evaluation of the performance of _Anopheles_ species identification with MALDI-TOF MS

In order to evaluate the performance of _Anopheles_ species identification with MALDI-TOF MS, 1049 mass spectra of the 105 PCR-identified specimens included in the test panel were queried against the reference database, yielding 2,659,215 pairwise comparisons. A simulation experiment was carried out by selecting at random a specified number of spots per specimen of the test panel varying from 1 to 9 with 1000 repeats. The MALDI-TOF identification result was defined as the reference spectrum giving the highest CCI value, considering all randomly selected spots in the analysis. An identification threshold was set for the CCI value below which no identification result was given. The MALDI-TOF identification result was categorized as true positive (test specimen matching with the same species in the reference mass spectra database with a CCI value above the identification threshold), true negative (specimen of a species not represented in the reference mass spectra database with a CCI value below the identification threshold), false positive (matches between different species with a CCI value above the identification threshold) and false negative (species represented in the reference mass spectra database matching with a CCI value below the identification threshold), considering the PCR identification result as the reference. The results of this simulation experiment were used to estimate the sensitivity, positive predictive value and accuracy. Noteworthy, as species coverage in the reference mass spectra database increases with the sample size, very few specimens of species not referenced in the reference mass spectra database are included in the test panel thereby impairing the estimation of specificity and negative predictive value (because there is no true negative result in the assessment). Therefore, a separate simulation experiment was carried out using the best match of each spectrum after excluding comparisons of the same species in the analysis (i.e., forcing the result to be either false positive or true negative by mimicking queries of unreferenced species), allowing estimation of the specificity for queries of unreferenced species (i.e., probability of a result being falsely positive if the query species is not referenced in the reference database) but not of the negative predictive value (because the negative predictive value also depends on the proportion of false negative). Performance metrics were estimated at varying values of identification threshold and numbers of technical replicates (Figure 3). Setting the identification threshold to -14 and considering one spot in the analysis, the sensitivity was 0.96 (95% credible interval [CrI]: 0.92 to 0.99), the positive predictive value was 0.94 (95%CrI: 0.92 to 0.96) and the accuracy was 0.90 (95%CrI: 0.87 to 0.94).

![__Figure 3. Evaluation of the performance of the reference mass spectra database for _Anopheles_ species identification using the test panel.__ (A) Sensibility and specificity determined at varying identification threshold considering one spot per specimen, (B) corresponding receiving operator characteristics curve, (C) sensibility and specificity determined at varying identification threshold considering four spots per specimen and (D) corresponding receiving operator characteristics curve. The shaded areas in panels A and C show the 95% credible interval around the median estimate of 1000 simulations. The dashed line in panels B and D shows the performance of a random classification.](analysis_files/figure-html/figure6_roc_curve-1.png)

<br/>
