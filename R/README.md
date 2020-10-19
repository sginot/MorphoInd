# MorphoInd
R package to compute various morphological indices (Log-Size Index, Variability Size Index, VSI*) and confidence intervals.

WARNING: THESE FUNCTIONS ARE STILL BEING MODIFIED, ALTHOUGH THEY SHOULD BE WORKING SOME OPTIONS MAY NOT YET BE OPERATIONAL. RESULTS ARE NOT GUARANTEED.

The 'morpho_indices' function computes the Log-Size Index (LSI) (Simpson 1941), Variability Size Index (VSI) (Uerpmann 1982) and VSI* (modified version of VSI, Escarguel 2008). Additionally, it can produce bootstraped values for confidence intervals, either with parametric or non-parametric approaches.

Input data should be formatted as a list of N samples to be compared, including one reference sample. Elements of the list can either be matrices with varying numbers of individuals as rows and j morphological variables as columns; or be matrices with three rows (first row = "N", number of individuals ; second row = "Mean", mean value for variable ; third row = "St-Dev", standard deviation for variable) and j morphological variables as columns. Both type of matrices (i.e. with individual measurements and with summary statistics) can be combined in the same list.

The 'plot_mo_ind' function allows flexible plotting of the output of 'morpho_indices'. The input of 'plot_mo_ind' should be the entire output of 'morpho_indices'. Filtering of index, samples and morphological variables to be plotted is done using various arguments of the function.
