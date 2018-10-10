# The Kara Davis Lab Developmental Classifier

The Developmental Classifier was made to assign cells of unknown origin to their closest defined population. One example could be to define populations using a healthy sample and then assigning individual cells of B-ALL patients to these populations. The current version of the code was implemented by [Christina Bligaard Pedersen](https://github.com/cbligaard) with great inspiration from [Dmitry Tebaykin](https://github.com/dtebaykin) - the original idea and code were developed by [Nikolay Samusik](https://github.com/nsamusik) and [Zinaida Good](https://github.com/zinagood) and presented [here](http://doi.org/10.1038/nm.4505 "Nature Medicine Article").

The current implementation is one, rather long, script, which contains all the functions - but don't worry, the guide below should take you through it and the script is made in one to make everything easier for you as well.

### The way it works
1. Perform manual gating of reference populations and save each population as a separate FCS. *This is done before using the classifier itself.*

2. Provide a list of the markers you used for manual gating to the classifier script. You will use all of these markers for the classification of new samples.

3. Read manually gated populations one at a time, ArcSinh-transform with a co-factor of 5, and calculate a) the mean expression of each of the included markers, and b) the covariance matrix. The covariance matrix is then corrected such that all off-diagonal values = 0, and all diagonal values are set to the max of their value and 0.19, and finally, it is inverted for computational efficiency.

4. For each new sample (one per FCS file) to classify: Read FCS, reduce to include the selected markers. Apply ArcSinh-transform with a co-factor of 5. Then for each cell:
	1. Calculate Mahalanobis distance to each of the reference populations.
	2. Assign cell to the population which has the shortest distance if that distance is less than the defined threshold, otherwise assign to ‘unclassified’ bin.

5. The output is placed in a subdirectory of the sample folder. There are a set of files:
	1. Run information file (run_info.txt).
	2. Excel file with two sheets (Classification_Summary.xlsx):
		1. The counts of each cell type in each input file.
		2. The percentage of each cell type in each input file.
	3. An FCS file per input file with the added information of the distances from each cell to each of the reference populations and the assigned population.
	4. A folder per input FCS file containing an FCS file per reference population (contains the cells assigned that population). No files are made for empty populations.


### How to run
1. You need to install the packages [flowCore](http://bioconductor.org/packages/release/bioc/html/flowCore.html), [parallel](https://www.rdocumentation.org/packages/parallel/versions/3.5.1), [openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html), [tools](https://www.rdocumentation.org/packages/tools/versions/3.5.1), and [matlib](https://cran.r-project.org/web/packages/matlib/index.html). The code should do this for you, but some users may have trouble. 

2. Go to the section 'User-defined parameters'. Here, set the location of your reference files (they should all be in the same directory) and define some 'prettier' population names if you wish to (but please don't use special characters such as '+' or '/'). **But do note, that your reference file names should be numbered 1-n<sub>pops</sub> (examples are 01\_HSC.fcs and 12\_immatureB\_38-.fcs)**. Do not use any strange characters in you file names such as '+'. 

3. Define the markers you want to use for classification (those used for gating of the reference). If in doubt about how to set the classifier.names, try to take a look at the header of one of your reference FCS files by running this code: ```read.FCS(list.files(pops.dir, full.names = T)[1], transformation = F)@parameters$desc```

4. Set the treshold for unclassified cells (see [below](#setting-a-threshold-for-classification)).

5. Set the number of cores to use for parallel processing. You may use 1 (always possible), but to speed up you may want to use 2-4 cores. If you run `detectCores()`, you can check the number of cores available - but you don't want to use all of them. If your computer gets really noisy during the process you may also want to reduce the number of cores. 

6. Specify directory with all the samples to run classification for. All samples in the folder will be included.

7. Specify the directory to place the new output directory in, and specify the name of the output directory to be generated.

8. Run the rest of the code (click 'Source' in the top right corner of the code field in R-studio) and wait for the results. 


##### Setting a threshold for classification
Another big question is how to set the threshold for classification. One imagine cases where we want to classify all cells, but in some cases we may also want to keep a bin for unclassified cells. But when is a cell far enough from its nearest population to be put in such a bin? Again, it may depend on the case - leukemic cells are generally 'weird-looking' and should perhaps be left unclassified relatively often. Some options could be:

1. Just set it and go! If you don’t like the output fraction of unclassified, rerun.

2. Classify all cells, then look at the minimum distances and unclassify a certain fraction of them again.

3. Calculate the threshold based on healthy reference populations. E.g. A) what is the Xth percentile of distances from the reference cells to their respective centroids, or B) an extrapolation - when a cell is Y times further from the assigned centroid than the farthest reference cell (c<sub>farthest</sub>), unclassify it.

4. Compare the distance from a cell to its closest reference population and to the second-closest. If the difference between these two is too small, the classification is uncertain.

5. Look at *N* (= 2) closest populations - if these are not neighbors, leave the cell unclassified. Possibly only good for developmental settings in which the concept of neighbors is defined.

Currently, you can choose between options 1, 3B, and 5. 

* For option 1, define the `threshold` variable as a vector: `c('simple', X)`, where X is your selected number. Suggestion: The number of markers used for classification.
* For option 3B, define the `threshold` variable as a vector: `c('farthest', X)`, where X is the factor to multiple the distance from c<sub>farthest</sub> to the centroid with, to get the threshold.
* For option 5, define the `threshold` variable as a string: `neighbors`. *N* will be set to 2. However, this will likely lead to many cells being unclassified...

### Benchmarking results
Using four pooled healthy donors subjected to a variety of stimuli in a 10-fold cross-validation setup, in which the reference sample was iteratively defined using nine-tenths of the data and tested on the remaining one-tenth lead to an overall accuracy of 91.9 % (as presented [here](http://doi.org/10.1038/nm.4505 "Nature Medicine Article")). In all testing, the threshold was set to 11 (the number of markers used in the classification).

We tried to run the same benchmarking using a code with several ajustments - including an uncorrected covariance matrix (acc. = 89.2 %), or covariance matrices in which only off-diagonal (acc. = 88.0 %) or diagonal elements (acc. = 92.4 %) were changed. Curiously, when only corecring the diagonal elements, we could actually improve performance slightly - and interestingly, when changing the adjustment value from 0.19 to 0.4-0.6, the accuracy was increased to 93.4 %. Using these results it looks like only correcting the diagonal is a good idea. *However*, the correction of the off-diagonal elements matters a lot when studying a case with fewer reference data points (i.e. only a single sample). Here, not correcting the off-diagonal elements, may lead to populations with a relatively small set of reference data points being less likely to receive events, when classifying new samples. The correction of the covariance matrix generally belongs to the concept of [shrinkage](https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices#Shrinkage_estimation "Wikipedia link"). Another [method](https://www.rdocumentation.org/packages/corpcor/versions/1.6.9/topics/cov.shrink "cov.shrink") for shrinkage was also tested, but this did not improve performance. 

The cut-off for correction of diagonal elements should also be discussed. 0.19 may not be the optimal solution in all cases, and according to the cross-validation, it is not even the best setting for our data. **The effect of this value could also be tested on leukemic data.**

<br></br>
*Last updated by [Christina Bligaard Pedersen](https://github.com/cbligaard) on October 9, 2018.*
