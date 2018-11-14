#### Script for The Kara Davis Lab Developmental Classifier
#### by Christina Bligaard Pedersen
#### Updated on October 10, 2018

#### Loading libraries ####
libs <- c('flowCore', 'parallel', 'matlib', 'openxlsx', 'tools')
for (L in libs){
  if(!require(L, character.only = T)) {
    install.packages(L, dep = T, quiet = T)
    
    if (!require(L, character.only = T)) {
        suppressWarnings(source('http://bioconductor.org/biocLite.R')(L))
    }
  }
}

rm(L, libs)


#### User-defined parameters ####

# Specify location of gated pops
pops_dir <- '~/Documents/My_CYTOF_Files/Gated_ref_pops/'

# Define custom population names
pop_names <- c('HSC', 'Progenitor I', 'Progenitor II', 'Progenitor III', 'Pre-Pro-B', 
               'Pro-BI', 'Pro-BII', 'Pre-BI', 'Pre-BII', 'Immature BI', 'Immature BII',
               'Immature BIII', 'Mature BI', 'Mature BII', 'Early non-BI', 'Early non-BII', 'Mature Non-B')

# Channels that will be used for classification. Include the markers that were used for gating. Don't include markers that are uniformly present or absent. 
classifier_names <- c('CD19', 'CD45', 'CD20', 'CD24', 'CD34', 'CD38', 'CD127', 'CD179b', 'IgMi', 'Tdt')

# Set the threshold for unclassified cells
threshold <- c('simple', 10)    # Options e.g. c('simple', 10), c('farthest', 1.2), or 'neighbors'

# Set cores for parallel processing
nCores <- 1

# Specify directory with all the samples to run classification for
sample_folder <- '~/Documents/My_CYTOF_Files/'

# Specify the basename of the output directory (output is placed in a subdirecvtory of the sample directory)
output_folder <- 'Results'



#### Defining necessary functions ####

# Making summary tables (excel) and printing information about run to a text file
write_summary <- function(classifications, pop_names, output_folder, start_time, stop_time, pops_dir, 
                          sample_folder, threshold, classifier_names) {
  
  # Making summary tables
  table1 <- NULL; table2 <- NULL
  for (i in 1:length(classifications)) {
    classes <- table(factor(classifications[[i]]$MahID, levels = 0:length(pop_names)))
    table1 <- rbind(table1, classes)
    table2 <- rbind(table2, classes/sum(classes)*100)
  }
  rownames(table1) <- rownames(table2) <- names(classifications); colnames(table1) <- colnames(table2) <- c('Unclassified', pop_names)
  write.xlsx(x = list('Counts' = table1, 'Percent' = table2), file = file.path(output_folder, 'Classification_Summary.xlsx'), row.names = T)
  print('Finished writing excel output file.')
  
  
  
  # Writing info about run to a file
  write(paste0('The Developmental Classifier started parsing the input cells at ', start_time, ', and finished at ', stop_time, '.\n\n',
               'A total of ', sum(table1), ' cells from ', length(classifications), ' files were parsed and ', round(sum(table1[,1])/sum(table1)*100, 2), ' % of these were left unclassified.\n\n',
               'The reference file directory was ', pops_dir, ' and the number of reference populations was ', length(pop_names), '.\n\n',
               'The input directory was ', sample_folder, ' and the output directory was ', output_folder, '.\n\n',
               'The threshold(s) for unclassified cells was ', threshold, ' and the number of markers used was ', length(classifier_names), '.\n\n',
               'List of reference files:\n', paste(list.files(pops_dir, pattern = paste0('^(', paste(sprintf('%02d', 1:length(pop_names)), collapse = '|'), ').+.fcs*')), collapse = ', '), '\n\n',
               'List of input files:\n', paste(names(classifications), collapse = ', '), '\n\n',
               'List of applied markers:\n', paste(classifier_names, collapse = ', ')), file = file.path(output_folder, 'run_info.txt'))
}



# Making a single FCS file per input sample containing the raw data and columns containing information about classification and make one file per population
write_out_FCS <- function(file_to_classify, fcs, class_stats, output_folder, pop_names) {
  
  # Making main FCS with all events
  fcs_edited <- cbind(fcs@exprs, class_stats)
  
  output_file <- paste0(output_folder, '/Classified_', basename(file_to_classify))
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  write_modified_FCS(as.matrix(fcs_edited), output_file, channel_descriptions = fcs@parameters$desc,
                     reference_description = flowCore::description(fcs))
  
  cat('Main output FCS files written for', basename(file_to_classify), '\n')
  
  
  # Making population FCSs
  population_folder <- paste0(output_folder, '/', file_path_sans_ext(basename(file_to_classify)), '_pops')
  if (file.exists(population_folder)) {
    unlink(population_folder, recursive = T)
  }
  dir.create(population_folder, showWarnings = F)
  
  pop_names <- c('Unclassified', pop_names)
  
  for (pop in unique(class_stats$MahID)) {
    output_file <- paste0(population_folder, '/', sprintf('%02d', pop), '_', gsub(' ', '_', pop_names[pop+1]), '_', basename(file_to_classify))
    
    fcs_pop <- fcs_edited[fcs_edited[,'MahID']==pop,]
    write_modified_FCS(as.matrix(fcs_pop), output_file, channel_descriptions = fcs@parameters$desc,
                       reference_description = flowCore::description(fcs))
  }
  cat('Wrote', length(unique(class_stats$MahID)), 'population FCS files for', basename(file_to_classify), '\n')
}


# Write a modified FCS File (Code by Dmitry Tebaykin)
write_modified_FCS <- function(fcs_exprs, fcs_name, channel_descriptions = NULL, 
                               reference_description = NULL, old_description = NULL) {
  if (requireNamespace('flowCore', quietly = TRUE)) {
    if (is.matrix(fcs_exprs)) {
      
      # Build metadata for FCS file
      pd <- c()  # 'params' phenoData
      dl <- list()  # 'description' list
      
      dl[['$DATATYPE']] <- 'F'
      
      if (!is.null(reference_description)) {
        if (!is.null(reference_description[['$DATE']])){
          dl[['$DATE']] <- reference_description[['$DATE']]
        }
        if (!is.null(reference_description[['$BTIM']])){
          dl[['$BTIM']] <- reference_description[['$BTIM']]
        }
        if (!is.null(reference_description[['$ETIM']])){
          dl[['$ETIM']] <- reference_description[['$ETIM']]
        }
        if (!is.null(reference_description[['$CYT']])){
          dl[['$CYT']] <- reference_description[['$CYT']]
        }
        if (!is.null(reference_description[['$CYTSN']])){
          dl[['$CYTSN']] <- reference_description[['$CYTSN']]
        }      
      }
      
      for (c in 1:ncol(fcs_exprs)) {
        c_name <- colnames(fcs_exprs)[c]
        c_desc <- colnames(fcs_exprs)[c]
        if (!is.null(channel_descriptions)){
          if (!is.na(channel_descriptions[c])){
            c_desc <- channel_descriptions[c]  
          }
        }
        
        c_min <- floor(min(0, min(fcs_exprs[, c])))  # Prevents flowCore from shifting range
        c_max <- ceiling(max(fcs_exprs[, c]))
        c_rng <- c_max - c_min + 1
        
        pl <- matrix(c(c_name, c_desc, c_rng, c_min, c_max), nrow = 1)
        colnames(pl) <- c('name', 'desc', 'range', 'minRange', 'maxRange')
        rownames(pl) <- paste('$P', c, sep = '') 
        pd <- rbind(pd, pl)
        
        if (!is.null(reference_description[[paste0('P', c, 'DISPLAY')]])){
          dl[[paste0('P', c, 'DISPLAY')]] <- reference_description[[paste0('P', c, 'DISPLAY')]]
        } 
        
        if (!is.null(reference_description[[paste0('$P', c, 'G')]])){
          dl[[paste0('$P', c, 'G')]] <- reference_description[[paste0('$P', c, 'G')]]
        } 
        
        if (!is.null(reference_description[[paste0('$P', c, 'R')]])){
          dl[[paste0('$P', c, 'R')]] <- reference_description[[paste0('$P', c, 'R')]]
        } else {
          dl[[paste('$P', c, 'R',sep = '')]] <- toString(c_rng); # Range
        }
        
        if (!is.null(reference_description[[paste0('$P', c, 'B')]])){
          dl[[paste0('$P', c, 'B')]] <- reference_description[[paste0('$P', c, 'B')]]
        } else {
          dl[[paste('$P', c, 'B', sep = '')]] <- '32';      # Number of bits
        }
        
        if (!is.null(reference_description[[paste0('$P', c, 'E')]])){
          dl[[paste0('$P', c, 'E')]] <- reference_description[[paste0('$P', c, 'E')]]
        } else { 
          dl[[paste('$P', c, 'E', sep = '')]] <- '0,0';      # Exponent
        }
        
        dl[[paste('$P', c, 'N', sep = '')]] <- c_name;	    # Name
        dl[[paste('$P', c, 'S', sep = '')]] <- c_desc;	    # Desc	
      }	
      
      if (!is.null(old_description)) {
        if (!is.null(old_description[['$CYT']])){
          dl[['$CYT']] <- old_description[['$CYT']]
        }		  
        if (!is.null(old_description[['$DATE']])){
          dl[['$DATE']] <- old_description[['$DATE']]
        }
        if (!is.null(old_description[['$BTIM']])){
          dl[['$BTIM']] <- old_description[['$BTIM']]
        }
        if (!is.null(old_description[['$ETIM']])){
          dl[['$ETIM']] <- old_description[['$ETIM']]
        }
      }
      
      fcs_exprs <- flowCore::flowFrame(fcs_exprs, 
                                       as(data.frame(pd), 'AnnotatedDataFrame'), 
                                       description = dl) 
    }
    
    suppressWarnings(flowCore::write.FCS(fcs_exprs, fcs_name, what = 'numeric'))
  }
}


# Function to calculate, adjust, and invert the covariance matrix for the reference populations
get.cov <- function(data) {
  Sx <- cov(data)
  
  for (i in 1:nrow(Sx)) {
    for(j in 1:ncol(Sx)) {
      ifelse(i == j, Sx[i, j] <- max(Sx[i, j], 0.19), Sx[i, j] <- 0)
    }
  }
  
  Sx <- solve(Sx)
  return(Sx)
}


# Function to assign cells to populations
classify.cell <- function(cell, threshold) {
  
  # Calculate Mahalanobis distance to each population
  mah_dist <- c()
  for (population in 1:length(pop_names)) {
    mah_dist <- c(mah_dist, sqrt(mahalanobis(cell, get(paste0('Pop', population, '_means')), get(paste0('Pop', population, '_cov')), inverted = TRUE)))
  }

  # Assign label based on Mahalanobis distance
  cell_dist <- min(mah_dist)
  mah_id <- which(mah_dist == cell_dist)

  # Checking if cell should be unclassified
  if (is.numeric(threshold) & length(threshold)==1) {
    if (cell_dist > threshold) mah_id <- 0
  } else if (is.numeric(threshold)) {
    if (cell_dist > threshold[mah_id]) mah_id <- 0
  } else {
    if (abs(diff(order(mah_dist)[1:2])) != 1) mah_id <- 0
  }

  mah.stat <- c(mah_dist, mah_id)

  # Return calculated distances and the final label
  return(mah.stat)
}



#### Actual data processing ####

# A vector for saving the distance from the most distant cell in each reference population
farthest <- c()

# Reading defined populations
# Calculate mean and covariance matrix for each population
for (population in 1:length(pop_names)) {
  
  # Extract expression data for the current population and transform
  fcs <- read.FCS(list.files(pops_dir, pattern = paste0('^', sprintf('%02d', population), '.+.fcs*'), full.names = T), transformation = F)
  population_data <- fcs@exprs; colnames(population_data) <- fcs@parameters$desc
  
  # Checking if channel names are in all manually gated files
  if (!all(classifier_names %in% colnames(population_data))) {
    stop(paste0('Please check the channel names. Not all channel names exist in the provided FCS file for population ', population))
  }
  
  population_data <- asinh(population_data[,classifier_names]/5)
  
  # Calculate population mean
  assign(paste0('Pop', population, '_means'), colMeans(population_data[,classifier_names]))
  
  # Calculate the population Sx
  assign(paste0('Pop', population, '_cov'), get.cov(population_data[,classifier_names]))
  
  if (threshold[1] == 'farthest') {
    farthest <- c(farthest, (as.numeric(threshold[2]) * max(unname(apply(population_data, 1, function(y) { 
      return(sqrt(mahalanobis(y, get(paste0('Pop', population, '_means')), get(paste0('Pop', population, '_cov')), inverted = TRUE))) })))))
  }
  
  cat(paste0('Done reading and processing the file "', list.files(pops_dir, pattern = '.fcs')[population], '". Data will be used to represent the population you have called "', pop_names[population], ' cells".\n'))
}

# Checking validity of threshold before moving on
if (threshold[1] == 'farthest') {
  threshold <- farthest
} else if (threshold[1] == 'simple') {
  threshold <- as.numeric(threshold[2])
} else if (threshold != 'neighbors') {
  stop('Your specified threshold is not allowed! Please correct this before re-running.')
}

# Create output directory in case it doesn't exist - and check if sample folder is valid
output_folder <- file.path(sample_folder, output_folder)
if (dir.exists(sample_folder)) {
  dir.create(output_folder, showWarnings = FALSE)
} else {
  stop(paste0('Please specify a proper sample folder. The folder ', sample_folder, ' does not exist.'))
}

# Performing the classification on a set of files
classifications <- list(); start_time <- Sys.time()
for (file_to_classify in list.files(sample_folder, full.names = T, pattern = '.fcs')) {
  
  # Read file and transform
  fcs <- read.FCS(file_to_classify, transformation = F)
  to_classify <- fcs@exprs; colnames(to_classify) <- gsub('.+?_(.+)', '\\1', fcs@parameters$desc)
  to_classify <- asinh(to_classify[,classifier_names]/5)

  # Make cluster and feed information
  cl <- parallel::makePSOCKcluster(nCores)
  clusterExport(cl, c(paste0('Pop', seq(1,length(pop_names)), '_means'), paste0('Pop', seq(1,length(pop_names)), '_cov'), 'threshold', 'pop_names'))
  
  # Classifying all the cells
  class_stats <- as.data.frame(t(parApply(cl, to_classify, 1, classify.cell, threshold)))
  colnames(class_stats) <- c(paste0('Mah-Pop', 1:length(pop_names)), 'MahID')

  # Stop parallel processing
  stopCluster(cl)
  
  # Save results
  classifications[[basename(file_to_classify)]] <- class_stats
  cat('Classification complete for', basename(file_to_classify), '\n')
  
  # Making output FCS files
  write_out_FCS(file_to_classify, fcs, class_stats, output_folder, pop_names)
  
}

# Making summary tables (excel) and printing info file
write_summary(classifications, pop_names, output_folder, start_time, Sys.time(), pops_dir, 
              sample_folder, threshold, classifier_names)
