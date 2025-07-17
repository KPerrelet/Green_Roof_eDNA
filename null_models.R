################################################################################
################## shared species composition and null models ##################
################################################################################
library(dplyr)
library(betapart)
library(plyr)
library(vegan)

set.seed(1) # Set random seed for reproducibility

# Function to create a null model by pairing a site with a randomly selected site from another type (random pairing)
create_null_model_pair <- function(full.df, site, n_iterations = 100) {
  
  # Initialize vectors to store results for the null model
  null_shared_sp <- numeric(n_iterations)
  null_jaccard <- numeric(n_iterations)
  
  # Loop over the number of iterations to generate random null models
  for (i in 1:n_iterations) {
    # Create a null dataframe by selecting the site data and pairing with random sites from the opposite type
    null.df <- rbind.fill(full.df %>%
                            filter(type == "Roof", site_no == site) %>%
                            dplyr::select(-any_of(meta)) %>%
                            dplyr::mutate(across(everything(), ~ ifelse(. > 1, 1, .))) %>%
                            dplyr::select(where(~ sum(.) > 0)), 
                          full.df %>%
                            filter(type == "Ground", site_no == sample(setdiff(full.df$site_no, site), 1)) %>%
                            dplyr::select(-any_of(meta)) %>%
                            dplyr::mutate(across(everything(), ~ ifelse(. > 1, 1, .))) %>%
                            dplyr::select(where(~ sum(.) > 0))) 
    
    # Replace NA values with zeros
    null.df[is.na(null.df)] <- 0
    
    # Compute Jaccard dissimilarity and the number of shared speciesfor the null model
    null_jaccard[i] <- as.numeric(beta.pair(null.df, "jaccard")[3])
    null.shared.df <- as.data.frame(null.df[, colSums(null.df) > 1])
    null_shared_sp[i] <- ncol(null.shared.df)
  }
  return(list(null_shared_sp, null_jaccard))
}

# Function to create a full null model by randomly selecting species based on habitat type (semi-randomization)
create_null_model_full <- function(green_roof, ground_level, full.df, n_iterations = 100) {
  
  # Prepare the dataframes for Roof and Ground type data
  roof.df <- full.df %>%
    dplyr::filter(type == "Roof") %>%
    dplyr::select(-any_of(meta)) %>%
    dplyr::select_if(colSums(.) > 0)
  
  ground.df <- full.df %>%
    dplyr::filter(type == "Ground") %>%
    dplyr::select(-any_of(meta)) %>%
    dplyr::select_if(colSums(.) > 0)
  
  # Remove metadata columns from the full dataframe
  full.df <- full.df %>% 
    dplyr::select(-any_of(meta))
  
  # Initialize vectors to store the results for the null model
  null_shared_sp <- numeric(n_iterations)
  null_jaccard <- numeric(n_iterations)
  
  # Loop over the number of iterations to generate random null models
  for (i in 1:n_iterations) {
    
    # Initialize a null dataframe
    null.df <- matrix(ncol = ncol(full.df), nrow = 2) 
    colnames(null.df) <- colnames(full.df)
    null.df[is.na(null.df)] <- 0
    
    # Randomly select species for the Roof and Ground-level types and assign them to the null dataframe
    null.df[1, which(colnames(null.df) %in% sample(colnames(roof.df), specnumber(green_roof)))] <- 1
    null.df[2, which(colnames(null.df) %in% sample(colnames(ground.df), specnumber(ground_level)))] <- 1
    
    # Compute Jaccard dissimilarity and the number of shared speciesfor the null model
    null_jaccard[i] <- as.numeric(beta.pair(null.df, "jaccard")[3])
    null.shared.df <- as.data.frame(null.df[, colSums(null.df) > 1])
    null_shared_sp[i] <- ncol(null.shared.df)
  }
  return(list(null_shared_sp, null_jaccard))
}

# Import data
taxo.df <- read.csv("skyline_taxo.csv")
spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
meta <- c("site_no", "type", "sp_richness")
shared.df <- read.csv("skyline_shareddist.csv")

# Initialize an empty vector to store shared species information
shared_sp <- c()

# Loop over each unique site and compute shared species and Jaccard index
for (site in unique(spmat.df$site_no )) {
  
  # Prepare a raw dataframe for the current site (filter and convert species data to binary)
  raw.df <- spmat.df %>%
    filter(site_no == site) %>%
    dplyr::select(-any_of(meta)) %>%
    dplyr::mutate(across(everything(), ~ ifelse(. > 1, 1, .))) %>%
    dplyr::select(where(~ sum(.) > 0))
  
  # Initialize a counter for shared species
  i <- 0
  for (sp in raw.df) {
    i <- i + 1
    if (sp[1] > 0 &  sp[2] > 0) {  # Check if species is present in both Roof and Ground
      shared_sp <- append(shared_sp, colnames(raw.df)[i]) # Append shared species to the list
    }
  }
  
  # Compute Jaccard dissimilarity and the number of shared species 
  shared.df[shared.df$site_no == site, "jaccard"] <- as.numeric(beta.pair(raw.df, "jaccard")[3])
  raw.shared <- as.data.frame(raw.df[, colSums(raw.df) > 1])
  shared.df[shared.df$site_no == site, "sp_shared"] <- ncol(raw.shared)
  
  # Compute the Jaccard dissimilarity for the null models by pairing the site with a random site (random pairing)
  null_jaccard_distances <- create_null_model_pair(spmat.df, site = site) # We keep the site with the type specified but we randomly pair it with a site from a different type
  shared.df[shared.df$site_no == site, "sp_shared.pair"] <- mean(null_jaccard_distances[[1]])
  shared.df[shared.df$site_no == site, "jaccard.pair"] <- mean(null_jaccard_distances[[2]])
  
  # Compute the Jaccard dissimilarity for the null models by using random species selections (semi-randomization)
  null_jaccard_distances <- create_null_model_full(raw.df[1, ], raw.df[2, ], full.df = spmat.df)
  shared.df[shared.df$site_no == site, "sp_shared.rand"] <- mean(null_jaccard_distances[[1]])
  shared.df[shared.df$site_no == site, "jaccard.rand"] <- mean(null_jaccard_distances[[2]])
}  

# Create a dataframe containing the species shared between Roof and Ground-level types
shared <- spmat.df[, which(colnames(spmat.df) %in% c("site_no", unique(shared_sp)))]
for (site in unique(shared$site_no)) {
  site_ij <- shared[shared$site_no == site, ] %>%
    dplyr::select(-site_no)
  for (i in 1:ncol(site_ij)) {
    if (site_ij[1, i] == 0 | site_ij[2, i] == 0) { # Keep only species present in both Roof and Ground-level
      shared[shared$site_no == site, i+1] <- 0
    }
  }
}

# Aggregate data per site (sum of species across columns)
shared <- shared %>%
  dplyr::select(-site_no) %>%
  aggregate(by = list(shared$site_no), FUN = sum)
colnames(shared)[1] <- "site_no"

# Merge the aggregated shared data with the original shared dataframe
shared.df <- left_join(shared.df, shared, by = "site_no")

# Save the final shared species matrix to a CSV file
write.csv(shared.df, "skyline_sharedmat.csv")
