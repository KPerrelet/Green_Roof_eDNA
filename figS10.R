################################################################################
#################################### fig s8 ####################################
################################################################################
library(ggplot2)
library(ggpubr)

# Import data
shared.df <- read.csv("skyline_sharedmat.csv", check.names = F)
taxo.df <- read.csv("skyline_taxo.csv")

# Select columns based on taxa in the taxonomy dataframe and remove columns with zero sum (species absent in all sites)
shared.df <- shared.df %>%
  dplyr::select(any_of(taxo.df$taxa)) %>%
  dplyr::select_if(colSums(.) > 0)
  
# Convert values greater than 1 in the shared species matrix to 1 (indicating presence or absence)
shared.df[shared.df > 1] <- 1

# Create a list of species and their corresponding site counts (i.e., how many sites they appear in)
shared.list <- data.frame(taxa = colnames(shared.df),
                          no_sites = colSums(shared.df))

# Add mobility type information for each species (Epigeic, Edaphic, or Flying)
shared.list$mobility <- taxo.df[taxo.df$taxa %in% shared.list$taxa, ]$mobility

ggplot(shared.list, aes(x = no_sites, fill = mobility)) +
  geom_histogram(aes(y = stat(density)), alpha = 0.3, bins = 54, color = "black", position = "identity") +
  geom_density(alpha = 0.5, linetype = "solid") +
  scale_fill_manual(name = "", labels = c("Epigeic", "Edaphic", "Flying"), values = c("#FE4A49", "#D4B483", "#2AB7CA")) +
  labs(x = "Number of sites", y = "Shared species density") +
  theme_pubr(base_size = 16) +
  facet_wrap(~ mobility, scales = "free_y") + 
  theme(strip.background = element_blank(), 
        strip.text = element_blank())

