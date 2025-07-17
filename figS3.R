################################################################################
#################################### fig s3 ####################################
################################################################################
library(ggpubr)

# Import data
spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
roof.covariates.df <- read.csv("skyline_covmat_greenroof.csv")
roof.covariates.df$pv <- ifelse(roof.covariates.df$pv == 0, FALSE, TRUE)
roof.covariates.df$substrate_type <- as.factor(roof.covariates.df$substrate_type)

# Subset to modelled green roofs
roof.covariates.df <- roof.covariates.df[roof.covariates.df$site_no %in% spmat.df$site_no, ]
roof.covariates.df$ground_level_richness <- spmat.df[spmat.df$site_no %in% roof.covariates.df$site_no & spmat.df$type == "Ground", ]$sp_richness

# Define the list of covariates (columns) to include in the correlation matrix
covariates <- c("vegetation", "height", "depth", "age_roof", "area", "substrate_type",
                "green_frac_50", "green_frac_500", "distance_gr", "pv", "ground_level_richness")

covariates.df <- roof.covariates.df[, covariates] 

# Rename the covariates for visualization purposes
covariates.df <- covariates.df %>% 
  dplyr::rename(
    "Roof age" = age_roof,
    "Depth" = depth,
    "Stony substrate" = substrate_type,
    "Vegetation" = vegetation,
    "Ground-level richness" = ground_level_richness,
    "% green (500 m)" = green_frac_500,
    "Solar panels" = pv,
    "% green (50 m)" = green_frac_50,
    "Height" = height,
    "Distance to roof" = distance_gr,
    "Area" = area
  )

plots <- list()
for (var in colnames(covariates.df)) {
  if (class(covariates.df[, var]) %in% c("factor", "logical")) {
    p <- ggplot(covariates.df, aes(x = as.factor(.data[[var]]))) +
      geom_bar() +
      labs(x = "Value", y = "Count", title = var) +
      theme_pubr(base_size = 12) +
      theme(title = element_text(size = 12))
  } else {
    p <- ggplot(covariates.df, aes(x = as.numeric(.data[[var]]))) +
      geom_histogram(bins = 10) +
      labs(x = "Value", y = "Count", title = var) +
      theme_pubr(base_size = 12) +
      theme(title = element_text(size = 12))
  }
  
  
  plots[[var]] <- p
}

ggarrange(plotlist = plots, ncol = 6, nrow = 2)