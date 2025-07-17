################################################################################
##################################### fig 4 ####################################
################################################################################
library(randomForest)
library(ggplot2)
library(ggpubr)
library(vegan)
library(dplyr)

# Set seed for reproducibility
set.seed(1)

# Import data
spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
taxo.df <- read.csv("skyline_taxo.csv")
meta <- c("site_no", "type", "sp_richness")

# Calculate species richness for different mobility modes
spmat.df$richness_flying <- specnumber(spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "flying", ]$taxa)])
spmat.df$richness_epigeic <- specnumber(spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "epigeic", ]$taxa)])
spmat.df$richness_edaphic <- specnumber(spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "edaphic", ]$taxa)])
# Calculate total species richness excluding aquatic species
spmat.df$richness_total <- spmat.df$richness_flying + spmat.df$richness_epigeic + spmat.df$richness_edaphic 

# Subset only green roof data and load green roof covariates
roof.df <- spmat.df[spmat.df$type == "Roof", ]
roof.covariates.df <- read.csv("skyline_covmat_greenroof.csv")
roof.covariates.df$pv <- ifelse(roof.covariates.df$pv == 0, FALSE, TRUE)
roof.covariates.df$substrate_type <- as.factor(roof.covariates.df$substrate_type)
roof.df <- merge(roof.df, roof.covariates.df, by = "site_no")

# Add ground-level species richness for comparison
roof.df$ground_richness <- spmat.df[spmat.df$type == "Ground", ]$richness_total

# Define predictor variables and standardization
covariates <- c("vegetation", "height", "depth", "age_roof", "area",
                "green_frac_50", "green_frac_500", "distance_gr", "ground_richness")
roof.df[, covariates] <- scale(roof.df[, covariates], center = T, scale = T) # scale to zero
covariates <- append(covariates, c("pv", "substrate_type"))

# Initialize empty data frames for GLM and RF results
res.lm <- data.frame()
res.rf <- data.frame()

# Loop through different species richness response variables
for (response in c("richness_flying",
                   "richness_epigeic", 
                   "richness_edaphic", 
                   "richness_total")) {
  for (covariate in covariates) {   # Perform univariate GLMs for each covariate
    response.lm <- glm(roof.df[, response] ~ roof.df[, covariate], family = "poisson", data = roof.df)
    res.lm <- rbind(res.lm, data.frame(response = response, 
                                       parameter = covariate, 
                                       coef = response.lm$coefficients[-1],
                                       CI2.5 = confint(response.lm)[2, 1], 
                                       CI97.5 = confint(response.lm)[2, 2], 
                                       pvalue = summary(response.lm)$coefficients[2, 4]))
    
  }
  
  # Fit a random forest model with all covariates
  formula <- as.formula(paste(response, paste(covariates, collapse = "+"), sep = " ~ "))
  rf <- randomForest(formula,
                     data = roof.df,
                     importance = TRUE,
                     ntree = 50000,
                     na.action = na.omit)
  rf.df <- as.data.frame(importance(rf, type = 1))
  rf.df$covariate <- rownames(rf.df)
  rf.df$response <- response
  res.rf <- rbind(res.rf, rf.df)
}

# Reorder factor levels for plotting
res.lm$response <- factor(res.lm$response, 
                          levels = c("richness_total", "richness_flying", "richness_epigeic", "richness_edaphic")) 
res.rf$response <- factor(res.rf$response, 
                          levels = c("richness_total", "richness_flying", "richness_epigeic", "richness_edaphic")) 

# Rename variables for better readability in plots
res.lm <- res.lm %>%
  mutate(parameter = recode(parameter, 
                            age_roof = "Roof age",
                            depth = "Depth",
                            substrate_type = "Stony substrate",
                            vegetation = "Vegetation",
                            ground_richness = "Ground-level richness",
                            green_frac_500 = "% green (500 m)",
                            pv = "Solar panels",
                            green_frac_50 = "% green (50 m)",
                            height = "Height",
                            distance_gr = "Distance to roof",
                            area = "Area"))

# Order covariates based on mean coefficients
ordered_covariates <- res.lm %>%
  group_by(parameter) %>%
  dplyr::summarize(mean_coef = mean(coef)) %>%
  dplyr::arrange(-mean_coef) %>%
  dplyr::pull(parameter)

coef.plot <- ggplot(res.lm, aes(y = coef, x = parameter, color = response)) +
  geom_tile(data = res.lm[res.lm$parameter %in% c("Depth",
                                                  "Vegetation",
                                                  "Distance to roof",
                                                  "% green (500 m)",
                                                  "Height",
                                                  "Stony substrate"), ], 
            aes(y = 0, x = parameter, height = Inf, width = 1), fill =  "grey", alpha = 0.1, color = NA) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = CI2.5, ymax = CI97.5), width = 0, linewidth = 1.07, position = position_dodge(width=0.7)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  labs(x = "", y = "Standardized coefficient") +
  scale_color_manual(name = "", labels = c("All", "Flying", "Epigeic", "Edaphic"), values = c("#000000", "#2AB7CA", "#FE4A49", "#D4B483")) +
  theme_pubr(base_size = 18) +
  scale_x_discrete(limits = ordered_covariates) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14), 
        legend.text = element_text(size = 18))

# Rename RF covariates for better readability
res.rf <- res.rf %>%
  mutate(covariate = recode(covariate, 
                            age_roof = "Roof age",
                            depth = "Depth",
                            substrate_type = "Stony substrate",
                            vegetation = "Vegetation",
                            ground_richness = "Ground-level richness",
                            green_frac_500 = "% green (500 m)",
                            pv = "Solar panels",
                            green_frac_50 = "% green (50 m)",
                            height = "Height",
                            distance_gr = "Distance to roof",
                            area = "Area"))

ordered_covariates_per_response <- res.rf %>%
  group_by(response, covariate) %>%
  dplyr::summarize(mean_inc_mse = mean(`%IncMSE`)) %>%
  dplyr::arrange(response, mean_inc_mse)

# Create variable importance plots from RF results
plots <- list()
for (response_group in levels(res.rf$response)) {
  ordered_covariates_for_response <- ordered_covariates_per_response %>%
    dplyr::filter(response == response_group) %>%
    dplyr::pull(covariate)
  
  response_data <- res.rf %>%
    dplyr::filter(response == response_group) %>%
    dplyr::mutate(covariate = factor(covariate, levels = ordered_covariates_for_response))
  
  
  color <- ifelse(response_group ==  "richness_total", "#000000", 
                  ifelse(response_group ==  "richness_flying", "#2AB7CA", 
                         ifelse(response_group ==  "richness_epigeic", "#FE4A49", "#D4B483")))
  name <- ifelse(response_group ==  "richness_total", "All", 
                 ifelse(response_group ==  "richness_flying", "Flying species", 
                        ifelse(response_group ==  "richness_epigeic", "Epigeic species", "Edaphic species")))
  
  
  plot <- ggplot(response_data, aes(y = covariate, x = `%IncMSE`)) +
    geom_vline(xintercept = 0) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6, fill = color) +
    geom_hline(yintercept = 0) +
    labs(y = "", x = ifelse(name %in% c("Epigeic species", "Edaphic species"), "Increase in accuracy", "")) +
    theme_pubr(base_size = 18) +
    theme(axis.text.x = element_text(size = 14), 
          plot.margin = unit(c(0.2,0.3,0.2,0.2), "cm")) + 
    ggtitle(name) 
  
  plots[[response_group]] <- plot
}

ggarrange(coef.plot, 
          ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                    ncol = 2, nrow = 2, common.legend = TRUE), 
          nrow = 2, common.legend = T, 
          heights = c(0.8, 1), 
          labels = c("A", "B"), 
          font.label = list(size = 24), vjust = 0)
