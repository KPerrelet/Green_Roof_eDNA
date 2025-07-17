################################################################################
#################################### fig s6 ####################################
################################################################################
library(ggplot2)
library(ggpubr)
library(vegan)

spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
taxo.df <- read.csv("skyline_taxo.csv")
shared.df <- read.csv("skyline_sharedmat.csv")
meta <- c("site_no", "type", "sp_richness")
roof.covariates.df <- read.csv("skyline_covmat_greenroof.csv")
ground.covariates.df <- read.csv("skyline_covmat_groundlevel.csv")

# Calculate species richness for different mobility types
spmat.df$richness_flying <- specnumber(spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "flying", ]$taxa)])
spmat.df$richness_epigeic <- specnumber(spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "epigeic", ]$taxa)])
spmat.df$richness_edaphic <- specnumber(spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "edaphic", ]$taxa)])

# Subset data for roof and ground types
roof.df <- spmat.df[spmat.df$type == "Roof", ]
ground.df <- spmat.df[spmat.df$type == "Ground", ]

# Import covariate data for roof and ground sites
roof.df <- merge(roof.df, roof.covariates.df, by = "site_no")
ground.df <- merge(ground.df, ground.covariates.df, by = "site_no")

# Calculate species richness for shared species between roof and ground
shared.df$richness_flying <- specnumber(shared.df[, which(colnames(shared.df) %in% taxo.df[taxo.df$mobility == "flying", ]$taxa)])
shared.df$richness_epigeic <- specnumber(shared.df[, which(colnames(shared.df) %in% taxo.df[taxo.df$mobility == "epigeic", ]$taxa)])
shared.df$richness_edaphic <- specnumber(shared.df[, which(colnames(shared.df) %in% taxo.df[taxo.df$mobility == "edaphic", ]$taxa)])
shared.df <- merge(shared.df, ground.covariates.df, by = "site_no")

# Loop through different combinations of site type and mobility type to perform linear models and extract p-values
summary.df <- data.frame()
for (typ in c("Roof", "Ground", "Shared")) { 
  for (mob in c("richness_flying", "richness_epigeic", "richness_edaphic")) {
    if (typ == "Roof") {
      # Perform linear regression for roof data: relationship between richness fraction and grey fraction
      pvalue <- summary(lm(roof.df[, mob] / roof.df$sp_richness ~ grey_frac_500, roof.df))$coefficients[2, 4]
      summary.df <- rbind(summary.df, data.frame(frac = roof.df[, mob] / roof.df$sp_richness, 
                                                 type = "Green roof",
                                                 mobility = mob, 
                                                 grey_frac_500 = roof.df$grey_frac_500,
                                                 pvalue = pvalue,
                                                 significance = ifelse(pvalue <= 0.05, "Significant", "Unsignificant")))
      
    } 
    if (typ == "Ground") {
      pvalue <- summary(lm(ground.df[, mob] / ground.df$sp_richness ~ grey_frac_500, ground.df))$coefficients[2, 4]
      summary.df <- rbind(summary.df, data.frame(frac = ground.df[, mob] / ground.df$sp_richness,
                                                 type = "Ground-level",
                                                 mobility = mob,
                                                 grey_frac_500 = ground.df$grey_frac_500,
                                                 pvalue = pvalue,
                                                 significance = ifelse(pvalue <= 0.05, "Significant", "Unsignificant")))
    }
    if (typ == "Shared") {
      pvalue <- summary(lm(shared.df[, mob] / shared.df$sp_shared ~ grey_frac_500, shared.df))$coefficients[2, 4]
      summary.df <- rbind(summary.df, data.frame(frac = shared.df[, mob] / shared.df$sp_shared, 
                                                 type = "Shared",
                                                 mobility = mob, 
                                                 grey_frac_500 = shared.df$grey_frac_500,
                                                 pvalue = pvalue,
                                                 significance = ifelse(pvalue <= 0.05, "Significant", "Unsignificant")))
    }
  }
}   

ggplot(summary.df, aes(x = grey_frac_500, y = frac, color = mobility, linetype = significance)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", lwd = 1.3) +
  facet_wrap(~ type) +
  labs(x = "% grey spaces (500 m buffer)", y = "Fraction of total species richness") + 
  scale_color_manual(name = "", labels = c("Edaphic", "Epigeic", "Flying"), values = c("#D4B483", "#FE4A49", "#2AB7CA")) +
  theme_pubr(base_size = 18) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 16)) + 
  stat_cor(size = 5)
