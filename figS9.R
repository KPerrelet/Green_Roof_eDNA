################################################################################
#################################### fig s7 ####################################
################################################################################
library(MASS)
library(ggplot2)
library(ggpubr)

# Import data
spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
shared.df <- read.csv("skyline_sharedmat.csv", check.names = F)

# Order data by site number (important for matching data correctly)
spmat.df <- spmat.df[order(spmat.df$site_no), ]
shared.df <- shared.df[order(shared.df$site_no), ]

# Subset data for roof and ground types
roof.df <- spmat.df[spmat.df$type == "Roof", ]
ground.df <- spmat.df[spmat.df$type == "Ground", ]

# Create a dataframe containing richness for roof, ground, and shared species
all.alpha <- data.frame(roof.df$sp_richness,
                        ground.df$sp_richness,
                        shared.df$sp_shared)
colnames(all.alpha) <- c("richness_roof",
                         "richness_ground", 
                         "richness_shared")

# Perform a negative binomial regression to analyze the relationship between richness on roof and ground sites
summary(glm.nb(richness_roof ~ richness_ground, data = all.alpha))

ggplot(all.alpha, aes(x = richness_ground, y = richness_roof)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 1, alpha = 0.2) +
  geom_point(aes( color = richness_shared, size = richness_shared)) +
  geom_smooth(method = "lm", color = "red", alpha = 0.2) +
  labs(y = "Green Roof richness", x = "Ground-level richness") +
  scale_size_continuous(range = c(2, 6), name = "Shared species richness") +
  scale_color_gradient(high = "black", low = "grey90", name = "Shared species richness") +
  theme_pubr(base_size = 10) +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend("Shared species richness"),
    size = guide_legend("Shared species richness")
  )
