################################################################################
##################################### fig5 #####################################
################################################################################
library(tidyr)
library(ggplot2)
library(ggpubr)

# Function to generate density distribution plots
get_distribution_plot <- function(df, 
                                  covariates, 
                                  colors = c("#C85200FF", "#006BA4FF", "#595959FF", "#FFBC79FF", "#A2C8ECFF"),
                                  labels = c("Observed", "Semi-randomization", "Random pairing"),
                                  x_axis, 
                                  y_axis,
                                  show_ticks = T,
                                  xlim) {
  
  # Convert the selected columns to long format for ggplot compatibility
  df.long <- pivot_longer(df[, colnames(df) %in% covariates], 
                          cols = covariates)
  
  # Compute density distributions for each covariate
  density.df <- data.frame()
  for (cov in covariates) {
    d <- data.frame(x = density(df[, cov])$x, y = density(df[, cov])$y)
    d$y[abs(cumsum(d$y * mean(diff(d$x))) - 0.5) > 0.45] <- 0     # Find where 95% of the data is
    density.df <- rbind(density.df, data.frame(d, type = cov))
  }
  
  if (show_ticks == TRUE) {
    axis_text_x_theme <- element_text(angle = 90, size = 16, hjust = 1)
  } else {
    axis_text_x_theme <- element_blank()
  }
  
  plot <- ggplot(df.long, aes(x = value, fill = name, color = name)) +
    geom_area(data = density.df, aes(x = x, y = y, fill = type, color = type), alpha = 0.5, position = "dodge", inherit.aes = F) +
    geom_density(alpha = 0, linetype = "solid") +
    scale_fill_manual(name = "Model",
                      values = colors,
                      labels = labels) +
    scale_color_manual(name = "Model",
                       values = colors,
                       labels = labels) +
    labs(x = x_axis, y = y_axis) + 
    theme_pubr(base_size = 18) + 
    coord_flip(xlim = xlim) + 
    guides(fill = guide_legend(title = "Model", override.aes = list(alpha = 1))) + 
    theme(
      plot.background = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = axis_text_x_theme,  
      plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
      plot.title = element_blank(), 
      plot.caption = element_blank()  
    )
  
  return(plot)  
}

# Function to determine line type based on statistical significance
add_significance_linetype <- function(df, x_var, y_var) {
  
  # Use a Poisson model if the response variable is species richness, otherwise use linear regression
  if (grepl(y_var, "sp_shared") == T) { 
    model <- glm(df[, x_var] ~ df[, y_var], family = "poisson", data = df)
  } else{
    model <- lm(df[, x_var] ~ df[, y_var], data = df)
  }
  
  # Extract the p-value for the slope
  p_value <- summary(model)$coefficients[2, 4] 
  
  # Use a solid line if significant (p < 0.05), otherwise use dashed
  if (p_value < 0.05) {
    return("solid")
  } else {
    return("dashed")
  }
}

# Import data
spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
shared.df <- read.csv("skyline_sharedmat.csv", check.names = F)
taxo.df <- read.csv("skyline_taxo.csv")
meta <- c("site_no", "type", "sp_richness")

# Add species richness values for roof and ground habitats
shared.df$sp_roof <- spmat.df[spmat.df$type == "Roof", ]$sp_richness
shared.df$sp_ground <- spmat.df[spmat.df$type == "Ground", ]$sp_richness

# Generate shared species richness density plot
sp_shared.plot <- get_distribution_plot(shared.df, covariates = c("sp_shared", 
                                                                  "sp_shared.pair",
                                                                  "sp_shared.rand"), x_axis = "", y_axis = "", xlim = c(0, 22), show_ticks = T, 
                                        labels = c("Observed", "Random pairing", "Semi-randomization"),
                                        colors = c("#C85200FF", "#006BA4FF", "#595959FF")) 

# Generate Jaccard dissimilarity density plot
jaccard.plot <- get_distribution_plot(shared.df, covariates = c("jaccard", 
                                                                "jaccard.pair",
                                                                "jaccard.rand"), x_axis = "", y_axis = "Density", xlim = c(0.83, 1), show_ticks = T,
                                      labels = c("Observed", "Random pairing", "Semi-randomization"),
                                      colors = c("#C85200FF", "#006BA4FF", "#595959FF", "#FFBC79FF")) 

# Scatter plot: shared species richness vs. height with regression line
sp_shared.plot1 <- ggplot(shared.df, aes(x = height)) + 
  geom_point(aes(y = sp_shared), color = "#C85200FF") +
  geom_smooth(aes(y = sp_shared), linetype = add_significance_linetype(shared.df, "height", "sp_shared"), color = "#C85200FF", method = "lm") +
  geom_point(aes(y = sp_shared.pair), color = "#006BA4FF") +
  geom_smooth(aes(y = sp_shared.pair), linetype = add_significance_linetype(shared.df, "height", "sp_shared.pair"), color = "#006BA4FF", method = "lm") +
  geom_point(aes(y = sp_shared.rand), color = "#595959FF") +
  geom_smooth(aes(y = sp_shared.rand), linetype = add_significance_linetype(shared.df, "height", "sp_shared.rand"), color = "#595959FF", method = "lm") +
  coord_cartesian(ylim = c(0, 22)) + 
  labs(x = "", y = "Shared species richness") + 
  theme_pubr(base_size = 18) + 
  theme(
    plot.background = element_blank(),
    axis.text.y = element_text(angle = 90, size = 16),
    axis.text.x = element_text(size = 16),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank(),  
    plot.caption = element_blank()   
  )

# Scatter plot: shared species richness vs. log-transformed paired distance
sp_shared.plot2 <- ggplot(shared.df, aes(x = log(paired_distance))) + 
  geom_point(aes(y = sp_shared), color = "#C85200FF") +
  geom_smooth(aes(y = sp_shared), linetype = add_significance_linetype(shared.df, "paired_distance", "sp_shared"), color = "#C85200FF", method = "lm") +
  geom_point(aes(y = sp_shared.pair), color = "#006BA4FF") +
  geom_smooth(aes(y = sp_shared.pair), linetype = add_significance_linetype(shared.df, "paired_distance", "sp_shared.pair"), color = "#006BA4FF", method = "lm") +
  geom_point(aes(y = sp_shared.rand), color = "#595959FF") +
  geom_smooth(aes(y = sp_shared.rand), linetype = add_significance_linetype(shared.df, "paired_distance", "sp_shared.rand"), color = "#595959FF", method = "lm") +
  coord_cartesian(ylim = c(0, 22)) + 
  labs(x = "", y = "") + 
  theme_pubr(base_size = 18) + 
  theme(
    plot.background = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank(),  
    plot.caption = element_blank()  
  )

# Scatter plot: shared species richness vs. species richness on roofs
sp_shared.plot3 <- ggplot(shared.df, aes(x = sp_roof)) + 
  geom_point(aes(y = sp_shared), color = "#C85200FF") +
  geom_smooth(aes(y = sp_shared), linetype = add_significance_linetype(shared.df, "sp_roof", "sp_shared"), color = "#C85200FF", method = "lm") +
  geom_point(aes(y = sp_shared.pair), color = "#006BA4FF") +
  geom_smooth(aes(y = sp_shared.pair), linetype = add_significance_linetype(shared.df, "sp_roof", "sp_shared.pair"), color = "#006BA4FF", method = "lm") +
  geom_point(aes(y = sp_shared.rand), color = "#595959FF") +
  geom_smooth(aes(y = sp_shared.rand), linetype = add_significance_linetype(shared.df, "sp_roof", "sp_shared.rand"), color = "#595959FF", method = "lm") +
  coord_cartesian(ylim = c(0, 22)) + 
  labs(x = "", y = "") + 
  theme_pubr(base_size = 18) + 
  theme(
    plot.background = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank(), 
    plot.caption = element_blank()
  )

# Scatter plot: shared species richness vs. species richness on the ground
sp_shared.plot4 <- ggplot(shared.df, aes(x = sp_ground)) + 
  geom_point(aes(y = sp_shared), color = "#C85200FF") +
  geom_smooth(aes(y = sp_shared), linetype = add_significance_linetype(shared.df, "sp_ground", "sp_shared"), color = "#C85200FF", method = "lm") +
  geom_point(aes(y = sp_shared.pair), color = "#006BA4FF") +
  geom_smooth(aes(y = sp_shared.pair), linetype = add_significance_linetype(shared.df, "sp_ground", "sp_shared.pair"), color = "#006BA4FF", method = "lm") +
  geom_point(aes(y = sp_shared.rand), color = "#595959FF") +
  geom_smooth(aes(y = sp_shared.rand), linetype = add_significance_linetype(shared.df, "sp_ground", "sp_shared.rand"), color = "#595959FF", method = "lm") +
  coord_cartesian(ylim = c(0, 22)) + 
  labs(x = "", y = "") + 
  theme_pubr(base_size = 18) + 
  theme(
    plot.background = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank(),  
    plot.caption = element_blank()
  )

# Repeat the scatter plots but for Jaccard dissimilarity
jaccard.plot1 <- ggplot(shared.df, aes(x = height)) + 
  geom_point(aes(y = jaccard), color = "#C85200FF") +
  geom_smooth(aes(y = jaccard), linetype = add_significance_linetype(shared.df, "height", "jaccard"), color = "#C85200FF", method = "lm") +
  geom_point(aes(y = jaccard.pair), color = "#006BA4FF") +
  geom_smooth(aes(y = jaccard.pair), linetype = add_significance_linetype(shared.df, "height", "jaccard.pair"), color = "#006BA4FF", method = "lm") +
  geom_point(aes(y = jaccard.rand), color = "#595959FF") +
  geom_smooth(aes(y = jaccard.rand), linetype = add_significance_linetype(shared.df, "height", "jaccard.rand"), color = "#595959FF", method = "lm") +
  coord_cartesian(ylim = c(0.83, 1)) + 
  labs(x = "Height", y = "Jaccard dissimilarity") + 
  theme_pubr(base_size = 18) + 
  theme(
    plot.background = element_blank(),
    axis.text.y = element_text(angle = 90, size = 16),
    axis.text.x = element_text(angle = 90, size = 16, hjust = 1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank(),
    plot.caption = element_blank() 
  )

jaccard.plot2 <- ggplot(shared.df, aes(x = log(paired_distance))) + 
  geom_point(aes(y = jaccard), color = "#C85200FF") +
  geom_smooth(aes(y = jaccard), linetype = add_significance_linetype(shared.df, "paired_distance", "jaccard"), color = "#C85200FF", method = "lm") +
  geom_point(aes(y = jaccard.pair), color = "#006BA4FF") +
  geom_smooth(aes(y = jaccard.pair), linetype = add_significance_linetype(shared.df, "paired_distance", "jaccard.pair"), color = "#006BA4FF", method = "lm") +
  geom_point(aes(y = jaccard.rand), color = "#595959FF") +
  geom_smooth(aes(y = jaccard.rand), linetype = add_significance_linetype(shared.df, "paired_distance", "jaccard.rand"), color = "#595959FF", method = "lm") +
  coord_cartesian(ylim = c(0.83, 1)) + 
  labs(x = "Distance (log)", y = "") + 
  theme_pubr(base_size = 18) + 
  theme(
    plot.background = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 16, hjust = 1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank(), 
    plot.caption = element_blank() 
  )

jaccard.plot3 <- ggplot(shared.df, aes(x = sp_roof)) + 
  geom_point(aes(y = jaccard), color = "#C85200FF") +
  geom_smooth(aes(y = jaccard), linetype = add_significance_linetype(shared.df, "sp_roof", "jaccard"), color = "#C85200FF", method = "lm") +
  geom_point(aes(y = jaccard.pair), color = "#006BA4FF") +
  geom_smooth(aes(y = jaccard.pair), linetype = add_significance_linetype(shared.df, "sp_roof", "jaccard.pair"), color = "#006BA4FF", method = "lm") +
  geom_point(aes(y = jaccard.rand), color = "#595959FF") +
  geom_smooth(aes(y = jaccard.rand), linetype = add_significance_linetype(shared.df, "sp_roof", "jaccard.rand"), color = "#595959FF", method = "lm") +
  coord_cartesian(ylim = c(0.83, 1)) + 
  labs(x = "Green roof\nrichness", y = "") + 
  theme_pubr(base_size = 18) + 
  theme(
    plot.background = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 16, hjust = 1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank(), 
    plot.caption = element_blank() 
  )

jaccard.plot4 <- ggplot(shared.df, aes(x = sp_ground)) + 
  geom_point(aes(y = jaccard), color = "#C85200FF") +
  geom_smooth(aes(y = jaccard), linetype = add_significance_linetype(shared.df, "sp_ground", "jaccard"), color = "#C85200FF", method = "lm") +
  geom_point(aes(y = jaccard.pair), color = "#006BA4FF") +
  geom_smooth(aes(y = jaccard.pair), linetype = add_significance_linetype(shared.df, "sp_ground", "jaccard.pair"), color = "#006BA4FF", method = "lm") +
  geom_point(aes(y = jaccard.rand), color = "#595959FF") +
  geom_smooth(aes(y = jaccard.rand), linetype = add_significance_linetype(shared.df, "sp_ground", "jaccard.rand"), color = "#595959FF", method = "lm") +
  coord_cartesian(ylim = c(0.83, 1)) + 
  labs(x = "Ground-level\nrichness", y = "") + 
  theme_pubr(base_size = 18) + 
  theme(
    plot.background = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 16, hjust = 1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_blank(),
    plot.caption = element_blank() 
  )

ggarrange(sp_shared.plot1, sp_shared.plot2, sp_shared.plot3, sp_shared.plot4, sp_shared.plot, 
          jaccard.plot1, jaccard.plot2, jaccard.plot3, jaccard.plot4, jaccard.plot,
          common.legend = T,
          align = "h",
          widths = c(1, 1, 1, 1, 0.8),
          ncol = 5, nrow = 2)
